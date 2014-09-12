do.pop.predict.balance <- function(inp, outdir, nr.traj, ages, pred=NULL, keep.vital.events=FALSE, function.inputs=NULL, 
									rebalance=TRUE, start.time.index=1, 
									verbose=FALSE, parallel=FALSE, nr.nodes=NULL, .countries=NULL, .migration.pars=NULL, ...) {
	countries.idx <- if(is.null(.countries)) which(UNlocations$location_type==4) else which(UNlocations$country_code %in% .countries)
	country.codes <- UNlocations$country_code[countries.idx]
	ncountries <- length(country.codes)
	nr_project <- length(inp$proj.years)
	nages <- length(ages)
	mx.ages <- c(0,1,ages[2:nages])
	ages21 <- ages[ages<=100]
	nest <- length(inp$estim.years)
	if(!file.exists(outdir)) 
		dir.create(outdir, recursive=TRUE)
	present.and.proj.years <- c(inp$estim.years[nest], inp$proj.years)
	
	status.for.gui <- paste('out of', nr_project, 'time periods.')
	gui.options <- list()
	inp.to.save <- list()
	# remove big or redundant items from inputs to be saved
	for(item in ls(inp)[!grepl('^migMpred$|^migFpred$|^TFRpred$|^e0Fpred$|^e0Mpred$|^estim.years$|^proj.years$|^wpp.years$', ls(inp))]) 
		inp.to.save[[item]] <- get(item, inp)
	npred <- nr_project
	
	mig.rate.prev <- NULL
	if(!is.null(.migration.pars)) {
		inp[['migration.parameters']] <- .migration.pars$posterior
		mig.rate.prev <- .migration.pars$ini.rates#/5.
	}
	outdir.tmp <- file.path(outdir, '_tmp_')
	if(file.exists(outdir.tmp) && start.time.index==1) unlink(outdir.tmp, recursive=TRUE)
	if(start.time.index==1) dir.create(outdir.tmp)
	
	if(start.time.index > 1) { # reload initial population
		env.tmp <- new.env()
		load(file.path(outdir.tmp, paste0('pop_time_', start.time.index-1, '.rda')), envir=env.tmp)
		popM.prev <- env.tmp$totpm
		popF.prev <- env.tmp$totpf
		popM.hch.prev <- env.tmp$totpm.hch
		popF.hch.prev <- env.tmp$totpf.hch
		mig.rate.prev <- env.tmp$mig.rate
	}
	UNnames <- UNlocations[countries.idx,'name']
	if(parallel) {
		if(is.null(nr.nodes)) nr.nodes <- getOption("cl.cores", detectCores())
		nr.nodes <- min(nr.nodes, ncountries)
	}
	countries.input <- new.env()
	# Extract the country-specific stuff from the inputs
	if(verbose) cat('\nLoading inputs for ', ncountries, ' countries ')
	if(parallel) {
		if(verbose) cat('(in parallel on ', nr.nodes, ' nodes).')
		cl <- create.pop.cluster(nr.nodes, ...)
		clusterExport(cl, c("inp", "nr.traj", "country.codes", "UNnames"), envir=environment())
		inpc.list <- parLapplyLB(cl, 1:ncountries, 
						function(cidx) get.country.inputs(country.codes[cidx], inp, nr.traj, UNnames[cidx]))
		stopCluster(cl)
		for(cidx in 1:ncountries) {	
			if(is.null(inpc.list[[cidx]]) || length(inpc.list[[cidx]]$POPm0)==0) next
			countries.input[[as.character(country.codes[cidx])]] <- inpc.list[[cidx]]
		}
	} else { # load inputs sequentially
		if(verbose) cat('(sequentially).')
		for(cidx in 1:ncountries) {
			inpc <- get.country.inputs(country.codes[cidx], inp, nr.traj, UNnames[cidx])
			if(is.null(inpc) || length(inpc$POPm0)==0) next
			countries.input[[as.character(country.codes[cidx])]] <- inpc
		}
	}
	nr.traj <- min(nr.traj, sapply(ls(countries.input), function(ccode) return(ncol(countries.input[[ccode]]$TFRpred))))
	npred <- min(npred, sapply(ls(countries.input), function(ccode) return(nrow(countries.input[[ccode]]$TFRpred))))
	npredplus1 <- npred+1
	npasfr <- nrow(countries.input[[as.character(country.codes[1])]]$PASFR)
	nvariants <- nrow(countries.input[[as.character(country.codes[1])]]$TFRhalfchild)
	if(verbose) cat('\nProcessing ', nr.traj, ' trajectories for each country ')
	if(length(countries.input) < ncountries) {
		ncountries <- length(countries.input)
		country.codes <- as.integer(ls(countries.input))
		countries.idx <- sapply(country.codes, function(x) which(UNlocations[,'country_code'] == x))
		UNnames <- UNlocations[countries.idx,'name']
		if(parallel) nr.nodes <- min(nr.nodes, ncountries)
	}
	if(!is.null(mig.rate.prev)) {
		mig.rate.prev <- mig.rate.prev[as.character(country.codes)]
		mig.rate.prev <- matrix(mig.rate.prev, nrow=length(mig.rate.prev), ncol=nr.traj)
	}
	kannisto <- list()
	observed <- new.env()
	for(country in ls(countries.input)) {
		kannisto[[country]] <- runKannisto(nest, countries.input[[country]])
		if(keep.vital.events) 
			observed[[country]] <- compute.observedVE(countries.input[[country]], inp$pop.matrix, 
										countries.input[[country]]$MIGtype, kannisto[[country]], 
										as.integer(country), inp$estim.years)
	}
	res.env <- new.env()
	with(res.env, {
		totp <- mig.rate <- matrix(0, nrow=ncountries, ncol=nr.traj, dimnames=list(country.codes, NULL))
		totpm <- totpf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
		migrationm <- migrationf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
		totp.hch <- matrix(0, nrow=ncountries, ncol=nvariants, dimnames=list(country.codes, NULL))
		totpm.hch <- totpf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
		migrationm.hch <- migrationf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
		if(keep.vital.events) {
			btm <- btf <- asfert <- array(0, dim=c(7, ncountries, nr.traj), dimnames=list(NULL, country.codes, NULL))
			deathsm <- deathsf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
			btm.hch <- btf.hch <- asfert.hch <- array(0, dim=c(7, ncountries, nvariants), dimnames=list(NULL, country.codes, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
			mxm <- mxf <- array(0, dim=c(28, ncountries, nr.traj), dimnames=list(mx.ages, country.codes, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(28, ncountries, nvariants), dimnames=list(mx.ages, country.codes, NULL))
			migntraj <- countries.input[[as.character(country.codes[1])]][['mig.nr.traj']]
			migm <- array(0, dim=c(27, ncountries, migntraj), dimnames=list(ages, country.codes, NULL))
			migf <- array(0, dim=c(27, ncountries, migntraj), dimnames=list(ages, country.codes, NULL))
		}
	})
	debug <- FALSE
	if(parallel) {
		if(verbose) cat(' (in parallel on ', nr.nodes, ' nodes).')
		cl <- create.pop.cluster(nr.nodes, ...)
		clusterExport(cl, c("nr.traj", "country.codes",  "UNnames", "countries.input", "kannisto", "npasfr", 
								"ages", "nvariants", 
								"keep.vital.events", "verbose"), envir=environment())
	} else if(verbose) cat(' (sequentially).')
	for(time in start.time.index:npred) {
		unblock.gtk.if.needed(paste('finished', time, status.for.gui), gui.options)
		if(verbose)
			cat('\nProcessing time period ', time)
		.ini.pop.res.env(res.env, keep.vital.events)	
		if(parallel) {
			clusterExport(cl, c("time", "mig.rate.prev"), envir=environment())
			if(time > 1) clusterExport(cl, c("popM.prev", "popF.prev", "popM.hch.prev", "popF.hch.prev"), envir=environment())
			res.list <- parLapplyLB(cl, 1:ncountries, 
							function(cidx) do.pop.predict.one.country.no.migration(time, 
												UNnames[cidx], countries.input[[as.character(country.codes[cidx])]], 
												kannisto[[as.character(country.codes[cidx])]], nr.traj, npasfr, 
												if(time>1)popM.prev[,cidx, ,drop=FALSE] else NULL, 
												if(time>1)popF.prev[,cidx, ,drop=FALSE] else NULL, 
												if(time>1)popM.hch.prev[,cidx, ,drop=FALSE] else NULL, 
												if(time>1)popF.hch.prev[,cidx, ,drop=FALSE] else NULL,
												#list(M=migpred$M[cidx,,, drop=FALSE], F=migpred$F[cidx,,, drop=FALSE]),
												ages, nvariants, 
												keep.vital.events=keep.vital.events, verbose=verbose))
			
		} else { # process sequentially
			res.list <- list()
			for(cidx in 1:ncountries) {
				res.list[[cidx]] <- do.pop.predict.one.country.no.migration(time, 
												UNnames[cidx], countries.input[[as.character(country.codes[cidx])]], 
												kannisto[[as.character(country.codes[cidx])]], nr.traj, npasfr, 
												if(time>1)popM.prev[,cidx, ,drop=FALSE] else NULL, 
												if(time>1)popF.prev[,cidx, ,drop=FALSE] else NULL, 
												if(time>1)popM.hch.prev[,cidx, ,drop=FALSE] else NULL, 
												if(time>1)popF.hch.prev[,cidx, ,drop=FALSE] else NULL, 
												#list(M=migpred$M[cidx,,, drop=FALSE], F=migpred$F[cidx,,, drop=FALSE]),
												ages, nvariants, 
												keep.vital.events=keep.vital.events, verbose=verbose)
			}
		}
		for(cidx in 1:ncountries) { # collect results
			for(par in c('totp', 'totp.hch'))
				res.env[[par]][cidx,] <- res.list[[cidx]][[par]]
			for(par in c('totpm', 'totpf', 'totpm.hch', 'totpf.hch'))#'migrationm', 'migrationf', 'migrationm.hch', 'migrationf.hch'
				res.env[[par]][,cidx,] <- res.list[[cidx]][[par]]
			if(keep.vital.events) {
				for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'mxm', 'mxf', 
							'btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'mxm.hch','mxf.hch')) # 'migm', 'migf',
					res.env[[par]][,cidx,] <- res.list[[cidx]][[par]]
			}
		}
		migpred <- get.balanced.migration(time, country.codes, countries.input, nr.traj, rebalance, mig.rate.prev,
							if(time>1)popM.prev[,, ,drop=FALSE] else NULL, 
							if(time>1)popF.prev[,, ,drop=FALSE] else NULL, ages, res.env, verbose=verbose)
		lage <- dim(migpred$M)[1]
		res.env$migrationm[1:lage,,] <- migpred$M
		res.env$migrationf[1:lage,,] <- migpred$F
		# TODO: what should be migrationm.hch?
		res.env$migrationm.hch[1:lage,,] <- apply(migpred$M, c(1,2), mean)
		res.env$migrationf.hch[1:lage,,] <- apply(migpred$F, c(1,2), mean)
		res.env$mig.rate <- migpred$rate
		# New population counts
		res.env$totpm <- res.env$totpm + res.env$migrationm
		res.env$totpf <- res.env$totpf + res.env$migrationf
		res.env$totp <- apply(res.env$totpm + res.env$totpf, c(2,3), sum)
		res.env$totpm.hch <- res.env$totpm.hch + res.env$migrationm.hch
		res.env$totpf.hch <- res.env$totpf.hch + res.env$migrationf.hch
		res.env$totp.hch <- apply(res.env$totpm.hch + res.env$totpf.hch, c(2,3), sum)
		
		if(keep.vital.events) {
			# TODO: Are these really the same? Do we need both?
			res.env$migm <- res.env$migrationm
			res.env$migf <- res.env$migrationf
		}
		# print accummulated warnings
		all.warns <- unique(unlist(migpred$warns))
		print.warn <- function(key, value, warn) if(warn %in% value) return(key) else return(c())
		for(warn in all.warns) {
			cntrs <- mapply(print.warn, names(migpred$warns), migpred$warns, warn=warn)
			warning(warn, ': ', paste(unlist(cntrs), collapse=', '), immediate.=TRUE)
		}
		if (any(res.env$totpm < -1e-4) || any(res.env$totpf < -1e-4))
			warning('Final population negative for some countries and age groups.', immediate.=TRUE)
		popM.prev <- res.env$totpm
		popF.prev <- res.env$totpf
		popM.hch.prev <- res.env$totpm.hch
		popF.hch.prev <- res.env$totpf.hch
		mig.rate.prev <- res.env$mig.rate
		
		with(res.env, {
			save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch, mig.rate,
				 file = file.path(outdir.tmp, paste0('pop_time_', time, '.rda')))
			if(keep.vital.events) 
				save(btm, btf, deathsm, deathsf, asfert, mxm, mxf, migm, migf,
						btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, 
						mxm.hch, mxf.hch, 
						file=file.path(outdir.tmp, paste0('vital_events_time_', time, '.rda')))
		})
		gc()
	} # end time
	if(parallel) stopCluster(cl)
	if(verbose) cat('\nRe-formatting data ')
	quant.env <- restructure.pop.data.and.compute.quantiles(outdir.tmp, outdir, npred, countries.input, observed, kannisto, 
					present.and.proj.years, keep.vital.events, parallel=parallel, nr.nodes=nr.nodes, verbose=verbose)
	if(verbose) cat(' done.\n')
	unlink(outdir.tmp, recursive=TRUE)
	#save meta file
	country.rows <- UNlocations[countries.idx,c('country_code', 'name')]
	colnames(country.rows) <- c('code', 'name')
	bayesPop.prediction <- structure(list(
							nr.traj = nr.traj,	
							quantiles = quant.env$PIs_cqp,
               				traj.mean.sd = quant.env$mean_sd,
               				quantilesM = quant.env$quantM, 
               				traj.mean.sdM = quant.env$mean_sdM,
               				quantilesF = quant.env$quantF, 
               				traj.mean.sdF = quant.env$mean_sdF,
               				quantilesMage = quant.env$quantMage, 
               				quantilesFage = quant.env$quantFage, 
               				quantilesPropMage = quant.env$quantPropMage, 
               				quantilesPropFage = quant.env$quantPropFage,
               				estim.years=inp$estim.years, 
               				proj.years=present.and.proj.years, # includes present period (middle of periods)
               				proj.years.pop=present.and.proj.years+2, # end of periods
               				wpp.year = inp$wpp.year,
			   				inputs = inp.to.save, # save as list because environment takes much more space
			   				function.inputs=function.inputs,
			   				countries=country.rows,
			   				ages=ages), class='bayesPop.prediction')

	prediction.file <- file.path(outdir, 'prediction.rda')
	save(bayesPop.prediction, file=prediction.file)
	cat('\nPrediction stored into', outdir, '\n')
	return(bayesPop.prediction)
}

get.balanced.migration <- function(time, country.codes, inputs, nr.traj, rebalance, mig.rate.prev,
												popM.prev, popF.prev, ages, env, verbose=FALSE) {
	nr.countries <- length(country.codes)
	e <- new.env()
	labor.codes <- labor.countries()
	ages <- ages[1:21]
	lages <- length(ages)
	e$migrm <- e$migrf <- matrix(NA, ncol=nr.countries, nrow=lages, dimnames=list(ages, country.codes))
	e$migrm.labor <- e$migrf.labor <- matrix(0, ncol=nr.countries, nrow=lages, dimnames=list(ages, country.codes))
	migrationm <- migrationf <- array(NA, c(lages, nr.countries, nr.traj), dimnames=list(ages, country.codes, NULL))
	mig.rate <- matrix(NA, nrow=nr.countries, ncol=nr.traj, dimnames=list(country.codes,NULL))
	
	pop <- rep(NA, nr.countries)
	names(pop) <- country.codes
	data(land_area_wpp2012)	
	warns <- list()
	for(itraj in 1:nr.traj){
		for(cidx in 1:nr.countries) {
			inpc <- inputs[[as.character(country.codes[cidx])]]
			#pop.ini.by.group <- if(time > 1) list(M=popM.prev[,cidx,itraj], F=popF.prev[,cidx,itraj])
			#					else list(M=inpc$POPm0, F=inpc$POPf0)
			pop.ini.by.group <- list(M=env$totpm[,cidx,itraj], F=env$totpf[,cidx,itraj])
			pop.ini <- sum(pop.ini.by.group$M + pop.ini.by.group$F)
			migpred <- .get.migration.one.trajectory(inpc, itraj, time, 
								if(!is.null(mig.rate.prev)) mig.rate.prev[cidx, itraj] else NULL, pop.ini, 
								pop.group=pop.ini.by.group, country.code=country.codes[cidx])
			e$migrm[,cidx] <- migpred$M
			e$migrf[,cidx] <- migpred$F
			if(!is.null(migpred$laborM)) {
				e$migrm.labor[,cidx] <- migpred$laborM
				e$migrf.labor[,cidx] <- migpred$laborF
			}
			mig.rate[cidx, itraj] <- migpred$rate
			pop[cidx] <- pop.ini
			warns[[as.character(country.codes[cidx])]] <- unique(c(warns[[as.character(country.codes[cidx])]], migpred$warns))
		}
		e$pop.by.age <- list(m=env$totpm[,,itraj], f=env$totpf[,,itraj])
		if(rebalance) {
			rebalance.migration2groups(e, pop, itraj)			
			negatives <- country.codes[unique(e$negatives)]
			for(country in negatives)
				warns[[as.character(country)]] <- unique(c(warns[[as.character(country)]], 'population negative while balancing'))
		}
		migrationm[,,itraj] <- e$migrm + e$migrm.labor
		migrationf[,,itraj] <- e$migrf + e$migrf.labor
	}
	return(list(M=migrationm, F=migrationf, rate=mig.rate, warns=warns))
}

do.pop.predict.one.country.no.migration <- function(time, country.name, inpc, kannisto, nr.traj, npasfr, 
												popM.prev, popF.prev, popM.hch.prev, popF.hch.prev, 
												ages, nvariants, keep.vital.events, verbose=FALSE) {
	pop.ini <- list(M=inpc$POPm0, F=inpc$POPf0)
	res.env <- new.env()
	mx.ages <- c(0,1,ages[2:length(ages)])
	with(res.env, {
		totp <-  rep(NA, nr.traj) #mig.rate <-
		totpm <- totpf <- matrix(NA, nrow=27, ncol=nr.traj, dimnames=list(ages, NULL)) # migrationm <- migrationf <- 
		totp.hch <- rep(NA, nvariants)
		totpm.hch <- totpf.hch <- matrix(NA, nrow=27, ncol=nvariants, dimnames=list(ages, NULL)) #migrationm.hch <- migrationf.hch <- 
		if(keep.vital.events) {
			btm <- btf <- asfert <- matrix(0, nrow=7, ncol=nr.traj)
			deathsm <- deathsf <- matrix(0, nrow=27, ncol=nr.traj, dimnames=list(ages, NULL))
			btm.hch <- btf.hch <- asfert.hch <- matrix(0, nrow=7, ncol=nvariants)
			deathsm.hch <- deathsf.hch <- matrix(0, nrow=27, ncol=nvariants, dimnames=list(ages, NULL))
			mxm <- mxf <- matrix(0, nrow=28, ncol=nr.traj, dimnames=list(mx.ages, NULL))
			mxm.hch <- mxf.hch <- matrix(0, nrow=28, ncol=nvariants, dimnames=list(mx.ages, NULL))
			#migntraj <- inpc$mig.nr.traj
			#migm <- matrix(0, nrow=27, ncol=migntraj, dimnames=list(ages, NULL))
			#migf <- matrix(0, nrow=27, ncol=migntraj, dimnames=list(ages, NULL))
		}
	})
	# if migration model used then migration at the end of each interval for all countries
	#mig.type <- if(is.null(inpc$migration.parameters)) inpc$MIGtype else 9 
	for(itraj in 1:nr.traj) {
		asfr <- inpc$PASFR[,time,drop=FALSE]/100.
		for(i in 1:npasfr) asfr[i,] <- inpc$TFRpred[time,itraj] * asfr[i,]
		LTres <- modifiedLC(1, kannisto, inpc$e0Mpred[time,itraj], 
									inpc$e0Fpred[time,itraj], verbose=verbose)		
		if(time > 1) { # reset initial population to the one at the previous time step
			pop.ini$M <- popM.prev[,1,itraj]
			pop.ini$F <- popF.prev[,1,itraj]
		}
		#this.migpred <- list(M=migpred$M[1,,itraj], F=migpred$F[1,,itraj])
		#popres <- StoPopProj(1, pop.ini, LTres, asfr, inpc$SRB[time], this.migpred, mig.type, country.name=country.name,
		#							keep.vital.events=keep.vital.events)
		popres <- PopProjNoMigr(1, pop.ini, LTres, asfr, inpc$SRB[time], country.name=country.name,
									keep.vital.events=keep.vital.events)
		with(res.env, {
			totp[itraj] <- popres$totpop[2]
			totpm[,itraj] <- popres$mpop[,2]
			totpf[,itraj] <- popres$fpop[,2]
			#migrationm[,itraj] <- popres$mmigr
			#migrationf[,itraj] <- popres$fmigr
			if(keep.vital.events) {
				btm[,itraj] <- popres$mbt
				btf[,itraj] <- popres$fbt
				deathsm[,itraj] <- popres$mdeaths
				deathsf[,itraj] <- popres$fdeaths
				asfert[,itraj] <- asfr
				mxm[,itraj] <- LTres$mx[[1]]
				mxf[,itraj] <- LTres$mx[[2]]
				#migtraj <- min(itraj, ncol(migm))
				#migm[,migtraj] <- popres$mmigr
				#migtraj <- min(itraj, ncol(migf))
				#migf[,migtraj] <- popres$fmigr 
			}
		})
	}
	pop.ini <- list(M=inpc$POPm0, F=inpc$POPf0)
	for (variant in 1:nvariants) { # compute the two half child variants
		asfr <- inpc$PASFR[,time,drop=FALSE]/100.
		for(i in 1:npasfr) asfr[i,] <- inpc$TFRhalfchild[variant,time] * asfr[i,]
		LTres <- modifiedLC(1, kannisto, inpc$e0Mmedian[time], 
								inpc$e0Fmedian[time], verbose=verbose, debug=debug)
		if(time > 1) {
			pop.ini$M <- popM.hch.prev[,1, variant]
			pop.ini$F <- popF.hch.prev[,1, variant]
		}
		#migpred.hch <- list(M=apply(res.env$migrationm, 1, median), F=apply(res.env$migrationf, 1, median))
		#popres <- StoPopProj(1, pop.ini, LTres, asfr, inpc$SRB[time], migpred.hch, mig.type, 
		#					country.name=country.name, keep.vital.events=keep.vital.events)
		popres <- PopProjNoMigr(1, pop.ini, LTres, asfr, inpc$SRB[time], 
							country.name=country.name, keep.vital.events=keep.vital.events)
		with(res.env, {
			totp.hch[variant] <- popres$totpop[2]
			totpm.hch[,variant] <- popres$mpop[,2]
			totpf.hch[,variant] <- popres$fpop[,2]
			#migrationm.hch[,variant] <- popres$mmigr
			#migrationf.hch[,variant] <- popres$fmigr 
			if(keep.vital.events) {
				btm.hch[,variant] <- popres$mbt
				btf.hch[,variant] <- popres$fbt
				deathsm.hch[,variant] <- popres$mdeaths
				deathsf.hch[,variant] <- popres$fdeaths
				asfert.hch[,variant] <- asfr
				mxm.hch[,variant] <- LTres$mx[[1]]
				mxf.hch[,variant] <- LTres$mx[[2]]
			}
		})
	}
	res <- as.list(res.env)
	rm(list=ls(res.env), envir=res.env)
	return(res)
}

.get.migration.traj <- function(pred, par, country) {
		cidx <- pred$inputs[[par]][,'country_code'] == country 
		idx <- cidx & is.element(pred$inputs[[par]][,'year'], pred$inputs$proj.years)
		if(sum(idx) == 0) return(NULL)
		migdf <- pred$inputs[[par]][idx,-1]
		utrajs <- sort(unique(migdf$trajectory))
		ntrajs <- length(utrajs)
		migdf$age <- gsub("^\\s+|\\s+$", "", migdf$age) # trim leading and trailing whitespace
		lyears <- length(pred$inputs$proj.years)
		sorted.df <- data.frame(year=rep(pred$inputs$proj.years, each=ntrajs*21), trajectory=rep(rep(utrajs, each=21), times=lyears),
									age=c(paste(seq(0,95,by=5), seq(4,99,by=5), sep='-'), '100+'))
		# this is to get rows of the data frame in a particular order
		migdf <- merge(sorted.df, migdf, sort=FALSE)
		res <- array(migdf$value, dim=c(21, ntrajs, lyears))
		dimnames(res) <- list(1:21,  NULL, pred$inputs$proj.years)
		return(res)	
}

.get.migration.one.trajectory <- function(inpc, itraj=NULL, time=NULL, mig.rate.prev=NULL, pop=NULL, ...) {
	migpred <- list(M=NULL, F=NULL)
	if(!is.null(inpc$migration.parameters) && !is.null(itraj)) #TODO: what should it be for half child variants?
		return(sample.migration.trajectory.from.model(inpc, itraj, time, mig.rate.prev, pop, ...))
	for(sex in c('M', 'F')) {
		par <- paste0('mig', sex, 'pred')
		if(is.null(inpc[[par]]))
			migpred[[sex]] <- inpc[[paste0('MIG', tolower(sex))]]
		else {
			if(is.null(itraj)) {
				if(length(dim(inpc[[par]])) > 2) # has trajectory dimension, therefore need to take median
					migpred[[sex]] <- apply(inpc[[par]], c(1,3), 'median')
				else migpred[[sex]] <- inpc[[par]]
			} else migpred[[sex]] <- inpc[[par]][,itraj,]
		}
		migpred[[sex]] <- as.matrix(migpred[[sex]])
		if(!is.null(time)) migpred[[sex]] <- migpred[[sex]][,time]
	}
	return(migpred)
}


project.migration.one.country.one.step <- function(mu,phi,sigma,oldRate){
# Based on Jon Azose code 
#######################
#Project migration for a single country one time point into the future
#######################

  #Goal: Given a country's migration rate, project its migration RATE one step into the future
  #Inputs: mu, phi, and sigma parameters for an AR(1) model
  #        oldRate --- the country's migration rate at the "current" time point
  #        pop --- the projected total number of person years in this country in the next time period
  #Output: the projected total migration count over the next time period
  #
  #***NOTE: Units on migration rates and corresponding parameters are assumed to be the
  #         "natural" units, i.e. (total net # of migrants [not thousands] across 5 years)/(total person-years across 5 years)
  #         We can convert these to the scale reported by the UN (net *annual* migrants *per thousand*) by multiplying by 200.
  #
  #         The units on the pop inputs are similarly assumed to be total person-years across the five-year period [not thousands]
  #         The units on the output will projection will also be a total five-year count [not thousands]
  newRate <- 1
  while(newRate < -0.33 || newRate > 0.665)
  	newRate <- mu + phi*(oldRate-mu) + rnorm(n=1,mean=0,sd=sigma)
  return(newRate)
}

is.gcc <- function(country)
	return(country %in% c(634, 784, 414, 48, 512, 682)) # Qatar, UAE, Kuwait, Bahrain, Oman, SA

sample.migration.trajectory.from.model <- function(inpc, itraj=NULL, time=NULL, mig.rate.prev=NULL, pop=NULL, 
													pop.group=NULL, country.code=NULL) {													
	pars <- inpc$migration.parameters[itraj,]
	land.area <- NA
	if(country.code %in% land_area_wpp2012$country_code)
		land.area <- land_area_wpp2012[land_area_wpp2012$country_code==country.code,'land_area']
	i <- 1
	k <- 1
	warns <- c()
	while(TRUE) { 
		rate <- project.migration.one.country.one.step(pars$mu, pars$phi, pars$sigma, mig.rate.prev)
		if(is.na(rate)) stop('Migration rate is NA')
		#mig.count <- ((1+rate)^5 - 1) * pop.prev # instanteneous rate
		mig.count <- rate * pop
		if(is.gcc(country.code)) { # cap to historical maximum
			mig.count <- min(mig.count, max(colSums(inpc$observed$MIGm + inpc$observed$MIGf)))
		}
		if(!is.na(land.area) && (pop + mig.count)/land.area > 44) { # check density
			if(k<100) {
				k <- k+1
				next 
			}
			warns <- c(warns, 'migration truncated due to high density')
			mig.count <- 44 * land.area - pop
			rate <- mig.count/pop
			#warning('Density too high for ', country.code, ' (', (pop + mig.count)/land.area, ')', immediate.=TRUE)
		}
		schedMname <- 'M'
		schedFname <- 'F'
		if(rate < 0 && !is.null(inpc$migration.age.schedule[['Mnegative']])) {
			schedMname <- 'Mnegative'
			schedFname <- 'Fnegative'
		}
		msched <- inpc$migration.age.schedule[[schedMname]][,time]
		fsched <- inpc$migration.age.schedule[[schedFname]][,time]
		# age-specific migration counts		
		migM <- mig.count*msched
		migF <- mig.count*fsched

		if(all(pop.group$M[1:21] + migM >= 0) && all(pop.group$F[1:21] + migF >= -1e-4))  break # assure positive count
		i <- i+1
		if(i>100 && ((sum(pop.group$M) + sum(pop.group$F) + mig.count) > -1e-4) && (
				(sum(pop.group$M[1:21][msched>0]) + sum(pop.group$F[1:21][fsched>0]) + mig.count) > -1e-4)) { # adjust age schedules
			prev.isneg <- rep(FALSE, 42)
			j <- 1
			sample.new.rate <- FALSE
			while(any(c(pop.group$M[1:21] + migM, pop.group$F[1:21] + migF) < -1e-4)) { 
				isneg <- prev.isneg | (c(pop.group$M[1:21] + migM, pop.group$F[1:21] + migF) < -1e-4)
				shifts <- -c(migM + pop.group$M[1:21], migF + pop.group$F[1:21])
				shifts[!isneg] <- 0
				shifts[prev.isneg] <- 0
				if(sum(shifts)==0) {sample.new.rate <- TRUE; break}
				msched <- msched + shifts[1:21]/mig.count
				fsched <- fsched + shifts[22:length(shifts)]/mig.count
				tsched.neg <- sum(c(msched,fsched)[isneg])
				tsched.notneg <- sum(c(msched,fsched)[!isneg])
				msched[!isneg[1:21]] <- (1-tsched.neg) * msched[!isneg[1:21]]/tsched.notneg
				fsched[!isneg[22:length(shifts)]] <- (1-tsched.neg) * fsched[!isneg[22:length(shifts)]]/tsched.notneg
				migM <- mig.count * msched
				migF <- mig.count * fsched
				prev.isneg <- isneg
				j <- j+1
				#if(j > 21) stop('Age schedule cannot be adjusted.')
			}
			if(!sample.new.rate) {
				warns <- c(warns, 'migration age-schedule modified')
				break
			}
		}
		if(i > 10000) {
			warning('Unable to modify age schedule to get positive population for country ', country.code, immediate.=TRUE)
			break
		}
	}
	migM.labor <- migF.labor <- NULL
	if(country.code %in% labor.countries()) {
		prop <- prop.labor.migration.for.country(country.code)
		mig.count.labor <- mig.count * prop
		migM.labor <- mig.count.labor*msched
		migF.labor <- mig.count.labor*fsched
		mig.count.rest <- mig.count * (1-prop)
		migM <- mig.count.rest*msched
		migF <- mig.count.rest*fsched
	}
	return(list(M=migM, F=migF, rate=rate, laborM=migM.labor, laborF=migF.labor, warns=warns))
}

rebalance.population.by.migration <- function(e) { # not used
	nr.traj <- ncol(e$totp)
	for(itraj in 1:nr.traj) {
		countries.totals <- e$totp[,itraj]
		sumpop <- sum(countries.totals)
		sumadj <- difarray <- rep(0, nrow(e$totp))
		#if(itraj == 2) stop('')
		for(sex in c('m', 'f')) {
			par <- paste0('migration',sex)
			parpop <- paste0('totp',sex)
			for(age in dimnames(e$migrationm)[[1]]) {			
				dif <- sum(e[[par]][age,,itraj])
				difarray <- difarray + e[[par]][age,,itraj]
				dif.countries <- dif/sumpop * countries.totals
				e[[par]][age,,itraj] <- e[[par]][age,,itraj] - dif.countries
				e[[parpop]][age,,itraj] <- e[[parpop]][age,,itraj] + dif.countries
				sumadj <- sumadj + dif.countries
			}
		}
		if(sum(sumadj) > 1e+6) stop('')
		e$totp[,itraj] <- e$totp[,itraj] + sumadj		
		print(c(itraj, 'rebalancing migration', sum(sumadj)))
	}
}

rebalance.migration.with.mask <- function(e, pop, what='', pop.by.age, mask=NULL) {
	sumpop <- sum(pop)
	which.negative <- c()
	if(is.null(mask)) mask <- rep(TRUE, dim(e$migrm)[2])	
	for(sex in c('m', 'f')) {
		par <- paste0('migr',sex, what)
		for(age in dimnames(e$migrm)[[1]]) {
			i <- 1
			this.sumpop <- sumpop
			this.pop <- pop
			while(i < 100) {		
				dif <- sum(e[[par]][age,])
				dif.countries <- dif/this.sumpop * this.pop
				e[[par]][age,] <- e[[par]][age,] - dif.countries
				#if(what != '') break
				wneg <- which(e[[par]][age,] + pop.by.age[[sex]][age,] < -1e-4)
				this.mask <- mask[wneg]
				if(sum(this.mask)>0) {
					#stop('')
					e[[par]][age,mask] <- pmax(e[[par]][age,mask], -pop.by.age[[sex]][age,mask])
					this.pop[wneg[this.mask]] <- 0
					this.sumpop <- sum(this.pop)
					i <- i+1
					which.negative <- c(which.negative, wneg[this.mask])
					if(sum(this.pop[mask]) == 0) stop('')
				} else break
			}
			if(i>=100) stop('')
		}
	}
	e$negatives <- which.negative
}

rebalance.migration <- function(e, pop, what='', pop.by.age, check.negatives=FALSE) {
	sumpop <- sum(pop)
	which.negative <- c()	
	for(sex in c('m', 'f')) {
		par <- paste0('migr',sex, what)
		for(age in dimnames(e[[par]])[[1]]) {
			i <- 1
			this.sumpop <- sumpop
			this.pop <- pop
			while(i < 100) {		
				dif <- sum(e[[par]][age,])
				dif.countries <- dif/this.sumpop * this.pop
				e[[par]][age,] <- e[[par]][age,] - dif.countries
				if(!check.negatives) break
				wneg <- which(e[[par]][age,] + pop.by.age[[sex]][age,] < -1e-4)
				if(length(wneg) > 0) {
					#stop('')
					e[[par]][age,] <- pmax(e[[par]][age,], -pop.by.age[[sex]][age,])
					this.pop[wneg] <- 0
					this.sumpop <- sum(this.pop)
					i <- i+1
					which.negative <- c(which.negative, wneg)
					if(sum(this.pop) == 0) stop('')
				} else break
			}
			if(i>=100) stop('')
		}
	}
	e$negatives <- which.negative
}


rebalance.migration2groups <- function(e, pop, itraj) {
	# save these for debugging purposes
	e$migrm.before <- e$migrm
	e$migrf.before <- e$migrf
	e$migrm.labor.before <- e$migrm.labor
	e$migrf.labor.before <- e$migrf.labor
	before <- sum(e$migrm + e$migrf)
	before.labor <- sum(e$migrm.labor + e$migrf.labor)

	#mask <- rep(TRUE, dim(e$migrm)[2])
	slabor <- 0
	if(!is.null(e$migrm.labor))
		slabor <- sum(e$migrm.labor)
	#if(slabor != 0)
	#	mask[as.integer(names(pop)) %in% labor.countries()] <- FALSE 
	rebalance.migration(e, pop, pop.by.age=e$pop.by.age) # rest of the world
	if(slabor != 0) { # balancing for labor countries
		pop.labor <- pop * (as.integer(names(pop)) %in% labor.countries())
		#pop.by.age <- e$pop.by.age
		#pop.by.age$m[1:21,] <- pop.by.age$m[1:21,] + e$migrm
		#pop.by.age$f[1:21,] <- pop.by.age$f[1:21,] + e$migrf
		rebalance.migration(e, pop.labor, pop.by.age=e$pop.by.age, what='.labor')
	}
	e$migrm.after <- e$migrm
	e$migrf.after <- e$migrf
	e$migrm.labor.after <- e$migrm.labor
	e$migrf.labor.after <- e$migrf.labor
	after <- sum(e$migrm + e$migrf)
	after.labor <- sum(e$migrm.labor + e$migrf.labor)
	if(abs(after + after.labor - (before+before.labor)) > 1e6) 
				cat('\n', itraj, ' rebalancing migration: total=', after + after.labor - (before+before.labor),
							', labor=', after.labor-before.labor)
	# set migr to sum in order to deal with negatives
	e$migrm <- e$migrm + e$migrm.labor
	e$migrf <- e$migrf + e$migrf.labor
	rebalance.migration(e, pop, pop.by.age=e$pop.by.age, check.negatives=TRUE)
	dif <- abs(e$migrm - (e$migrm.after + e$migrm.labor.after))
	if(max(dif) > 1000) {
		wc <- ceiling(which.max(dif)/21)
		cat('\n', itraj, ' assuring positive pop: dif=', max(dif), ' country ', names(pop)[wc])
		#if(max(dif) > 100000) stop('')
	}
	e$migrm.labor[] <- 0
	e$migrf.labor[] <- 0
}

labor.countries <- function()
	return(prop.labor.migration()$country_code)
	
prop.labor.migration <- function()
	return(data.frame(country_code=c(
			 784, #United Arab Emirates
			 48,  # Bahrain
			 414, #Kuwait
			 512, # Oman
			 634, # Qatar
			 682, # Saudi Arabia
			 50,  # Bangladesh
			 818, # Egypt
			 360, # Indonesia
			 356, # India
			 586, # Pakistan
			 608 # Philippines
			),
			proportion=c(
				 0.933289201047455,   
				0.911919749810222,  
				1.04330264032715,      
				1.40418731176528,   
				0.898226910319697,  
				0.967835511646786,    
				0.744625843641641,    
				0.634207391400709,   
				0.337688415314329,    
				0.325522552720403,  
				0.372455661022053,   
				0.277753819108149
				))
		)
prop.labor.migration.for.country <- function(code){
	props <- prop.labor.migration()
	return(props[props$country_code==code,'proportion'])
}



restructure.pop.data.and.compute.quantiles.old <- function(source.dir, dest.dir, npred, inputs, observed, kannisto, 
									present.and.proj.years, keep.vital.events=FALSE, verbose=FALSE){
	envs <- list()
	for(time in 1:npred) {
		envs[[time]] <- new.env()
		load(file.path(source.dir, paste0('pop_time_', time, '.rda')), envir=envs[[time]])
		if(keep.vital.events) 
			load(file.path(source.dir, paste0('vital_events_time_', time, '.rda')), envir=envs[[time]])
	}
	ncountries <- nrow(envs[[1]]$totp)
	country.codes <- rownames(envs[[1]]$totp)
	nr.traj <- ncol(envs[[1]]$totp)
	nvariants <- ncol(envs[[1]]$totp.hch)
	ages <- dimnames(envs[[1]]$totpm)[[1]]
	nages <- length(ages)
	mx.ages <- c(0,1,ages[2:nages])
	npredplus1 <- npred + 1
	present.and.proj.years.pop <- present.and.proj.years + 2
	if(keep.vital.events) {
		migMntraj <- dim(envs[[1]]$migm)[length(dim(envs[[1]]$migm))]
		migFntraj <- dim(envs[[1]]$migf)[length(dim(envs[[1]]$migf))]
	}
	quantiles.to.keep <- get.quantiles.to.keep()
	nquant <- length(quantiles.to.keep)
	quant.env <- new.env()
	with(quant.env, {
		PIs_cqp <- quantM <- quantF <- array(NA, c(ncountries, nquant, npredplus1),
						dimnames=list(country.codes, quantiles.to.keep, present.and.proj.years.pop))
		quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(ncountries, nages, nquant, npredplus1),
						dimnames=list(country.codes, ages, quantiles.to.keep, present.and.proj.years.pop))
		mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(ncountries, 2, npredplus1), 
						dimnames=list(country.codes, c('mean', 'sd'), present.and.proj.years.pop))
	})
	for(cidx in 1:ncountries) {
		if(verbose) cat(cidx, ', ')
		inpc <- inputs[[country.codes[cidx]]]
		totp <- matrix(NA, nrow=npredplus1, ncol=nr.traj, dimnames=list(present.and.proj.years.pop, NULL))
		totpm <- totpf <- array(NA, dim=c(27, npredplus1, nr.traj), dimnames=list(ages, present.and.proj.years.pop, NULL))
		totp.hch <- matrix(NA, nrow=npredplus1, ncol=nvariants, dimnames=list(present.and.proj.years.pop, NULL))
		totpm.hch <- totpf.hch <- array(NA, dim=c(27, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years.pop, NULL))
		# values from current year
		totp[1,] <- sum(inpc$POPm0) + sum(inpc$POPf0)
		totpm[1:length(inpc$POPm0),1,] <- as.matrix(inpc$POPm0)[,rep(1,nr.traj)]
		totpf[1:length(inpc$POPf0),1,] <- as.matrix(inpc$POPf0)[,rep(1,nr.traj)]
		totp.hch[1,] <- totp[1,1]
		totpm.hch[,1,] <- totpm[,1,rep(1,nvariants)]
		totpf.hch[,1,] <- totpf[,1,rep(1,nvariants)]
		if(keep.vital.events) {
			btm <- btf <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm <- deathsf <- array(0, dim=c(27, npredplus1, nr.traj), dimnames=list(ages, present.and.proj.years, NULL))
			asfert <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			btm.hch <- btf.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years, NULL))
			asfert.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			mxm <- mxf <- array(0, dim=c(28, npredplus1, nr.traj), dimnames=list(mx.ages, present.and.proj.years, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(28, npredplus1, nvariants), dimnames=list(mx.ages, present.and.proj.years, NULL))
			migm <- array(0, dim=c(27, npredplus1, migMntraj), dimnames=list(ages, present.and.proj.years, NULL))
			migf <- array(0, dim=c(27, npredplus1, migFntraj), dimnames=list(ages, present.and.proj.years, NULL))
			# values from current year
			obs <- observed[[country.codes[cidx]]]
			MxKan <- kannisto[[country.codes[cidx]]]
			repi <- rep(1,nr.traj) # index for repeating columns
			btm[1:dim(obs$btm)[1],1,] <- obs$btm[,dim(obs$btm)[2],repi]
			btf[1:dim(obs$btf)[1],1,] <- obs$btf[,dim(obs$btf)[2],repi]
			deathsm[1:dim(obs$deathsm)[1],1,] <- obs$deathsm[,dim(obs$deathsm)[2],repi]
			deathsf[1:dim(obs$deathsf)[1],1,] <- obs$deathsf[,dim(obs$deathsf)[2],repi]
			asfert[1:dim(obs$asfert)[1],1,] <- obs$asfert[,dim(obs$asfert)[2],repi]
			mxm[1:dim(MxKan[[1]]$mx)[1],1,] <- as.matrix(MxKan[[1]]$mx[,dim(MxKan[[1]]$mx)[2]])[repi]
			mxf[1:dim(MxKan[[2]]$mx)[1],1,] <- as.matrix(MxKan[[2]]$mx[,dim(MxKan[[2]]$mx)[2]])[repi]
			migm[1:dim(inpc$observed$MIGm)[1],1,] <- as.matrix(inpc$observed$MIGm[,dim(inpc$observed$MIGm)[2]])[,rep(1,migMntraj)]
			migf[1:dim(inpc$observed$MIGf)[1],1,] <- as.matrix(inpc$observed$MIGf[,dim(inpc$observed$MIGf)[2]])[,rep(1,migFntraj)]
			btm.hch[,1,] <- btm[,1,rep(1,nvariants)]
			btf.hch[,1,] <- btf[,1,rep(1,nvariants)]
			deathsm.hch[,1,] <- deathsm[,1,rep(1,nvariants)]
			deathsf.hch[,1,] <- deathsf[,1,rep(1,nvariants)]
			asfert.hch[,1,] <- asfert[,1,rep(1,nvariants)]
			mxm.hch[,1,] <- mxm[,1,rep(1,nvariants)]
			mxf.hch[,1,] <- mxf[,1,rep(1,nvariants)]
		}
		for(time in 1:npred) {
			totp[time+1,] <- envs[[time]]$totp[cidx,]
			totpm[,time+1,] <- envs[[time]]$totpm[,cidx,]
			totpf[,time+1,] <- envs[[time]]$totpf[,cidx,]
			totp.hch[time+1,] <- envs[[time]]$totp.hch[cidx,]
			totpm.hch[,time+1,] <- envs[[time]]$totpm.hch[,cidx,]
			totpf.hch[,time+1,] <- envs[[time]]$totpf.hch[,cidx,]
			if(keep.vital.events) {
				btm[,time+1,] <- envs[[time]]$btm[,cidx,]
				btf[,time+1,] <- envs[[time]]$btf[,cidx,]
				deathsm[,time+1,] <- envs[[time]]$deathsm[,cidx,]
				deathsf[,time+1,] <- envs[[time]]$deathsf[,cidx,]
				asfert[,time+1,] <- envs[[time]]$asfert[,cidx,]
				mxm[,time+1,] <- envs[[time]]$mxm[,cidx,]
				mxf[,time+1,] <- envs[[time]]$mxf[,cidx,]
				migm[,time+1,] <- envs[[time]]$migm[,cidx,]
				migf[,time+1,] <- envs[[time]]$migf[,cidx,]
				btm.hch[,time+1,] <- envs[[time]]$btm.hch[,cidx,]
				btf.hch[,time+1,] <- envs[[time]]$btf.hch[,cidx,]
				deathsm.hch[,time+1,] <- envs[[time]]$deathsm.hch[,cidx,]
				deathsf.hch[,time+1,] <- envs[[time]]$deathsf.hch[,cidx,]
				asfert.hch[,time+1,] <- envs[[time]]$asfert.hch[,cidx,]
				mxm.hch[,time+1,] <- envs[[time]]$mxm.hch[,cidx,]
				mxf.hch[,time+1,] <- envs[[time]]$mxf.hch[,cidx,]
			}
		}
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch,
			 file = file.path(dest.dir, paste('totpop_country', country.codes[cidx], '.rda', sep='')))
		if(keep.vital.events) 
			save(btm, btf, deathsm, deathsf, asfert, mxm, mxf, migm, migf,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, 
				mxm.hch, mxf.hch, observed,
				file=file.path(dest.dir, paste('vital_events_country', country.codes[cidx], '.rda', sep='')))
		with(quant.env, {
			PIs_cqp[cidx,,] = apply(totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			mean_sd[cidx,1,] <- apply(totp, 1, mean, na.rm = TRUE)
			mean_sd[cidx,2,] = apply(totp, 1, sd, na.rm = TRUE)
			for (i in 1:nages) {
				if(nr.traj == 1) {
					quantMage[cidx,i,,] <- matrix(rep(totpm[i,,1],nquant) , nrow=nquant, byrow=TRUE)
					quantFage[cidx,i,,] <- matrix(rep(totpf[i,,1],nquant) , nrow=nquant, byrow=TRUE)
					quantPropMage[cidx,i,,] <- matrix(rep(totpm[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
					quantPropFage[cidx,i,,] <- matrix(rep(totpf[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
				} else {
					quantMage[cidx,i,,] <- apply(totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
					quantFage[cidx,i,,] <- apply(totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
					quantPropMage[cidx,i,,] <- apply(totpm[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
					quantPropFage[cidx,i,,] <- apply(totpf[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
				}
			}
			stotpm <- colSums(totpm, na.rm=TRUE)
			quantM[cidx,,] = apply(stotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			mean_sdM[cidx,1,] <- apply(stotpm, 1, mean, na.rm = TRUE)
			mean_sdM[cidx,2,] = apply(stotpm, 1, sd, na.rm = TRUE)
			stotpf <- colSums(totpf, na.rm=TRUE)
			quantF[cidx,,] = apply(stotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			mean_sdF[cidx,1,] <- apply(stotpf, 1, mean, na.rm = TRUE)
			mean_sdF[cidx,2,] = apply(stotpf, 1, sd, na.rm = TRUE)
		})
	}
	return(quant.env)
}

restructure.pop.data.and.compute.quantiles <- function(source.dir, dest.dir, npred, inputs, observed, kannisto, 
									present.and.proj.years, keep.vital.events=FALSE, parallel=FALSE, nr.nodes=NULL, verbose=FALSE, ...){
	envs <- list()
	for(time in 1:npred) {
		envs[[time]] <- new.env()
		load(file.path(source.dir, paste0('pop_time_', time, '.rda')), envir=envs[[time]])
		if(keep.vital.events) 
			load(file.path(source.dir, paste0('vital_events_time_', time, '.rda')), envir=envs[[time]])
	}
	ncountries <- nrow(envs[[1]]$totp)
	country.codes <- rownames(envs[[1]]$totp)
	nr.traj <- ncol(envs[[1]]$totp)
	nvariants <- ncol(envs[[1]]$totp.hch)
	ages <- dimnames(envs[[1]]$totpm)[[1]]
	nages <- length(ages)
	mx.ages <- c(0,1,ages[2:nages])
	npredplus1 <- npred + 1
	present.and.proj.years.pop <- present.and.proj.years + 2
	migMntraj <- migFntraj <- NULL
	if(keep.vital.events) {
		migMntraj <- dim(envs[[1]]$migm)[length(dim(envs[[1]]$migm))]
		migFntraj <- dim(envs[[1]]$migf)[length(dim(envs[[1]]$migf))]
	}
	quantiles.to.keep <- get.quantiles.to.keep()
	nquant <- length(quantiles.to.keep)
	quant.env <- new.env()
	with(quant.env, {
		PIs_cqp <- quantM <- quantF <- array(NA, c(ncountries, nquant, npredplus1),
						dimnames=list(country.codes, quantiles.to.keep, present.and.proj.years.pop))
		quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(ncountries, nages, nquant, npredplus1),
						dimnames=list(country.codes, ages, quantiles.to.keep, present.and.proj.years.pop))
		mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(ncountries, 2, npredplus1), 
						dimnames=list(country.codes, c('mean', 'sd'), present.and.proj.years.pop))
	})
	if(parallel) {
		if(verbose) cat('(in parallel on ', nr.nodes, ' nodes).')
		cl <- create.pop.cluster(nr.nodes, ...)
		clusterExport(cl, c("time", "nr.traj", "dest.dir", "country.codes", "inputs", "kannisto", 
								 "ages", "mx.ages", "npred", "nvariants", "observed", "present.and.proj.years",
								"migMntraj", "migFntraj", "keep.vital.events", "verbose"), envir=environment())
		res.list <- parLapplyLB(cl, 1:ncountries, 
							function(cidx) restructure.pop.data.and.compute.quantiles.one.country(cidx, country.codes[cidx], 
												envs, dest.dir, nr.traj, npred, nvariants, inputs[[country.codes[cidx]]], 
												observed[[country.codes[cidx]]], kannisto[[country.codes[cidx]]], 
												present.and.proj.years, ages, mx.ages, migMntraj, migFntraj,
												keep.vital.events=keep.vital.events, verbose=verbose))
		stopCluster(cl)
	} else { # process sequentially
		if(verbose) cat('(sequentially) ... ')
		res.list <- list()
		for(cidx in 1:ncountries) {
				if(verbose) cat(cidx, ', ')
				res.list[[cidx]] <- restructure.pop.data.and.compute.quantiles.one.country(cidx, country.codes[cidx], 
												envs, dest.dir, nr.traj, npred, nvariants, inputs[[country.codes[cidx]]], 
												observed[[country.codes[cidx]]], kannisto[[country.codes[cidx]]], 
												present.and.proj.years, ages, mx.ages, migMntraj, migFntraj,
												keep.vital.events=keep.vital.events, verbose=verbose)
		}
	}
	
	for(cidx in 1:ncountries) {
		for(par in c('PIs_cqp', 'mean_sd', 'quantM', 'quantF', 'mean_sdM', 'mean_sdF'))
			quant.env[[par]][cidx,,] <- res.list[[cidx]][[par]]
		for(par in c('quantMage', 'quantFage', 'quantPropMage', 'quantPropFage'))
			quant.env[[par]][cidx,,,] <- res.list[[cidx]][[par]]
	}
	return(quant.env)
}


restructure.pop.data.and.compute.quantiles.one.country <- function(cidx, country, envs, dest.dir, nr.traj, npred, nvariants, 
																	inpc, obs, MxKan, present.and.proj.years, 
																   ages, mx.ages, migMntraj, migFntraj,
																   keep.vital.events=FALSE, verbose=FALSE){
	npredplus1 <- npred + 1
	present.and.proj.years.pop <- present.and.proj.years + 2
	nages <- length(ages)
	quantiles.to.keep <- get.quantiles.to.keep()
	nquant <- length(quantiles.to.keep)
	quant.env <- new.env()
	with(quant.env, {
		PIs_cqp <- quantM <- quantF <- matrix(NA, nrow=nquant, ncol=npredplus1,
						dimnames=list(quantiles.to.keep, present.and.proj.years.pop))
		quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(nages, nquant, npredplus1),
						dimnames=list(ages, quantiles.to.keep, present.and.proj.years.pop))
		mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(2, npredplus1), 
						dimnames=list(c('mean', 'sd'), present.and.proj.years.pop))
	})
	repi <- rep(1,nr.traj) # index for repeating columns
	res.env <- new.env()
	with(res.env, {
		totp <- matrix(NA, nrow=npredplus1, ncol=nr.traj, dimnames=list(present.and.proj.years.pop, NULL))
		totpm <- totpf <- array(NA, dim=c(27, npredplus1, nr.traj), dimnames=list(ages, present.and.proj.years.pop, NULL))
		totp.hch <- matrix(NA, nrow=npredplus1, ncol=nvariants, dimnames=list(present.and.proj.years.pop, NULL))
		totpm.hch <- totpf.hch <- array(NA, dim=c(27, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years.pop, NULL))
		# values from current year
		totp[1,] <- sum(inpc$POPm0) + sum(inpc$POPf0)
		totpm[1:length(inpc$POPm0),1,] <- as.matrix(inpc$POPm0)[,rep(1,nr.traj)]
		totpf[1:length(inpc$POPf0),1,] <- as.matrix(inpc$POPf0)[,rep(1,nr.traj)]
		totp.hch[1,] <- totp[1,1]
		totpm.hch[,1,] <- totpm[,1,rep(1,nvariants)]
		totpf.hch[,1,] <- totpf[,1,rep(1,nvariants)]
		if(keep.vital.events) {
			btm <- btf <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm <- deathsf <- array(0, dim=c(27, npredplus1, nr.traj), dimnames=list(ages, present.and.proj.years, NULL))
			asfert <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			btm.hch <- btf.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years, NULL))
			asfert.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			mxm <- mxf <- array(0, dim=c(28, npredplus1, nr.traj), dimnames=list(mx.ages, present.and.proj.years, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(28, npredplus1, nvariants), dimnames=list(mx.ages, present.and.proj.years, NULL))
			migm <- array(0, dim=c(27, npredplus1, migMntraj), dimnames=list(ages, present.and.proj.years, NULL))
			migf <- array(0, dim=c(27, npredplus1, migFntraj), dimnames=list(ages, present.and.proj.years, NULL))
			# values from current year		
			btm[1:dim(obs$btm)[1],1,] <- obs$btm[,dim(obs$btm)[2],repi]
			btf[1:dim(obs$btf)[1],1,] <- obs$btf[,dim(obs$btf)[2],repi]
			deathsm[1:dim(obs$deathsm)[1],1,] <- obs$deathsm[,dim(obs$deathsm)[2],repi]
			deathsf[1:dim(obs$deathsf)[1],1,] <- obs$deathsf[,dim(obs$deathsf)[2],repi]
			asfert[1:dim(obs$asfert)[1],1,] <- obs$asfert[,dim(obs$asfert)[2],repi]
			mxm[1:dim(MxKan[[1]]$mx)[1],1,] <- as.matrix(MxKan[[1]]$mx[,dim(MxKan[[1]]$mx)[2]])[repi]
			mxf[1:dim(MxKan[[2]]$mx)[1],1,] <- as.matrix(MxKan[[2]]$mx[,dim(MxKan[[2]]$mx)[2]])[repi]
			migm[1:dim(inpc$observed$MIGm)[1],1,] <- as.matrix(inpc$observed$MIGm[,dim(inpc$observed$MIGm)[2]])[,rep(1,migMntraj)]
			migf[1:dim(inpc$observed$MIGf)[1],1,] <- as.matrix(inpc$observed$MIGf[,dim(inpc$observed$MIGf)[2]])[,rep(1,migFntraj)]
			btm.hch[,1,] <- btm[,1,rep(1,nvariants)]
			btf.hch[,1,] <- btf[,1,rep(1,nvariants)]
			deathsm.hch[,1,] <- deathsm[,1,rep(1,nvariants)]
			deathsf.hch[,1,] <- deathsf[,1,rep(1,nvariants)]
			asfert.hch[,1,] <- asfert[,1,rep(1,nvariants)]
			mxm.hch[,1,] <- mxm[,1,rep(1,nvariants)]
			mxf.hch[,1,] <- mxf[,1,rep(1,nvariants)]
		}
	})
	for(time in 1:npred) {
		for(par in c('totp', 'totp.hch'))
			res.env[[par]][time+1,] <- envs[[time]][[par]][cidx,]
		for(par in c('totpm', 'totpf', 'totpm.hch', 'totpf.hch'))
			res.env[[par]][,time+1,] <- envs[[time]][[par]][,cidx,]
		if(keep.vital.events) {
			for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'mxm', 'mxf', 'migm', 'migf',
							'btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'mxm.hch', 'mxf.hch'))
				res.env[[par]][,time+1,] <- envs[[time]][[par]][,cidx,]
		}
	}
	observed <- obs
	with(res.env, {
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch,
			 file = file.path(dest.dir, paste('totpop_country', country, '.rda', sep='')))
		if(keep.vital.events) 
			save(btm, btf, deathsm, deathsf, asfert, mxm, mxf, migm, migf,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, 
				mxm.hch, mxf.hch, observed,
				file=file.path(dest.dir, paste('vital_events_country', country, '.rda', sep='')))
	})
	stotpm <- colSums(res.env$totpm, na.rm=TRUE)
	stotpf <- colSums(res.env$totpf, na.rm=TRUE)
	with(quant.env, {
		PIs_cqp[,] <- apply(res.env$totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sd[1,] <- apply(res.env$totp, 1, mean, na.rm = TRUE)
		mean_sd[2,] <- apply(res.env$totp, 1, sd, na.rm = TRUE)
		for (i in 1:nages) {
			if(nr.traj == 1) {
				quantMage[i,,] <- matrix(rep(res.env$totpm[i,,1],nquant) , nrow=nquant, byrow=TRUE)
				quantFage[i,,] <- matrix(rep(res.env$totpf[i,,1],nquant) , nrow=nquant, byrow=TRUE)
				quantPropMage[i,,] <- matrix(rep(res.env$totpm[i,,1]/res.env$totp,nquant) , nrow=nquant, byrow=TRUE)
				quantPropFage[i,,] <- matrix(rep(res.env$totpf[i,,1]/res.env$totp,nquant) , nrow=nquant, byrow=TRUE)
			} else {
				quantMage[i,,] <- apply(res.env$totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantFage[i,,] <- apply(res.env$totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropMage[i,,] <- apply(res.env$totpm[i,,]/res.env$totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropFage[i,,] <- apply(res.env$totpf[i,,]/res.env$totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			}
		}
		quantM[,] = apply(stotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdM[1,] <- apply(stotpm, 1, mean, na.rm = TRUE)
		mean_sdM[2,] = apply(stotpm, 1, sd, na.rm = TRUE)
		quantF[,] = apply(stotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdF[1,] <- apply(stotpf, 1, mean, na.rm = TRUE)
		mean_sdF[2,] = apply(stotpf, 1, sd, na.rm = TRUE)
	})
	return(quant.env)
}

migration.age.schedule <- function(country, npred, inputs) {
	# original code by Jon Azose
	nAgeGroups <- 21
	#Initialize male and female matrices.
	maleArray <- matrix(0, nrow=nAgeGroups, ncol=npred)
	femaleArray <- matrix(0, nrow=nAgeGroups, ncol=npred)

	#Handle Croatia separately because of some bad data
	#Use 2010-2015 schedules, which aren't messed up.
	if(country == 191) {
		maleVec <- inputs$MIGm[inputs$MIGm$country_code==191,"2010-2015"]
		femaleVec <- inputs$MIGf[inputs$MIGf$country_code==191,"2010-2015"];
		tot <- sum(maleVec+femaleVec)
		#Set all future migration schedules for Croatia to match that one.
		croatiaM <- matrix(rep(maleVec/tot, npred), nrow=nAgeGroups)
		croatiaF <- matrix(rep(femaleVec/tot, npred), nrow=nAgeGroups)		
		return(list(M=croatiaM, F=croatiaF))
	}
	cidxM <- which(inputs$MIGm$country_code==country)
	cidxF <- which(inputs$MIGf$country_code==country)
	start.idx <- which(substr(colnames(inputs$MIGm), 1,4)==as.character(inputs$present.year))
	maleVec <- as.matrix(inputs$MIGm[cidxM,start.idx:ncol(inputs$MIGm)])
    femaleVec <- as.matrix(inputs$MIGf[cidxF,start.idx:ncol(inputs$MIGf)])
    colnames(maleVec) <- colnames(femaleVec) <- colnames(inputs$MIGm)[start.idx:ncol(inputs$MIGm)]
    tot <- colSums(maleVec+femaleVec)
    if(any(tot == 0)) {
		#Pull a model schedule to use in scenarios where the projection is 0
		#Use China's 2010-2015 data as the model
		modelmaleVec <- inputs$MIGm[inputs$MIGm$country_code==156, "2010-2015"]
		modelfemaleVec <- inputs$MIGf[inputs$MIGf$country_code==156, "2010-2015"]
		modeltot <- sum(modelmaleVec+modelfemaleVec)
		modelM <- modelmaleVec/modeltot
		modelF <- modelfemaleVec/modeltot
	}
	non.zero <- tot != 0
	if(any(non.zero)) {
		maleArray[,which(non.zero)] <- t(apply(maleVec[,which(non.zero), drop=FALSE], 1, '/', tot[which(non.zero)]))
    	femaleArray[,which(non.zero)] <- t(apply(femaleVec[,which(non.zero), drop=FALSE], 1, '/', tot[which(non.zero)]))
    }
    if(any(!non.zero)) {
    	maleArray[,which(!non.zero)] <- matrix(modelM, nrow=nAgeGroups, ncol=sum(!non.zero))
    	femaleArray[,which(!non.zero)] <- matrix(modelF, nrow=nAgeGroups, ncol=sum(!non.zero))
    }
    # special handling for negative rates
    negM <- negF <- NULL
    if(is.gcc(country)) { # set negatives to zero
    	negM <- maleArray
    	negM[negM<0] <- 0
    	negF <- femaleArray
    	negF[negF<0] <- 0
    	tot <- sum(negM+negF)
    	negM <- negM/tot
    	negF <- negF/tot
    }
    if(country == 818) { # Egypt (special handling when rate positive)
    	negM <- maleArray
    	maleArray[maleArray<0] <- 0
    	negF <- femaleArray
    	femaleArray[femaleArray<0] <- 0
    	tot <- sum(maleArray + femaleArray)
    	maleArray <- maleArray/tot
    	femaleArray <- femaleArray/tot
    }
    if(country %in% c(528, 756)) { # Netherlands and Switzerland (get Czech schedule)
    	maleVec <- inputs$MIGm[inputs$MIGm$country_code==203, "2010-2015"]
    	femaleVec <- inputs$MIGf[inputs$MIGf$country_code==203, "2010-2015"]
    	tot <- sum(maleVec+femaleVec)
    	negM <- matrix(maleVec/tot, nrow=nAgeGroups, ncol=ncol(maleArray))
    	negF <- matrix(femaleVec/tot, nrow=nAgeGroups, ncol=ncol(femaleArray))
    }
	return(list(M=maleArray, F=femaleArray, Mnegative=negM, Fnegative=negF))
}

PopProjNoMigr <- function(npred, pop0, LT, asfr, srb, country.name=NULL, keep.vital.events=FALSE) {
	popm <- popf <- matrix(0, nrow=27, ncol=npred+1)
	popm[,1] <- c(pop0$M, rep(0, 27-length(pop0$M)))
	popf[,1] <- c(pop0$F, rep(0, 27-length(pop0$F)))
	totp <- c(sum(popm[,1]+popf[,1]), rep(0, npred))
	btageM <- btageF <- matrix(0, nrow=7, ncol=npred) # births by age of mother and sex of child
	deathsM <- deathsF <- matrix(0, nrow=27, ncol=npred)
	nproj <- npred
	returnIfNegative <- 0
	isNegative <- 0
	#while(TRUE) {
		res <- .C("PopProjNoMigration", as.integer(nproj), 
			srm=LT$sr[[1]], srf=LT$sr[[2]], asfr=as.numeric(as.matrix(asfr)), 
			srb=as.numeric(as.matrix(srb)), returnNothingIfNegative=as.integer(returnIfNegative), 
			popm=popm, popf=popf, totp=totp,
			Lm=as.numeric(LT$LLm[[1]]), Lf=as.numeric(LT$LLm[[2]]), 
			lxm=as.numeric(LT$lx[[1]]), lxf=as.numeric(LT$lx[[2]]), 
			btagem=as.numeric(btageM), btagef=as.numeric(btageF), 
			deathsm=as.numeric(deathsM), deathsf=as.numeric(deathsF),
			isNegative=as.integer(isNegative)
			)
		#if(returnIfNegative==1 && res$isNegative < 0) {
		#	warning('Negative population for ', country.name, '. Counts adjusted.', immediate.=TRUE)
		#	returnIfNegative <- 0
		#} else break
	#}
	vital.events <- list()
	if(keep.vital.events) {
		vital.events$mbt <- res$btagem
		vital.events$fbt <- res$btagef
		vital.events$mdeaths <- res$deathsm
		vital.events$fdeaths <- res$deathsf
	}
	return(c(list(totpop=res$totp, mpop=res$popm, fpop=res$popf), vital.events))
}


migration.age.schedules.not.used <- function(countries, npred, inputs) {
	# original code by Jon Azose
	nAgeGroups <- 21
	nCountries <- length(countries)
	countries.idx <- 1:nCountries
	#Initialize male and female matrices.
	maleArray <- array(0,dim=c(nAgeGroups,nCountries,npred))
	femaleArray <- array(0,dim=c(nAgeGroups,nCountries,npred))

	#First pull a model schedule to use in scenarios where the projection is 0
	#Use China's 2010-2015 data as the model
	country <- if(!(156 %in% countries)) countries[1] else 156 # if China not available, take the first country
	maleVec <- inputs$MIGm[inputs$MIGm$country_code==country, "2010-2015"]
	femaleVec <- inputs$MIGf[inputs$MIGf$country_code==country, "2010-2015"]
	tot <- sum(maleVec+femaleVec)
	modelM <- maleVec/tot
	modelF <- femaleVec/tot

	#Handle Croatia separately because of some bad data
	#Use 2010-2015 schedules, which aren't messed up.
	if(191 %in% countries) {
		maleVec <- inputs$MIGm[inputs$MIGm$country_code==191,"2010-2015"]
		femaleVec <- inputs$MIGf[inputs$MIGf$country_code==191,"2010-2015"];
		tot <- sum(maleVec+femaleVec)
		croatiaM <- maleVec/tot
		croatiaF <- femaleVec/tot
		#Now set all future migration schedules for Croatia to match that one.
		idx <- which(countries == 191)
  		maleArray[,idx,] <- croatiaM
  		femaleArray[,idx,] <- croatiaF
  		countries.idx <- countries.idx[-idx]
	}
	#Now process the rest of the countries (not Croatia)
	start.idx <- which(as.integer(substr(colnames(inputs$MIGm), 1,4))==inputs$present.year)
	for(i in countries.idx){
		for(j in 1:npred){
    		maleVec <- inputs$MIGm[inputs$MIGm$country_code==countries[i],start.idx+j]
    		femaleVec <- inputs$MIGf[inputs$MIGf$country_code==countries[i],start.idx+j]
		    tot <- sum(maleVec+femaleVec)
		    #If the projection is for non-zero migration:
    		if(tot != 0){ #Normalize to sum to 1			
      			maleArray[,i,j] <- maleVec/tot
      			femaleArray[,i,j] <- femaleVec/tot
		    } else { #Otherwise, just use model schedules
      			maleArray[,i,j] <- modelM
      			femaleArray[,i,j] <- modelF
    		}
  		}
	}
	return(list(male=maleArray, female=femaleArray))
}



do.pop.predict.balance.old <- function(inp, outdir, nr.traj, ages, pred=NULL, keep.vital.events=FALSE, function.inputs=NULL, start.time.index=1, 
									verbose=FALSE, parallel=FALSE, nr.nodes=NULL, ...) {
	countries.idx <- which(UNlocations$location_type==4)
	#countries.idx <- which(UNlocations$country_code %in% c(716, 250))
	country.codes <- UNlocations$country_code[countries.idx]
	ncountries <- length(country.codes)
	nr_project <- length(inp$proj.years)
	nages <- length(ages)
	mx.ages <- c(0,1,ages[2:nages])
	ages21 <- ages[ages<=100]
	nest <- length(inp$estim.years)
	if(!file.exists(outdir)) 
		dir.create(outdir, recursive=TRUE)
	present.and.proj.years <- c(inp$estim.years[nest], inp$proj.years)
	
	status.for.gui <- paste('out of', nr_project, 'time periods.')
	gui.options <- list()
	inp.to.save <- list()
	# remove big or redundant items from inputs to be saved
	for(item in ls(inp)[!grepl('^migMpred$|^migFpred$|^TFRpred$|^e0Fpred$|^e0Mpred$|^estim.years$|^proj.years$|^wpp.years$', ls(inp))]) 
		inp.to.save[[item]] <- get(item, inp)
	npred <- nr_project
	
	outdir.tmp <- file.path(outdir, '_tmp_')
	if(file.exists(outdir.tmp) && start.time.index==1) unlink(outdir.tmp, recursive=TRUE)
	if(start.time.index==1) dir.create(outdir.tmp)
	
	if(start.time.index > 1) { # reload initial population
		env.tmp <- new.env()
		load(file.path(outdir.tmp, paste0('pop_time_', start.time.index-1, '.rda')), envir=env.tmp)
		popM.prev <- env.tmp$totpm
		popF.prev <- env.tmp$totpf
		popM.hch.prev <- env.tmp$totpm.hch
		popF.hch.prev <- env.tmp$totpf.hch
	}
	if(parallel && is.null(nr.nodes)) nr.nodes <- getOption("cl.cores", detectCores())
	countries.input <- new.env()
	# Extract the country-specific stuff from the inputs
	if(verbose) cat('\nLoading inputs for ', ncountries, ' countries ')
	if(parallel) {
		if(verbose) cat('(in parallel on ', nr.nodes, ' nodes).')
		cl <- create.pop.cluster(nr.nodes, ...)
		clusterExport(cl, c("inp", "nr.traj", "country.codes", "countries.idx", "UNlocations"), envir=environment())
		inpc.list <- parLapplyLB(cl, 1:ncountries, 
						function(cidx) get.country.inputs(country.codes[cidx], inp, nr.traj, UNlocations[countries.idx[cidx],'name']))
		stopCluster(cl)
		for(cidx in 1:ncountries) {	
			if(is.null(inpc.list[[cidx]]) || length(inpc.list[[cidx]]$POPm0)==0) next
			countries.input[[as.character(country.codes[cidx])]] <- inpc.list[[cidx]]
		}
	} else { # load inputs sequentially
		if(verbose) cat('(sequentially).')
		for(cidx in 1:ncountries) {
			inpc <- get.country.inputs(country.codes[cidx], inp, nr.traj, UNlocations[countries.idx[cidx],'name'])
			if(is.null(inpc) || length(inpc$POPm0)==0) next
			countries.input[[as.character(country.codes[cidx])]] <- inpc
		}
	}
	nr.traj <- min(nr.traj, sapply(ls(countries.input), function(ccode) return(ncol(countries.input[[ccode]]$TFRpred))))
	npred <- min(npred, sapply(ls(countries.input), function(ccode) return(nrow(countries.input[[ccode]]$TFRpred))))
	npredplus1 <- npred+1
	npasfr <- nrow(countries.input[[as.character(country.codes[1])]]$PASFR)
	nvariants <- nrow(countries.input[[as.character(country.codes[1])]]$TFRhalfchild)
	if(verbose) cat('\nProcessing ', nr.traj, ' trajectories for each country.')
	if(length(countries.input) < ncountries) {
		ncountries <- length(countries.input)
		country.codes <- as.integer(ls(countries.input))
		countries.idx <- sapply(country.codes, function(x) which(UNlocations[,'country_code'] == x))
	}
	kannisto <- list()
	observed <- new.env()
	for(country in ls(countries.input)) {
		kannisto[[country]] <- runKannisto(nest, countries.input[[country]])
		if(keep.vital.events) 
			observed[[country]] <- compute.observedVE(countries.input[[country]], inp$pop.matrix, 
										countries.input[[country]]$MIGtype, kannisto[[country]], 
										as.integer(country), inp$estim.years)
	}
	debug <- FALSE
	for(time in start.time.index:npred) {	
		unblock.gtk.if.needed(paste('finished', time, status.for.gui), gui.options)
		if(verbose)
			cat('\nProcessing time period ', time)
		totp <- matrix(NA, nrow=ncountries, ncol=nr.traj, dimnames=list(country.codes, NULL))
		totpm <- totpf <- array(NA, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
		migrationm <- migrationf <- array(NA, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
		totp.hch <- matrix(NA, nrow=ncountries, ncol=nvariants, dimnames=list(country.codes, NULL))
		totpm.hch <- totpf.hch <- array(NA, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
		migrationm.hch <- migrationf.hch <- array(NA, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
		if(keep.vital.events) {
			btm <- btf <- asfert <- array(0, dim=c(7, ncountries, nr.traj), dimnames=list(NULL, country.codes, NULL))
			deathsm <- deathsf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
			btm.hch <- btf.hch <- asfert.hch <- array(0, dim=c(7, ncountries, nvariants), dimnames=list(NULL, country.codes, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
			mxm <- mxf <- array(0, dim=c(28, ncountries, nr.traj), dimnames=list(mx.ages, country.codes, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(28, ncountries, nvariants), dimnames=list(mx.ages, country.codes, NULL))
			migMntraj <- if(is.null(countries.input[[as.character(country.codes[1])]][['migMpred']])) 1 
						else dim(countries.input[[as.character(country.codes[1])]][['migMpred']])[2]
			migFntraj <- if(is.null(countries.input[[as.character(country.codes[1])]][['migFpred']])) 1 
						else dim(countries.input[[as.character(country.codes[1])]][['migFpred']])[2]  
			migm <- array(0, dim=c(27, ncountries, migMntraj), dimnames=list(ages, country.codes, NULL))
			migf <- array(0, dim=c(27, ncountries, migFntraj), dimnames=list(ages, country.codes, NULL))
		}
		for(cidx in 1:ncountries) {
			country <- country.codes[cidx]
			country.idx <- countries.idx[cidx]
			ccountry <- as.character(country)
			inpc <- countries.input[[ccountry]]
			pop.ini <- list(M=inpc$POPm0, F=inpc$POPf0)
			for(itraj in 1:nr.traj) {
				asfr <- inpc$PASFR[,time,drop=FALSE]/100.
				for(i in 1:npasfr) asfr[i,] <- inpc$TFRpred[time,itraj] * asfr[i,]
				#if(country==716 && itraj==61) debug <- TRUE
				#cat('\n\n\nTrajectory: ', itraj, '\n========\n\n')
				LTres <- modifiedLC(1, kannisto[[ccountry]], inpc$e0Mpred[time,itraj], 
									inpc$e0Fpred[time,itraj], verbose=verbose, debug=debug)
				migpred <- .get.migration.one.trajectory(inpc, itraj, time)			
				if(time > 1) { # reset initial population to the one at the previous time step
					pop.ini$M <- popM.prev[,cidx, itraj]
					pop.ini$F <- popF.prev[,cidx, itraj]
				}
				popres <- StoPopProj(1, pop.ini, LTres, asfr, inpc$SRB[time], migpred, inpc$MIGtype, country.name=UNlocations[country.idx,'name'],
									keep.vital.events=keep.vital.events)
				totp[cidx,itraj] <- popres$totpop[2]
				totpm[,cidx,itraj] <- popres$mpop[,2]
				totpf[,cidx,itraj] <- popres$fpop[,2]
				migrationm[,cidx,itraj] <- popres$mmigr # migpred$M
				migrationf[,cidx,itraj] <- popres$fmigr # migpred$F
				if(keep.vital.events) {
					btm[,cidx,itraj] <- popres$mbt
					btf[,cidx,itraj] <- popres$fbt
					deathsm[,cidx,itraj] <- popres$mdeaths
					deathsf[,cidx,itraj] <- popres$fdeaths
					asfert[,cidx,itraj] <- asfr
					mxm[,cidx,itraj] <- LTres$mx[[1]]
					mxf[,cidx,itraj] <- LTres$mx[[2]]
					migtraj <- min(itraj, migMntraj)
					migm[,cidx,migtraj] <- popres$mmigr # migpred[['M']]
					migtraj <- min(itraj, migFntraj)
					migf[,cidx,migtraj] <- popres$fmigr # migpred[['F']]
				}
			}
			pop.ini <- list(M=inpc$POPm0, F=inpc$POPf0)
			for (variant in 1:nvariants) { # compute the two half child variants
				asfr <- inpc$PASFR[,time,drop=FALSE]/100.
				for(i in 1:npasfr) asfr[i,] <- inpc$TFRhalfchild[variant,time] * asfr[i,]
				LTres <- modifiedLC(1, kannisto[[ccountry]], inpc$e0Mmedian[time], 
									inpc$e0Fmedian[time], verbose=verbose, debug=debug)
				migpred.hch <- .get.migration.one.trajectory(inpc, itraj=NULL, time=time)
				if(time > 1) {
					pop.ini$M <- popM.hch.prev[,cidx, variant]
					pop.ini$F <- popF.hch.prev[,cidx, variant]
				}
				popres <- StoPopProj(1, pop.ini, LTres, asfr, inpc$SRB[time], migpred.hch, inpc$MIGtype, 
							country.name=UNlocations[country.idx,'name'], 
							keep.vital.events=keep.vital.events)
				totp.hch[cidx,variant] <- popres$totpop[2]
				totpm.hch[,cidx,variant] <- popres$mpop[,2]
				totpf.hch[,cidx,variant] <- popres$fpop[,2]
				migrationm.hch[,cidx,variant] <- popres$mmigr # migpred.hch$M
				migrationf.hch[,cidx,variant] <- popres$fmigr # migpred.hch$F
				if(keep.vital.events) {
					btm.hch[,cidx,variant] <- popres$mbt
					btf.hch[,cidx,variant] <- popres$fbt
					deathsm.hch[,cidx,variant] <- popres$mdeaths
					deathsf.hch[,cidx,variant] <- popres$fdeaths
					asfert.hch[,cidx,variant] <- asfr
					mxm.hch[,cidx,variant] <- LTres$mx[[1]]
					mxf.hch[,cidx,variant] <- LTres$mx[[2]]
				}
			}
		} # end countries
		# Rebalancing (by trajectories)
		if(verbose) cat('\nRe-balance migration.')
		balance.env <- new.env()
		vars <- c('totp', 'totpm', 'totpf', 'migrationm', 'migrationf')
		for(par in vars) balance.env[[par]] <- get(par)
		rebalance.population.by.migration(balance.env)
		for(par in vars) # overwrite the objects in this environment
			assign(par, balance.env[[par]])
		for(par in vars) balance.env[[par]] <- get(paste0(par, '.hch'))
		rebalance.population.by.migration(balance.env)	
		for(par in vars) assign(paste0(par, '.hch'), balance.env[[par]])
		#stop('')
		popM.prev <- totpm
		popF.prev <- totpf
		popM.hch.prev <- totpm.hch
		popF.hch.prev <- totpf.hch
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch,
				 file = file.path(outdir.tmp, paste0('pop_time_', time, '.rda')))
		if(keep.vital.events) 
			save(btm, btf, deathsm, deathsf, asfert, mxm, mxf, migm, migf,
					btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, 
					mxm.hch, mxf.hch, 
					file=file.path(outdir.tmp, paste0('vital_events_time_', time, '.rda')))
	} # end time
	if(verbose) cat('\nRe-formatting data ...')
	quant.env <- restructure.pop.data.and.compute.quantiles(outdir.tmp, outdir, npred, countries.input, observed, kannisto, 
					present.and.proj.years, keep.vital.events, verbose=verbose)
	if(verbose) cat(' done.\n')
	unlink(outdir.tmp, recursive=TRUE)
	#save meta file
	country.rows <- UNlocations[countries.idx,c('country_code', 'name')]
	colnames(country.rows) <- c('code', 'name')
	bayesPop.prediction <- structure(list(
							nr.traj = nr.traj,	
							quantiles = quant.env$PIs_cqp,
               				traj.mean.sd = quant.env$mean_sd,
               				quantilesM = quant.env$quantM, 
               				traj.mean.sdM = quant.env$mean_sdM,
               				quantilesF = quant.env$quantF, 
               				traj.mean.sdF = quant.env$mean_sdF,
               				quantilesMage = quant.env$quantMage, 
               				quantilesFage = quant.env$quantFage, 
               				quantilesPropMage = quant.env$quantPropMage, 
               				quantilesPropFage = quant.env$quantPropFage,
               				estim.years=inp$estim.years, 
               				proj.years=present.and.proj.years, # includes present period (middle of periods)
               				proj.years.pop=present.and.proj.years+2, # end of periods
               				wpp.year = inp$wpp.year,
			   				inputs = inp.to.save, # save as list because environment takes much more space
			   				function.inputs=function.inputs,
			   				countries=country.rows,
			   				ages=ages), class='bayesPop.prediction')

	prediction.file <- file.path(outdir, 'prediction.rda')
	save(bayesPop.prediction, file=prediction.file)
	cat('\nPrediction stored into', outdir, '\n')
	return(bayesPop.prediction)
}
