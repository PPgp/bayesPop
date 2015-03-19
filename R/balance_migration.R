if(getRversion() >= "2.15.1") utils::globalVariables(c("land_area_wpp2012"))

do.pop.predict.balance <- function(inp, outdir, nr.traj, ages, pred=NULL, keep.vital.events=FALSE, function.inputs=NULL, 
									rebalance=TRUE, use.migration.model=TRUE, start.time.index=1, 
									verbose=FALSE, parallel=FALSE, nr.nodes=NULL, .countries=NULL, .migration.pars=NULL, ...) {
	# if .migration.pars is NULL, don't use the migration model
	if (use.migration.model && is.null(.migration.pars))
		stop('Argument .migration.pars has to be given if use.migration.model is TRUE.')
	countries.idx <- if(is.null(.countries)) which(UNlocations$location_type==4) else which(UNlocations$country_code %in% .countries)
	country.codes <- UNlocations$country_code[countries.idx]
	country.codes.char <- as.character(country.codes)
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
		
	mig.rate.prev <- mig.rate <- NULL
	if(use.migration.model) {
		if(is.null(.migration.pars$posterior))
			stop('Argument .migration.pars must have an element "posterior".')
		inp[['migration.parameters']] <- .migration.pars$posterior
		if(!is.null(.migration.pars$ini.rates)) {
			mig.rate.prev <- .migration.pars$ini.rates#/5.
		} else {
			mig.rate.prev <- inp$migration.rates[,ncol(inp$migration.rates)]
			names(mig.rate.prev) <- inp$migration.rates$country_code
		}
		inp$migration.rates <- inp$migration.rates[,-ncol(inp$migration.rates)] # remove the present year column as it is in mig.rate.prev
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
		if(!is.null(env.tmp$mig.rate)) {
			mig.rate <- env.tmp$mig.rate
			mig.rate.prev <- mig.rate[start.time.index-1,,]
		}
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
			countries.input[[country.codes.char[cidx]]] <- as.environment(inpc.list[[cidx]])
		}
	} else { # load inputs sequentially
		if(verbose) cat('(sequentially).')
		for(cidx in 1:ncountries) {
			inpc <- get.country.inputs(country.codes[cidx], inp, nr.traj, UNnames[cidx])
			if(is.null(inpc) || length(inpc$POPm0)==0) next
			countries.input[[country.codes.char[cidx]]] <- as.environment(inpc)
		}
	}
	nr.traj <- min(nr.traj, sapply(ls(countries.input), function(ccode) return(ncol(countries.input[[ccode]]$TFRpred))))
	npred <- min(npred, sapply(ls(countries.input), function(ccode) return(nrow(countries.input[[ccode]]$TFRpred))))
	npredplus1 <- npred+1
	npasfr <- nrow(countries.input[[country.codes.char[1]]]$PASFR)
	nvariants <- nrow(countries.input[[country.codes.char[1]]]$TFRhalfchild)
	if(verbose) cat('\nProcessing ', nr.traj, ' trajectories for each country ')
	if(length(countries.input) < ncountries) {
		ncountries <- length(countries.input)
		country.codes.char <- ls(countries.input)
		country.codes <- as.integer(country.codes.char)
		countries.idx <- sapply(country.codes, function(x) which(UNlocations[,'country_code'] == x))
		UNnames <- UNlocations[countries.idx,'name']
		if(parallel) nr.nodes <- min(nr.nodes, ncountries)
	}
	if(!is.null(mig.rate.prev)) {
		mig.rate.prev <- mig.rate.prev[country.codes.char]
		mig.rate.prev <- matrix(mig.rate.prev, nrow=length(mig.rate.prev), ncol=nr.traj)
		if(is.null(mig.rate)) 
			mig.rate <- array(NA, c(npred+1, ncountries, nr.traj), dimnames=list(NULL, country.codes, NULL))
		mig.rate[start.time.index,,] <- mig.rate.prev
	}
	kannisto <- list()
	observed <- new.env()
	kantor.pasfr <- new.env()
	
	for(country in country.codes.char) {
		kannisto[[country]] <- runKannisto(countries.input[[country]], inp$start.year, npred=npred)
		if(keep.vital.events) 
			observed[[country]] <- compute.observedVE(countries.input[[country]], inp$pop.matrix, 
										countries.input[[country]]$MIGtype, kannisto[[country]], 
										as.integer(country), inp$estim.years)
		tfr.med <- apply(countries.input[[country]]$TFRpred, 1, median)[nrow(countries.input[[country]]$TFRpred)]
		kantor.pasfr[[country]] <- list()
		for(itraj in 1:nr.traj)		
			kantor.pasfr[[country]][[itraj]] <- kantorova.pasfr(c(countries.input[[country]]$observed$TFRpred, countries.input[[country]]$TFRpred[,itraj]), 
										countries.input[[country]], norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)	
		for (variant in 1:nvariants) 
			kantor.pasfr[[country]][[variant + nr.traj]] <- kantorova.pasfr(
						c(countries.input[[country]]$observed$TFRpred, countries.input[[country]]$TFRhalfchild[variant,]), 
										countries.input[[country]], norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)					
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
			btm <- btf <- asfert <- pasfert <- array(0, dim=c(7, ncountries, nr.traj), dimnames=list(NULL, country.codes, NULL))
			deathsm <- deathsf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
			btm.hch <- btf.hch <- asfert.hch <- pasfert.hch <- array(0, dim=c(7, ncountries, nvariants), dimnames=list(NULL, country.codes, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
			mxm <- mxf <- array(0, dim=c(28, ncountries, nr.traj), dimnames=list(mx.ages, country.codes, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(28, ncountries, nvariants), dimnames=list(mx.ages, country.codes, NULL))
			migntraj <- countries.input[[as.character(country.codes[1])]][['mig.nr.traj']]
			migm <- array(0, dim=c(27, ncountries, migntraj), dimnames=list(ages, country.codes, NULL))
			migf <- array(0, dim=c(27, ncountries, migntraj), dimnames=list(ages, country.codes, NULL))
		}
	})
	debug <- FALSE
	res.env$mig.rate <- mig.rate
	if(parallel) {
		if(verbose) cat(' (in parallel on ', nr.nodes, ' nodes).')
		cl <- create.pop.cluster(nr.nodes, ...)
		clusterExport(cl, c("nr.traj", "country.codes",  "UNnames", "countries.input", "kannisto", "npasfr", 
								"ages", "nvariants", 
								"keep.vital.events", "verbose"), envir=environment())
	} else if(verbose) cat(' (sequentially).')
	
	wrapper.pop.predict.one.country.no.migration <- function(cidx) {
			if(time > 1) 
				do.pop.predict.one.country.no.migration(time, 
												UNnames[cidx], countries.input[[country.codes.char[cidx]]], 
												kannisto[[country.codes.char[cidx]]], kantor.pasfr[[country.codes.char[cidx]]], nr.traj, 
												popM.prev[,cidx, ,drop=FALSE], popF.prev[,cidx, ,drop=FALSE],
												popM.hch.prev[,cidx, ,drop=FALSE], popF.hch.prev[,cidx, ,drop=FALSE],
												ages, nvariants, 
												keep.vital.events=keep.vital.events, verbose=verbose)
			else do.pop.predict.one.country.no.migration(time, 
												UNnames[cidx], countries.input[[country.codes.char[cidx]]], 
												kannisto[[country.codes.char[cidx]]], kantor.pasfr[[country.codes.char[cidx]]], nr.traj,  
												ages=ages, nvariants=nvariants, keep.vital.events=keep.vital.events, verbose=verbose)
	} 
		
	for(time in start.time.index:npred) {
		unblock.gtk.if.needed(paste('finished', time, status.for.gui), gui.options)
		if(verbose)
			cat('\nProcessing time period ', time)
		.ini.pop.res.env(res.env, keep.vital.events)
		if(parallel) {
			clusterExport(cl, c("time", "mig.rate.prev"), envir=environment())
			if(time > 1) clusterExport(cl, c("popM.prev", "popF.prev", "popM.hch.prev", "popF.hch.prev"), envir=environment())
			res.list <- parLapplyLB(cl, 1:ncountries, wrapper.pop.predict.one.country.no.migration)
			
		} else { # process sequentially
			res.list <- list()
			for(cidx in 1:ncountries) res.list[[cidx]] <- wrapper.pop.predict.one.country.no.migration(cidx)
		}
		for(cidx in 1:ncountries) { # collect results
			for(par in c('totp', 'totp.hch'))
				res.env[[par]][cidx,] <- res.list[[cidx]][[par]]
			for(par in c('totpm', 'totpf', 'totpm.hch', 'totpf.hch'))#'migrationm', 'migrationf', 'migrationm.hch', 'migrationf.hch'
				res.env[[par]][,cidx,] <- res.list[[cidx]][[par]]
			if(keep.vital.events) {
				for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'pasfert', 'mxm', 'mxf', 
							'btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'pasfert.hch', 'mxm.hch','mxf.hch')) # 'migm', 'migf',
					res.env[[par]][,cidx,] <- res.list[[cidx]][[par]]
			}
		}
		migpred <- get.balanced.migration(time, country.codes, countries.input, nr.traj, rebalance, use.migration.model,
							if(time>1)popM.prev[,, ,drop=FALSE] else NULL, 
							if(time>1)popF.prev[,, ,drop=FALSE] else NULL, ages, res.env, verbose=verbose)
		lage <- dim(migpred$M)[1]
		res.env$migrationm[1:lage,,] <- migpred$M
		res.env$migrationf[1:lage,,] <- migpred$F
		res.env$migrationm.hch[1:lage,,] <- migpred$Mhch
		res.env$migrationf.hch[1:lage,,] <- migpred$Fhch
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
		
		with(res.env, {
			save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch, mig.rate,
				 file = file.path(outdir.tmp, paste0('pop_time_', time, '.rda')))
			if(keep.vital.events) 
				save(btm, btf, deathsm, deathsf, asfert, pasfert, mxm, mxf, migm, migf,
						btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, pasfert.hch, 
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

get.balanced.migration <- function(time, country.codes, inputs, nr.traj, rebalance, use.migration.model,
												popM.prev, popF.prev, ages, env, verbose=FALSE) {
	nr.countries <- length(country.codes)
	e <- new.env()
	labor.codes <- labor.countries()
	ages <- ages[1:21]
	lages <- length(ages)
	e$migrm <- e$migrf <- matrix(NA, ncol=nr.countries, nrow=lages, dimnames=list(ages, country.codes))
	e$migrm.labor <- e$migrf.labor <- matrix(0, ncol=nr.countries, nrow=lages, dimnames=list(ages, country.codes))
	migrationm <- migrationf <- array(NA, c(lages, nr.countries, nr.traj), dimnames=list(ages, country.codes, NULL))
	migrationm.hch <- migrationf.hch <- array(NA, c(lages, nr.countries, 2), dimnames=list(ages, country.codes, NULL))
	
	pop <- rep(NA, nr.countries)
	names(pop) <- country.codes
	country.codes.char <- as.character(country.codes)
	data(land_area_wpp2012)
	warns <- list()	
	for(itraj in 1:nr.traj){
		pop <- colSums(env$totpm[,,itraj] + env$totpf[,,itraj])
		for(cidx in 1:nr.countries) {
			inpc <- inputs[[country.codes.char[cidx]]]
			migpred <- .get.migration.one.trajectory(use.migration.model, inpc, itraj, time, pop[cidx], 
								popM=env$totpm[,cidx,itraj], popF=env$totpf[,cidx,itraj], country.code=country.codes[cidx], 
								mig.rates=if(!is.null(env$mig.rate)) env$mig.rate[,cidx,itraj] else NULL)
			e$migrm[,cidx] <- migpred$M
			e$migrf[,cidx] <- migpred$F
			if(!is.null(migpred$laborM)) {
				e$migrm.labor[,cidx] <- migpred$laborM
				e$migrf.labor[,cidx] <- migpred$laborF
			}
			env$mig.rate[time+1, cidx, itraj] <- migpred$rate
			warns[[country.codes.char[cidx]]] <- c(warns[[country.codes.char[cidx]]], migpred$warns)
		}
		e$popm <- env$totpm[,,itraj]
		e$popf <- env$totpf[,,itraj]
		if(rebalance) {
			rebalance.migration2groups(e, pop, itraj)			
			negatives <- as.character(country.codes[unique(e$negatives)])
			for(country in negatives)
				warns[[country]] <- c(warns[[country]], 'population negative while balancing')
		}
		migrationm[,,itraj] <- e$migrm + e$migrm.labor
		migrationf[,,itraj] <- e$migrf + e$migrf.labor
	}
	for(variant in 1:2){
		pop <- colSums(env$totpm[,,variant] + env$totpf[,,variant])
		for(cidx in 1:nr.countries) {
			inpc <- inputs[[as.character(country.codes[cidx])]]
			fixed.rate <- if(!is.null(env$mig.rate)) median(env$mig.rate[time,cidx,]) else NULL
			migpred <- .get.migration.one.trajectory(use.migration.model, inpc, variant, time, pop[cidx], 
								popM=env$totpm.hch[,cidx,variant], popF=env$totpf.hch[,cidx,variant],
								country.code=country.codes[cidx], fixed.rate=fixed.rate)
			e$migrm[,cidx] <- migpred$M
			e$migrf[,cidx] <- migpred$F
			if(!is.null(migpred$laborM)) {
				e$migrm.labor[,cidx] <- migpred$laborM
				e$migrf.labor[,cidx] <- migpred$laborF
			}
			warns[[country.codes.char[cidx]]] <- c(warns[[country.codes.char[cidx]]], migpred$warns)
		}
		e$popm <- env$totpm.hch[,,variant]
		e$popf <- env$totpf.hch[,,variant]
		if(rebalance) {
			rebalance.migration2groups(e, pop, variant)			
			negatives <- as.character(country.codes[unique(e$negatives)])
			for(country in negatives)
				warns[[country]] <- c(warns[[country]], 'population negative while balancing half child')
		}
		migrationm.hch[,,variant] <- e$migrm + e$migrm.labor
		migrationf.hch[,,variant] <- e$migrf + e$migrf.labor
	}
	for(country in names(warns)) warns[[country]] <- unique(warns[[country]])
	return(list(M=migrationm, F=migrationf, Mhch=migrationm.hch, Fhch=migrationf.hch, warns=warns))
}

do.pop.predict.one.country.no.migration <- function(time, country.name, inpc, kannisto, kantor.pasfr, nr.traj, 
												popM.prev=NULL, popF.prev=NULL, popM.hch.prev=NULL, popF.hch.prev=NULL, 
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
			btm <- btf <- asfert <- pasfert <- matrix(0, nrow=7, ncol=nr.traj)
			deathsm <- deathsf <- matrix(0, nrow=27, ncol=nr.traj, dimnames=list(ages, NULL))
			btm.hch <- btf.hch <- asfert.hch <- pasfert.hch <- matrix(0, nrow=7, ncol=nvariants)
			deathsm.hch <- deathsf.hch <- matrix(0, nrow=27, ncol=nvariants, dimnames=list(ages, NULL))
			mxm <- mxf <- matrix(0, nrow=28, ncol=nr.traj, dimnames=list(mx.ages, NULL))
			mxm.hch <- mxf.hch <- matrix(0, nrow=28, ncol=nvariants, dimnames=list(mx.ages, NULL))
		}
	})
	for(itraj in 1:nr.traj) {
		#asfr <- inpc$PASFR[,time,drop=FALSE]/100.
		pasfr <- kantor.pasfr[[itraj]][,time,drop=FALSE]
		asfr <- pasfr
		for(i in 1:nrow(asfr)) asfr[i,] <- inpc$TFRpred[time,itraj] * asfr[i,]
		#TODO: deal with fixed.mx
		LTres <- modifiedLC(1, kannisto, inpc$e0Mpred[time,itraj], 
									inpc$e0Fpred[time,itraj], verbose=verbose)		
		if(time > 1) { # reset initial population to the one at the previous time step
			pop.ini$M <- popM.prev[,1,itraj] # second dimension is equal 1 (country)
			pop.ini$F <- popF.prev[,1,itraj]
		}
		popres <- PopProjNoMigr(1, pop.ini, LTres, asfr, inpc$SRB[time], country.name=country.name,
									keep.vital.events=keep.vital.events)
		with(res.env, {
			totp[itraj] <- popres$totpop[2]
			totpm[,itraj] <- popres$mpop[,2]
			totpf[,itraj] <- popres$fpop[,2]
			if(keep.vital.events) {
				btm[,itraj] <- popres$mbt
				btf[,itraj] <- popres$fbt
				deathsm[,itraj] <- popres$mdeaths
				deathsf[,itraj] <- popres$fdeaths
				asfert[,itraj] <- asfr
				pasfert[,itraj] <- pasfr*100
				mxm[,itraj] <- LTres$mx[[1]]
				mxf[,itraj] <- LTres$mx[[2]]
			}
		})
	}
	pop.ini <- list(M=inpc$POPm0, F=inpc$POPf0)
	LTres <- modifiedLC(1, kannisto, inpc$e0Mmedian[time], 
								inpc$e0Fmedian[time], verbose=verbose)
	for (variant in 1:nvariants) { # compute the two half child variants
		#asfr <- inpc$PASFR[,time,drop=FALSE]/100.
		pasfr <- kantor.pasfr[[nr.traj+variant]][,time,drop=FALSE]
		asfr <- pasfr
		for(i in 1:nrow(asfr)) asfr[i,] <- inpc$TFRhalfchild[variant,time] * asfr[i,]		
		if(time > 1) {
			pop.ini$M <- popM.hch.prev[,1, variant]
			pop.ini$F <- popF.hch.prev[,1, variant]
		}
		popres <- PopProjNoMigr(1, pop.ini, LTres, asfr, inpc$SRB[time], 
							country.name=country.name, keep.vital.events=keep.vital.events)
		with(res.env, {
			totp.hch[variant] <- popres$totpop[2]
			totpm.hch[,variant] <- popres$mpop[,2]
			totpf.hch[,variant] <- popres$fpop[,2]
			if(keep.vital.events) {
				btm.hch[,variant] <- popres$mbt
				btf.hch[,variant] <- popres$fbt
				deathsm.hch[,variant] <- popres$mdeaths
				deathsf.hch[,variant] <- popres$fdeaths
				asfert.hch[,variant] <- asfr
				pasfert.hch[,variant] <- pasfr*100
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

.get.migration.one.trajectory <- function(use.migration.model, inpc, itraj=NULL, time=NULL,  pop=NULL, ...) {
	if(use.migration.model) #TODO: what should it be for half child variants?
		return(sample.migration.trajectory.from.model(inpc, itraj, time, pop, ...))
	# use migration predictions from inputs
	migpred <- list(M=NULL, F=NULL)
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


project.migration.one.country.one.step <- function(mu, phi, sigma, oldRates, country.code, rmax=NULL, relaxed.bounds=FALSE){
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
	nrates <- length(oldRates)
	oldRate <- oldRates[nrates]
	isGCC <- is.gcc(country.code)
	has.relaxedB <- has.relaxed.bounds(country.code)  
	relaxed <- has.relaxedB && relaxed.bounds
	#fun.max <- paste0("cummulative.max.rate", if(isGCC) "" else ".no.gcc")
	#fun.min <- paste0("cummulative.min.rate", if(isGCC) "" else ".no.gcc")
	fun.max <- paste0("max.multiplicative.pop.change", if(isGCC) "" else ".no.gcc", if(relaxed) ".small" else "")
	fun.min <- paste0("min.multiplicative.pop.change", if(relaxed) ".small" else "")
	#xmin <- .get.rate.limit(oldRates, nrates, fun.min, max, nperiods=6)
	#xmax <- .get.rate.limit(oldRates, nrates, fun.max, min, nperiods=6)
	xmin <- .get.rate.mult.limit(oldRates, nrates, fun.min, max, nperiods=6)
	xmax <- .get.rate.mult.limit(oldRates, nrates, fun.max, min, nperiods=6)
	if(!is.null(rmax)) xmax <- min(xmax, rmax)
	if(has.relaxedB && xmin > xmax) {
		avg <- (xmin + xmax)/2.
		xmin <- avg - 1e-3
		xmax <- avg + 1e-3 
	}
	
  #while(newRate < -0.33 || newRate > 0.665)
  	determ.part <- mu + phi*(oldRate-mu)
  	newRate <- rtruncnorm(n=1,a=xmin-determ.part, b=xmax-determ.part, mean=0, sd=sigma) + determ.part
  	#newRate <- rtruncnorm(n=1, b=xmax-determ.part, mean=0, sd=sigma) + determ.part
  	#newRate <- rnorm(n=1,mean=0, sd=sigma) + determ.part
  	#if (isGCC) stop('')
	# r <- 1
	# while(r < 1000000) {
		# newRate <- mu + phi*(oldRate-mu) + rnorm(n=1,mean=0,sd=sigma)
		# if(is.migrate.within.permissible.range(c(oldRates, newRate), nrates+1, isGCC, fun.min, fun.max)) return(newRate)
		# r <- r+1
	# }
	# stop("Can't simulate new migration rate for country ", country.code)
	return(newRate)
}

.get.rate.limit <- function(rates, n, cumfun, fun, nperiods=12) {
	res <- do.call(cumfun, list(1))
	for(i in 2:min(nperiods,n+1)) {
		s <- sum(rates[(n-i+2):n])
		res <- c(res, do.call(cumfun, list(i)) - s)
	}
	return(do.call(fun, list(res)))
}

.get.rate.mult.limit <- function(rates, n, cumfun, fun, nperiods=6) {
	res <- do.call(cumfun, list(1))
	for(i in 2:min(nperiods,n+1)) {
		p <- prod(1+rates[(n-i+2):n])
		res <- c(res, do.call(cumfun, list(i))/p)
	}
	return(do.call(fun, list(res))-1)
 }

is.gcc <- function(country)
	return(country %in% c(634, 784, 414, 48, 512, 682)) # Qatar, UAE, Kuwait, Bahrain, Oman, SA
	
has.relaxed.bounds <- function(country)
	return(country %in% c(136, 570, 584, 772, 796)) # Cayman Islands, Niue, Marshall Islands, Tokelau, Turks and Caicos Islands

cummulative.max.rate.no.gcc <- function(l)
	switch(l, 0.318445, 0.380350, 0.432385, 0.599260, 0.582940, 0.620895, 
			0.655665, 0.699100, 0.780100, 0.811745, 0.839445, 0.847995)

cummulative.max.rate.no.gcc.with.outliers <- function(l)
	switch(l, 0.538465, 0.628720, 0.927310, 1.195605, 1.266020, 1.586800, 
			1.648005, 1.701605, 1.745850, 1.841215, 2.086650, 2.187455)

cummulative.min.rate.no.gcc <- function(l)
	switch(l, -0.283625, -0.430245, -0.521440, -0.523135, -0.542200, -0.630610, 
			-0.712340, -0.792280, -0.831930, -0.836390, -0.914660, -0.946140)
	
cummulative.max.rate <- function(l)
	switch(l, 0.666770, 1.159820, 1.448865, 1.761840, 1.909545, 2.074115, 
			2.231875, 2.649090, 2.938135, 3.251110, 3.228510, 3.207015)
			
cummulative.min.rate <- function(l)
	switch(l, -0.330165, -0.430245, -0.521440, -0.523135, -0.542200, -0.630610, 
			-0.712340, -0.792280, -0.831930, -0.836390, -0.914660, -0.946140)
			
max.multiplicative.pop.change <- function(l)
	switch(l, 2.087025815, 3.606344662, 4.952147194, 6.992680749, 8.194692131, 9.773177262)
	
max.multiplicative.pop.change.no.gcc <- function(l)
	switch(l, 1.409065294, 1.509547891, 1.606376672, 1.926878976, 1.892288126, 1.932271285)
	
min.multiplicative.pop.change <- function(l)
	switch(l, 0.707796098, 0.62555766, 0.567081044, 0.56605145, 0.559550789, 0.508914659)
	
max.multiplicative.pop.change.no.gcc.small <- function(l)
	switch(l, 1.59473251994369, 2.12674802921728, 2.5940342528116, 3.27346148755602, 4.38710818083359, 6.07190344426719)
	
min.multiplicative.pop.change.small <- function(l)
	switch(l, 0.428948397185301, 0.373042042828842, 0.31526281780426, 0.24164343419914, 0.20427051444469, 0.172963052158649)
	
gcc.upper.threshold <- function(country) {
	switch(as.character(country),
		'634'= 1619, # Qatar
		'784'= 7632, # UAE
		'414'= 1330, # Kuwait
		'48'= 505,  # Bahrain
		'512'= 303, # Oman
		'682'= 3319, # SA
		NA)
}

sample.migration.trajectory.from.model <- function(inpc, itraj=NULL, time=NULL, pop=NULL, 
													popM=NULL, popF=NULL, country.code=NULL, mig.rates=NULL, fixed.rate=NULL) {													
	pars <- inpc$migration.parameters[itraj,]
	land.area <- NA
	if(country.code %in% land_area_wpp2012$country_code)
		land.area <- land_area_wpp2012[land_area_wpp2012$country_code==country.code,'land_area']
	i <- 1
	k <- 1
	zero.constant <- -1e-4
	warns <- c()
	while(TRUE) {
		if(is.null(fixed.rate)) {
			if(all(pars == 0)) rate <- 0
			else rate <- project.migration.one.country.one.step(pars$mu, pars$phi, pars$sigma, 
					c(as.numeric(inpc$migration.rates), mig.rates[1:time]), country.code, 
					rmax=if(pop>0) min(gcc.upper.threshold(country.code)/pop, if(!is.na(land.area)) 44*land.area/pop - 1 else NA, na.rm=TRUE) else NULL,
					relaxed.bounds=time < 4
					# max(colSums(inpc$observed$MIGm + inpc$observed$MIGf)
					)
		} else rate <- fixed.rate
		if(is.na(rate)) stop('Migration rate is NA')
		#mig.count <- ((1+rate)^5 - 1) * pop.prev # instanteneous rate
		mig.count <- rate * pop
		#if(is.gcc(country.code)) { # cap to historical maximum
		#	mig.count <- min(mig.count, max(colSums(inpc$observed$MIGm + inpc$observed$MIGf)))
		#}
		# if(!is.na(land.area) && (pop + mig.count)/land.area > 44) { # check density
			# if(k<100) {
				# k <- k+1
				# next 
			# }
			# warns <- c(warns, 'migration truncated due to high density')
			# mig.count <- 44 * land.area - pop
			# rate <- mig.count/pop
			# #warning('Density too high for ', country.code, ' (', (pop + mig.count)/land.area, ')', immediate.=TRUE)
		# }
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
		if(!is.null(fixed.rate) || rate == 0) break
		if(all(popM[1:21] + migM >= zero.constant) && all(popF[1:21] + migF >= zero.constant))  break # assure positive count
		i <- i+1
		#break
		if(i>1 && ((sum(popM[1:21][abs(msched)>0]) + sum(popF[1:21][abs(fsched)>0]) + mig.count) > zero.constant)
				) { # adjust age schedules
			prev.isneg <- rep(FALSE, 42)
			j <- 1
			sample.new.rate <- FALSE
			while(any(c(popM[1:21] + migM, popF[1:21] + migF) < zero.constant)) { 
				isneg <- prev.isneg | (c(popM[1:21] + migM, popF[1:21] + migF) < zero.constant)
				shifts <- -c(migM + popM[1:21], migF + popF[1:21])
				shifts[!isneg] <- 0
				shifts[prev.isneg] <- 0
				if(sum(shifts)==0) {sample.new.rate <- TRUE; break}
				sched.new <- c(msched, fsched) + shifts/mig.count
				delta <- sum(c(msched, fsched) - sched.new)
				sched.new[!isneg] <- sched.new[!isneg] + delta*abs(sched.new[!isneg])/sum(abs(sched.new[!isneg]))
				msched <- sched.new[1:21]
				fsched <- sched.new[22:length(sched.new)]	
				#tsched.neg <- sum(c(msched,fsched)[isneg])
				#tsched.notneg <- sum(c(msched,fsched)[!isneg])
				#msched[!isneg[1:21]] <- (1-tsched.neg) * msched[!isneg[1:21]]/tsched.notneg
				#fsched[!isneg[22:length(shifts)]] <- (1-tsched.neg) * fsched[!isneg[22:length(shifts)]]/tsched.notneg
				migM <- mig.count * msched
				migF <- mig.count * fsched
				prev.isneg <- isneg
				j <- j+1
				#if(j > 21) stop('Age schedule cannot be adjusted.')
			}
			if(!sample.new.rate) {
				#warns <- c(warns, 'migration age-schedule modified')
				break
			}
		}
		if(i > 1000) {
			warning('Unable to modify age schedule to get positive population for country ', country.code, immediate.=TRUE)
			break
		}
	}
	if(i>2) {
		wints <- c(10,50, seq(100, 10000, by=100))
		wints.ch <- c(3, wints)
		widx <- findInterval(i, wints)+1
		warns <- c(warns, paste('migration rate resampled', paste0(wints.ch[widx], '-', wints.ch[widx+1]), 'times'))	
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

rebalance.migration <- function(e, pop, what='', check.negatives=FALSE) {
	sumpop <- sum(pop)
	which.negative <- c()
	zero.constant <- -1e-4
	for(sex in c('m', 'f')) {
		par <- paste0('migr',sex, what)
		poppar <- paste0('pop', sex)
		for(age in dimnames(e[[par]])[[1]]) {
			i <- 1
			this.sumpop <- sumpop
			this.pop <- pop
			while(i < 100) {		
				dif <- sum(e[[par]][age,])
				dif.countries <- dif/this.sumpop * this.pop
				e[[par]][age,] <- e[[par]][age,] - dif.countries
				if(!check.negatives) break
				wneg <- which(e[[par]][age,] + e[[poppar]][age,] < zero.constant)
				if(length(wneg) > 0) {
					e[[par]][age,] <- pmax(e[[par]][age,], -e[[poppar]][age,])
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

	slabor <- 0
	if(!is.null(e$migrm.labor))
		slabor <- sum(e$migrm.labor)
	rebalance.migration(e, pop) # rest of the world; does not assure positve counts
	
	if(slabor != 0) { # balancing for labor countries
		pop.labor <- pop * (as.integer(names(pop)) %in% labor.countries())
		rebalance.migration(e, pop.labor, what='.labor') # does not assure positve counts
	}
	e$migrm.after <- e$migrm
	e$migrf.after <- e$migrf
	e$migrm.labor.after <- e$migrm.labor
	e$migrf.labor.after <- e$migrf.labor
	after <- sum(e$migrm + e$migrf)
	after.labor <- sum(e$migrm.labor + e$migrf.labor)
	if(abs(after + after.labor - (before+before.labor)) > 1e6) {
				cat('\n', itraj, ' rebalancing migration: total=', after + after.labor - (before+before.labor),
							', labor=', after.labor-before.labor)
				dif <- colSums(e$migrm.after + e$migrf.after + e$migrm.labor.after + e$migrf.labor.after - (e$migrm.before + e$migrf.before + e$migrm.labor.before + e$migrf.labor.before))
				cat('\nMax in ', names(pop)[which.max(abs(dif))], ':', dif[which.max(abs(dif))])
	}
	# set migr to sum in order to deal with negatives
	e$migrm <- e$migrm + e$migrm.labor
	e$migrf <- e$migrf + e$migrf.labor
	# assure positve counts
	rebalance.migration(e, pop, check.negatives=TRUE)
	dif <- abs(e$migrm - (e$migrm.after + e$migrm.labor.after))
	if(max(dif) > 1000) {
		wc <- ceiling(which.max(dif)/21)
		cat('\n', itraj, ' assuring positive pop: dif=', max(dif), ' country ', names(pop)[wc])
		if(max(dif) > 100000) stop('')
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
			asfert <- pasfert <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			btm.hch <- btf.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years, NULL))
			asfert.hch <- pasfert.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
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
			pasfert[1:dim(obs$pasfert)[1],1,] <- obs$pasfert[,dim(obs$pasfert)[2],repi]
			mxm[1:dim(MxKan[[1]]$mx)[1],1,] <- as.matrix(MxKan[[1]]$mx[,dim(MxKan[[1]]$mx)[2]])[repi]
			mxf[1:dim(MxKan[[2]]$mx)[1],1,] <- as.matrix(MxKan[[2]]$mx[,dim(MxKan[[2]]$mx)[2]])[repi]
			migm[1:dim(inpc$observed$MIGm)[1],1,] <- as.matrix(inpc$observed$MIGm[,dim(inpc$observed$MIGm)[2]])[,rep(1,migMntraj)]
			migf[1:dim(inpc$observed$MIGf)[1],1,] <- as.matrix(inpc$observed$MIGf[,dim(inpc$observed$MIGf)[2]])[,rep(1,migFntraj)]
			btm.hch[,1,] <- btm[,1,rep(1,nvariants)]
			btf.hch[,1,] <- btf[,1,rep(1,nvariants)]
			deathsm.hch[,1,] <- deathsm[,1,rep(1,nvariants)]
			deathsf.hch[,1,] <- deathsf[,1,rep(1,nvariants)]
			asfert.hch[,1,] <- asfert[,1,rep(1,nvariants)]
			pasfert.hch[,1,] <- pasfert[,1,rep(1,nvariants)]
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
			for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'pasfert', 'mxm', 'mxf', 'migm', 'migf',
							'btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'pasfert.hch', 'mxm.hch', 'mxf.hch'))
				res.env[[par]][,time+1,] <- envs[[time]][[par]][,cidx,]
		}
	}
	observed <- obs
	with(res.env, {
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch,
			 file = file.path(dest.dir, paste('totpop_country', country, '.rda', sep='')))
		if(keep.vital.events) 
			save(btm, btf, deathsm, deathsf, asfert, pasfert, mxm, mxf, migm, migf,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, pasfert.hch, 
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
	sched.country <- country
	first.year <- FALSE
	# Replace Bahrain and Saudi Arabia with schedule from Qatar 
	if(country %in% c(682, 48)) {
		   sched.country <- 634
	}
	if(country %in% c(28, 531, 364, 376, 462, 562, 630, 662)) { # Antigua and Barbuda, Curacao, Iran, Israel, Maldives, Niger, Puerto Rico, and Saint Lucia
		   sched.country <- 156 # China
		   first.year <- TRUE
	}
	cidxM <- which(inputs$MIGm$country_code==sched.country)
	cidxF <- which(inputs$MIGf$country_code==sched.country)
	col.idx <- which(substr(colnames(inputs$MIGm), 1,4)==as.character(inputs$present.year)):ncol(inputs$MIGm)
	if(first.year) {
		       cix <- rep(which(colnames(inputs$MIGm)=="2010-2015"), length(col.idx))
		       maleVec <- as.matrix(inputs$MIGm[cidxM,cix])
		       femaleVec <- as.matrix(inputs$MIGf[cidxF,cix])
	} else {
	  maleVec <- as.matrix(inputs$MIGm[cidxM,col.idx])
    	  femaleVec <- as.matrix(inputs$MIGf[cidxF,col.idx])
	}
    colnames(maleVec) <- colnames(femaleVec) <- colnames(inputs$MIGm)[col.idx]
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
    	tot <- apply(negM+negF, 2, sum)
    	negM <- t(apply(negM, 1, '/', tot))
    	negF <- t(apply(negF, 1, '/', tot))
    }
    if(country == 818) { # Egypt (special handling when rate positive)
    	negM <- maleArray
    	maleArray[maleArray<0] <- 0
    	negF <- femaleArray
    	femaleArray[femaleArray<0] <- 0
    	tot <- apply(maleArray + femaleArray, 2, sum)
    	maleArray <- t(apply(maleArray, 1, '/', tot))
    	femaleArray <- t(apply(femaleArray, 1, '/', tot))
    }
    if(country %in% c(528, 756)) { # Netherlands, Switzerland  (get Czech schedule for negative schedules)
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
	res <- .C("PopProjNoMigration", as.integer(nproj), 
			srm=LT$sr[[1]], srf=LT$sr[[2]], asfr=as.numeric(as.matrix(asfr)), 
			srb=as.numeric(as.matrix(srb)), returnNothingIfNegative=as.integer(returnIfNegative), 
			popm=popm, popf=popf, totp=totp,
			btagem=as.numeric(btageM), btagef=as.numeric(btageF), 
			deathsm=as.numeric(deathsM), deathsf=as.numeric(deathsF),
			isNegative=as.integer(isNegative)
			)
	vital.events <- list()
	if(keep.vital.events) {
		vital.events$mbt <- res$btagem
		vital.events$fbt <- res$btagef
		vital.events$mdeaths <- res$deathsm
		vital.events$fdeaths <- res$deathsf
	}
	return(c(list(totpop=res$totp, mpop=res$popm, fpop=res$popf), vital.events))
}


