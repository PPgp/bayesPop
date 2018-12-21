if(getRversion() >= "2.15.1") utils::globalVariables(c("land_area_wpp2015", "migration.thresholds"))

do.pop.predict.balance <- function(inp, outdir, nr.traj, ages, pred=NULL, countries=NULL, keep.vital.events=FALSE, function.inputs=NULL, 
									rebalance=TRUE, use.migration.model=TRUE, start.time.index=1, 
									verbose=FALSE, parallel=FALSE, nr.nodes=NULL, migration.settings=NULL, 
									chunk.size=100, reformat.only=FALSE, ...) {
	# if migration.settings is NULL, don't use the migration model
	if (use.migration.model && is.null(migration.settings))
		stop('Argument migration.settings has to be given if use.migration.model is TRUE.')
	countries.idx <- if(is.null(countries)) which(UNlocations$location_type==4) else which(UNlocations$country_code %in% countries)
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
	present.and.proj.years.pop <- present.and.proj.years + 2
	status.for.gui <- paste('out of', nr_project, 'time periods.')
	gui.options <- list()
	inp.to.save <- list()
	# remove big or redundant items from inputs to be saved
	for(item in ls(inp)[!grepl('^migMpred$|^migFpred$|^TFRpred$|^e0Fpred$|^e0Mpred$|^estim.years$|^proj.years$|^wpp.years$', ls(inp))]) 
		inp.to.save[[item]] <- get(item, inp)
	npred <- nr_project
		
	mig.rate.prev <- mig.rate <- NULL
	fixed.mig.rate <- FALSE
	adjust.mig <- FALSE
	if(use.migration.model) {
		if(is.null(migration.settings$posterior) && is.null(migration.settings[['projected.rates']]))
			stop('Argument migration.settings must either have an element "posterior" or an element "projected.rates".')
		if(!is.null(migration.settings$posterior)) {
			inp[['migration.parameters']] <- migration.settings$posterior
			if(!('country_code' %in% colnames(inp[['migration.parameters']])))
				stop('Column "country_code" must be included in the "migration.parameters" dataset.')
		} else { # fixed rates
			inp[['projected.migration.rates']] <- migration.settings[['projected.rates']]
			if(!('country_code' %in% colnames(inp[['projected.migration.rates']])))
				stop('Column "country_code" must be included in the "projected.rates" dataset.')
			fixed.mig.rate <- TRUE
		}
		if(!is.null(migration.settings$ini.rates)) {
			mig.rate.prev <- migration.settings$ini.rates#/5.
		} else {
			mig.rate.prev <- inp$migration.rates[,ncol(inp$migration.rates)]
			names(mig.rate.prev) <- inp$migration.rates$country_code
		}
		inp$migration.rates <- inp$migration.rates[,-ncol(inp$migration.rates)] # remove the present year column as it is in mig.rate.prev
		for(par in c('year.of.schedule', 'adjust.to', 'adjustM.to', 'adjustF.to'))
			if(!is.null(migration.settings[[par]])) inp[[paste0("migration.", par)]] <- migration.settings[[par]]
		if(any(c('adjust.to', 'adjustM.to', 'adjustF.to') %in% names(migration.settings))) adjust.mig <- TRUE
		# remove unwanted (redundant) columns
		remove.cols <- c("country_name", "name")
		for (par in c('migration.parameters', 'projected.migration.rates', 'migration.rates', 
		              'migration.adjust.to', 'migration.adjustM.to', 'migration.adjustF.to')) 
		    if(par %in% names(inp) && any(remove.cols %in% colnames(inp[[par]])))
		        inp[[par]] <- inp[[par]][, -which(colnames(inp[[par]]) %in% remove.cols)]
	} 
	outdir.tmp <- file.path(outdir, '_tmp_')
	if(file.exists(outdir.tmp) && start.time.index==1 && !reformat.only) unlink(outdir.tmp, recursive=TRUE)
	if(start.time.index==1 && !reformat.only) dir.create(outdir.tmp)
	#temp.file.name <- file.path(outdir.tmp, "temp_pop_ve.h5")
	#h5createFile(temp.file.name)
	#h5createGroup(temp.file.name, "pop")
	#h5createGroup(temp.file.name, "ve")
	
	if(start.time.index > 1) { # reload last rates and population
		env.tmp <- new.env()
		#load(file.path(outdir.tmp, paste0('pop_time_', start.time.index-1, '.rda')), envir=env.tmp)
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
		nr.nodes.cntry <- min(nr.nodes, ncountries)
	}
	countries.input <- new.env()
	# Extract the country-specific stuff from the inputs
	if(verbose) cat('\nLoading inputs for ', ncountries, ' countries ')
	if(parallel) {
		if(verbose) cat('(in parallel on ', nr.nodes.cntry, ' nodes).')
		cl <- create.pop.cluster(nr.nodes.cntry, ...)
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
	if(verbose) cat('\nProcessing ', npred, ' time periods for ', ncountries, ' countries for ', nr.traj, ' trajectories ')
	if(length(countries.input) < ncountries) {
		ncountries <- length(countries.input)
		country.codes.char <- ls(countries.input)
		country.codes <- as.integer(country.codes.char)
		countries.idx <- sapply(country.codes, function(x) which(UNlocations[,'country_code'] == x))
		UNnames <- UNlocations[countries.idx,'name']
		if(parallel) nr.nodes.cntry <- min(nr.nodes.cntry, ncountries)
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
	if(!reformat.only) {
	    init.res.env <- function(time, res.env) {
	#h5compression.level <- 0
	with(res.env, {
	    totp <- totpm <- totpf <- migrationm <- migrationf <- migm <- migf <- totp.hch <- totpm.hch <- totpf.hch <- list()
	    #tmpfiles <- list(totp = paste0("totp_", country.codes.char, ".bin")
	    for(country in country.codes.char) {
	        files <- list()
	        for(ind in c("totp", "totpm", "totpf", "migrationm", "migrationf", "migm", "migf", "totp.hch", "totpm.hch", "totpf.hch")) {
	            f <- files[[ind]] <- file.path(outdir.tmp, paste0(ind, "_", time, "_", country, ".bin"))
	            file.create(f)
	        }
	        totp[[country]] <- matter_vec(0, length = nr.traj, filemode = "rb+", paths = files$totp)
	        totpm[[country]] <- matter_mat(0, nrow = 27, ncol = nr.traj, dimnames=list(ages, NULL), filemode = "rb+", paths = files$totpm)
	        totpf[[country]] <- matter_mat(0, nrow = 27, ncol = nr.traj, dimnames=list(ages, NULL), filemode = "rb+", paths = files$totpf)
	        migrationm[[country]] <- matter_mat(0, nrow = 27, ncol = nr.traj, dimnames=list(ages, NULL), filemode = "rb+", paths = files$migrationm)
	        migrationf[[country]] <- matter_mat(0, nrow = 27, ncol = nr.traj, dimnames=list(ages, NULL), filemode = "rb+", paths = files$migrationf)
	        migm[[country]] <- matter_mat(0, nrow = 27, ncol = nr.traj, dimnames=list(ages, NULL), filemode = "rb+", paths = files$migm)
	        migf[[country]] <- matter_mat(0, nrow = 27, ncol = nr.traj, dimnames=list(ages, NULL), filemode = "rb+", paths = files$migf)
	        totp.hch[[country]] <- matter_vec(0, length = nvariants, filemode = "rb+", paths = files$totp.hch)
	        totpm.hch[[country]] <- matter_mat(0, nrow = 27, ncol = nvariants, dimnames=list(ages, NULL), filemode = "rb+", paths = files$totpm.hch)
	        totpf.hch[[country]] <- matter_mat(0, nrow = 27, ncol = nvariants, dimnames=list(ages, NULL), filemode = "rb+", paths = files$totpf.hch)
	    }
	    #totp <- matrix(0, nrow=ncountries, ncol=nr.traj, dimnames=list(country.codes, NULL))
        #totpm <- totpf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
        #migrationm <- migrationf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
        #totp.hch <- matrix(0, nrow=ncountries, ncol=nvariants, dimnames=list(country.codes, NULL))
        #totpm.hch <- totpf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
		if(keep.vital.events) {
		    btm <- btf <- asfert <- pasfert <- array(0, dim=c(7, ncountries, nr.traj), dimnames=list(NULL, country.codes, NULL))
            deathsm <- deathsf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
            btm.hch <- btf.hch <- asfert.hch <- pasfert.hch <- array(0, dim=c(7, ncountries, nvariants), dimnames=list(NULL, country.codes, NULL))
            deathsm.hch <- deathsf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
            mxm <- mxf <- array(0, dim=c(28, ncountries, nr.traj), dimnames=list(mx.ages, country.codes, NULL))
            mxm.hch <- mxf.hch <- array(0, dim=c(28, ncountries, nvariants), dimnames=list(mx.ages, country.codes, NULL))
		}
	})
	debug <- FALSE
	#res.env
	    }
	#migration.thresholds <- NULL
	#migration.thresholds <= get.migration.thresholds(inp$wpp.year)
	assign('migration.thresholds', get.migration.thresholds(inp$wpp.year), envir=.GlobalEnv)
	
	if(parallel) {
		nr.nodes.traj <- min(nr.nodes, nr.traj)
		if(verbose) cat(' (in parallel on ', nr.nodes.traj, ' nodes).')
		cl <- create.pop.cluster(nr.nodes.traj, ...)
		clusterExport(cl, c("nr.traj", "country.codes",  "UNnames", "countries.input", "kannisto", "npasfr", 
								"ages", "nvariants", "keep.vital.events", "verbose",
								"res.env", "npred", "country.codes.char", "kantor.pasfr", 
								"rebalance", "use.migration.model", "fixed.mig.rate", "outdir.tmp", 
								"migration.thresholds"), envir=environment())
	} else if(verbose) cat(' (sequentially).')
	
	adjust.half.child <- function(values, med) {
		res <- values
		idx.eq <- which(values[,1] == values[,2])
		res[idx.eq,] <- med[idx.eq]
		low <- ifelse(res[,1] > med, med, res[,1])
		dif <- res[,1] - low
		res[,1] <- low
		res[,2] <- res[,2] - dif
		high <- ifelse(res[,2] < med, med, res[,2])
		dif <- res[,2] - high
		res[,2] <- high
		res[,1] <- ifelse(res[,1] - dif <= med, res[,1] - dif, med)
		return(res)
	}
	pop.predict.half.child <- function() {
		.ini.pop.res.env(res.env, keep.vital.events)
		for(cidx in 1:ncountries) {
			country <- country.codes[cidx]
			file.name <- file.path(outdir, paste0('totpop_country', country, '.rda'))
			computed.env <- new.env()
			load(file.name, envir=computed.env)
			medmigm <- apply(computed.env$migm, c(1,2), "median")
			medmigf <- apply(computed.env$migf, c(1,2), "median")
			medpopm <- apply(computed.env$totpm, c(1,2), "median")
			medpopf <- apply(computed.env$totpf, c(1,2), "median")
			res.env$totpm.hch[,1,] <- computed.env$totpm[,1,1]
			res.env$totpf.hch[,1,] <- computed.env$totpf[,1,1]
			if(keep.vital.events) {
				file.name.ve <- file.path(outdir, paste0('vital_events_country', country, '.rda'))
				computed.env.ve <- new.env()
				load(file.name.ve, envir=computed.env.ve)
				for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'pasfert', 'mxm','mxf'))  
					res.env[[paste0(par, ".hch")]][,1,] <- computed.env.ve[[par]][,1,1]
			}
			popM.hch.prev <- res.env$totpm.hch[,1,]
			popF.hch.prev <- res.env$totpf.hch[,1,]
			popM.hch.prev[is.na(popM.hch.prev)] <- 0
			popF.hch.prev[is.na(popF.hch.prev)] <- 0
			for(time in 1:npred) {
				nomigpred <- do.pop.predict.one.country.no.migration.half.child(time, 
												UNnames[cidx], countries.input[[country.codes.char[cidx]]], 
												kannisto[[country.codes.char[cidx]]], kantor.pasfr[[country.codes.char[cidx]]], nr.traj,
												popM.hch.prev, popF.hch.prev, ages, 
												nvariants, keep.vital.events=keep.vital.events, verbose=verbose)

				res.env$totpm.hch[,time+1,] <- nomigpred[["totpm.hch"]] + medmigm[,time+1]
				res.env$totpf.hch[,time+1,] <- nomigpred[["totpf.hch"]] + medmigf[,time+1]
				# for ages where both variants are equal, they should be equal to the pop median
				res.env$totpm.hch[,time+1,] <- adjust.half.child(res.env$totpm.hch[,time+1,], medpopm[,time+1])
				res.env$totpf.hch[,time+1,] <- adjust.half.child(res.env$totpf.hch[,time+1,], medpopf[,time+1])

				if(keep.vital.events) {
					for(par in c('btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'pasfert.hch', 'mxm.hch','mxf.hch')) 
						res.env[[par]][,time+1,] <- nomigpred[[par]]	
				}
				popM.hch.prev <- res.env$totpm.hch[,time+1,]
				popF.hch.prev <- res.env$totpf.hch[,time+1,]
			}
			spop <- res.env$totpm.hch + res.env$totpf.hch
			res.env$totp.hch <- apply(spop, c(2,3), sum)
			for(par in c('totp.hch', 'totpm.hch', 'totpf.hch'))
				computed.env[[par]] <- res.env[[par]]
			save(list=ls(computed.env, all.names=TRUE), envir=computed.env, file=file.name)
			if(keep.vital.events) {
				for(par in c('btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'pasfert.hch', 'mxm.hch','mxf.hch')) 
					computed.env.ve[[par]] <- res.env[[par]]
				save(list=ls(computed.env.ve, all.names=TRUE), envir=computed.env.ve, file=file.name.ve)
			}
		}
	}
	
	wrapper.pop.predict.one.country.no.migration <- function(cidx) {
	    if(time > 1) 
            do.pop.predict.one.country.no.migration(time, UNnames[cidx], countries.input[[country.codes.char[cidx]]], 
                                                    kannisto[[country.codes.char[cidx]]], kantor.pasfr[[country.codes.char[cidx]]], nr.traj, 
                                                    popM.prev[[country.codes.char[cidx]]][], 
                                                    popF.prev[[country.codes.char[cidx]]][], ages, 
                                                    keep.vital.events=keep.vital.events, verbose=verbose)
	    else do.pop.predict.one.country.no.migration(time, UNnames[cidx], countries.input[[country.codes.char[cidx]]], 
	                                                 kannisto[[country.codes.char[cidx]]], kantor.pasfr[[country.codes.char[cidx]]], nr.traj,  
	                                                 ages=ages, keep.vital.events=keep.vital.events, verbose=verbose)
	}
	store.no.migration.results <- function(res, code, res.env) {
	    for(par in c('totp', 'totpm', 'totpf')) 
	        res.env[[par]][[code]][] <- res[[par]]
	    
	    if(keep.vital.events) {
	        for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'pasfert', 'mxm', 'mxf')) # 'migm', 'migf',
	            res.env[[par]][,cidx,] <- nomigpred[[cidx]][[par]]
	    }
	}
	res.env <- new.env()
	res.env$mig.rate <- mig.rate
	res.env$keep.vital.events <- keep.vital.events
	res.env$country.codes <- country.codes
	res.env$nr.traj <- nr.traj
	res.env$warns <- list()
	res.env$warns[["_template_"]] <- matrix(0, nrow=get.nr.warns(), ncol=npred)
	for(time in start.time.index:npred) {
	    unblock.gtk.if.needed(paste('finished', time, status.for.gui), gui.options)
	    if(verbose) cat('\nProcessing time period ', time)
		#.ini.pop.res.env(res.env, keep.vital.events)
	    init.res.env(time, res.env)
		if(parallel) {
		    clusterExport(cl, c("time", "mig.rate.prev"), envir=environment())
		    if(time > 1) clusterExport(cl, c("popM.prev", "popF.prev"), envir=environment())
		    nomigpred <- parLapplyLB(cl, 1:ncountries, wrapper.pop.predict.one.country.no.migration)
		    for(cidx in 1:ncountries) store.no.migration.results(cidx)
		} else { # process sequentially
		    for(cidx in 1:ncountries) store.no.migration.results(wrapper.pop.predict.one.country.no.migration(cidx), 
		                                                         country.codes.char[cidx], res.env)
		}
		get.balanced.migration(time, country.codes, countries.input, nr.traj, rebalance, use.migration.model,
								                ages, res.env,  use.fixed.rate=fixed.mig.rate, verbose=verbose)

		# Migration adjustments
		if(adjust.mig) {
		    # 1. adjust
		    #mig.before <- list(M = res.env$migrationm, F = res.env$migrationf)
		    adjust.migration.if.needed(time, present.and.proj.years.pop[time], country.codes, countries.input, res.env)
		    #mig.after <- list(M = res.env$migrationm, F = res.env$migrationf)
		    # 2. re-balance
		    if(rebalance)
		        rebalance.migration.for.all.trajectories(res.env)
		    # adjust migration rates
		    res.env$mig.rate[time + 1,,] <- apply(res.env$migrationm + res.env$migrationf, c(2,3), sum)/apply(res.env$totpm + res.env$totpf, c(2,3), sum)
		    #mig.after.balance <- list(M = res.env$migrationm, F = res.env$migrationf)
        }
		# New population counts
		for(cidx in 1:ncountries) {
		    res.env$totpm[[country.codes.char[cidx]]][] <- res.env$totpm[[country.codes.char[cidx]]][] + res.env$migrationm[[country.codes.char[cidx]]][]
            res.env$totpf[[country.codes.char[cidx]]][] <- res.env$totpf[[country.codes.char[cidx]]][] + res.env$migrationf[[country.codes.char[cidx]]][]
            res.env$totp[[country.codes.char[cidx]]][] <- colSums(res.env$totpm[[country.codes.char[cidx]]][] + res.env$totpf[[country.codes.char[cidx]]][])
            if (any(res.env$totpm[[country.codes.char[cidx]]][] < get.zero.constant()) || any(res.env$totpf[[country.codes.char[cidx]]][] < get.zero.constant())){
                add.pop.warn(country.codes.char[country], time, 5, res.env)  #'Final population negative for some age groups'
            }
		}
		#totpm <- drop(h5read(temp.file.name, "pop/totpm", index = list(time, NULL, NULL, NULL)) + h5read(temp.file.name, "pop/migrationm", index = list(time, NULL, NULL, NULL)))
		#totpf <- drop(h5read(temp.file.name, "pop/totpf", index = list(time, NULL, NULL, NULL)) + h5read(temp.file.name, "pop/migrationf", index = list(time, NULL, NULL, NULL)))
		#totpf <- res.env$totpf + res.env$migrationf
		#h5write(totpm, file=temp.file.name, name="pop/totpm", index = list(time, NULL, NULL, NULL))
		#h5write(totpf, file=temp.file.name, name="pop/totpf", index = list(time, NULL, NULL, NULL))
		#spop <- apply(totpm + totpf, c(2,3), sum)
		#h5write(spop, file=temp.file.name, name="pop/totp", index = list(time, NULL, NULL))
		#res.env$totp <- if(dim(res.env$totp)[1]==1) sum(spop) else apply(spop, c(2,3), sum) # distinction if there is only one country
		popM.prev <- res.env$totpm
		popF.prev <- res.env$totpf
		#if(is.null(dim(popM.prev))) # one country only; dimension dropped
		#	popM.prev <- abind(popM.prev, along=2)
		#if(is.null(dim(popF.prev))) 
	    #	popF.prev <- abind(popF.prev, along=2)
		
		# if (any(totpm < get.zero.constant()) || any(totpf < get.zero.constant())){
		# 	cntries.m <- which(apply(totpm, 2, function(x) any(x < get.zero.constant())))
		# 	cntries.f <- which(apply(totpf, 2, function(x) any(x < get.zero.constant())))
		# 	for(country in unique(c(cntries.m, cntries.f))) {
		# 		neg.times <- unique(which(apply(totpm[,country,], 2, function(x) any(x<0))),
		# 						which(apply(totpf[,country,], 2, function(x) any(x<0))))
		# 		add.pop.warn(country.codes.char[country], neg.times, 5, res.env)  #'Final population negative for some age groups'
		# 	}
		# }
		#with(res.env, {
		#	file.name <- file.path(outdir.tmp, paste0('pop_time_', time, '.rda'))
		#	save(totp, totpm, totpf, mig.rate, migm, migf, file = file.name)
		#	if(keep.vital.events) {
		#		file.name <- file.path(outdir.tmp, paste0('vital_events_time_', time, '.rda'))
		#		save(btm, btf, deathsm, deathsf, asfert, pasfert, mxm, mxf,  file=file.name)
		#	}
		#})	
		#rm(totpm, totpf, spop)
		#gc()	
	} # end time
	
	if(parallel) stopCluster(cl)
    } # end if(!reformat.only)
	
	# TODO: Are these really the same? Do we need both?
	#migm <- h5read(temp.file.name, "pop/migrationm")
	#migf <- h5read(temp.file.name, "pop/migrationf")
	#h5write(migm, file=temp.file.name, name="pop/migm")
	#h5write(migf, file=temp.file.name, name="pop/migf")
	
	if(verbose) cat('\nRe-formatting data ')
	quant.env <- restructure.pop.data.and.compute.quantiles(outdir.tmp, outdir, npred, country.codes.char, ages, nr.traj, countries.input, observed, kannisto, 
					present.and.proj.years, keep.vital.events, 
					#parallel=parallel,  # this can cause memory swapping
					parallel=FALSE, 
					nr.nodes=nr.nodes.cntry, 
					chunk.size=chunk.size, verbose=verbose)
	if(verbose) cat(' done.\n')
	unlink(outdir.tmp, recursive=TRUE)
	#pop.predict.half.child()
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
			   				ages=ages, warnings=res.env$warns), class='bayesPop.prediction')
	prediction.file <- file.path(outdir, 'prediction.rda')
	save(bayesPop.prediction, file=prediction.file)
	cat('\nPrediction stored into', outdir, '\n')
	print.pop.warnings(bayesPop.prediction, which.warns=c(2,3,5))
	return(bayesPop.prediction)
}

get.nr.warns <- function() length(get.pop.all.warns())
get.pop.all.warns <- function() {
	list(
		'Population negative while balancing', # 1
		'Unable to modify age schedule to get positive population', # 2
		'Migration rate resampled more than 500 times', # 3
		'Population negative while balancing half child', # 4
		'Final population negative for some age groups' # 5	
	)
}

get.pop.warn <- function(l) {
	get.pop.all.warns()[[l]]
}

print.pop.warnings <- function(pop.pred, which.warns=NULL) {
	# print accummulated warnings
	warns <- pop.pred$warnings
	all.warns <- get.pop.all.warns()
	if (is.null(which.warns)) which.warns <- 1:length(all.warns)
	cntry.warn <- function(country.warn, warn.code) country.warn[warn.code,]
	for(iwarn in which.warns) {
		cntries.warn <- t(sapply(warns, cntry.warn, iwarn))
		rownames(cntries.warn) <- names(warns)
		cntries.warn <- cntries.warn[apply(cntries.warn, 1, function(x) any(x>0)),,drop=FALSE]
		if(nrow(cntries.warn)>0) {
			cat("\n", get.pop.warn(iwarn), ":\n")
			rownames(cntries.warn) <- sapply(as.integer(rownames(cntries.warn)), function(x) UNlocations[which(UNlocations$country_code==x),'name'])
			colnames(cntries.warn) <- pop.pred$proj.years[-1]
			print(cntries.warn)
		}
	}
	cat("\n")
}

add.pop.warn <- function(country, time, code, env) {
	if(is.null(env$warns[[country]])) 
		env$warns[[country]] <- env$warns[["_template_"]]
	env$warns[[country]][code, time] <- env$warns[[country]][code, time] + 1
}

get.zero.constant <- function() -1e-4

get.balanced.migration <- function(time, country.codes, inputs, nr.traj, rebalance, use.migration.model,
                                   ages, env, use.fixed.rate=FALSE, verbose=FALSE) {
	nr.countries <- length(country.codes)
	e <- new.env()
	labor.codes <- labor.countries()
	mig.ages <- ages[1:21]
	lmig.ages <- length(mig.ages)
	lages <- length(ages)
	e$migrm <- e$migrf <- matrix(NA, ncol=nr.countries, nrow=lmig.ages, dimnames=list(mig.ages, country.codes))
	e$migrm.labor <- e$migrf.labor <- matrix(0, ncol=nr.countries, nrow=lmig.ages, dimnames=list(mig.ages, country.codes))
	
	totpm <- totpf <- matrix(NA, ncol=nr.countries, nrow=lages, dimnames=list(ages, country.codes))
	pop <- rep(NA, nr.countries)
	names(pop) <- country.codes
	country.codes.char <- as.character(country.codes)
	data(land_area_wpp2015)

	for(itraj in 1:nr.traj){
	    totpm[] <- NA
	    totpf[] <- NA
	    pop[] <- NA
	    #pop <- colSums(env$totpm[,,itraj] + env$totpf[,,itraj])
	    #dimnames(totpm) <- dimnames(totpf) <- list(c(ages, seq(105, by = 5, length = 6)), country.codes)
	    #pop <- colSums(totpm + totpf)
	    #names(pop) <- country.codes.char
	    for(cidx in 1:nr.countries) {
	        inpc <- inputs[[country.codes.char[cidx]]]
	        totpm[,cidx] <- env$totpm[[country.codes.char[cidx]]][, itraj]
            totpf[,cidx] <- env$totpf[[country.codes.char[cidx]]][, itraj]
            pop[cidx] <- sum(totpm[,cidx] + totpf[,cidx])
		    migpred <- .get.migration.one.trajectory(use.migration.model, inpc, itraj, time, pop[cidx], 
							popM=totpm[,cidx], popF=totpf[,cidx], country.code=country.codes[cidx], 
							mig.rates=if(!is.null(env$mig.rate)) env$mig.rate[,cidx,itraj] else NULL, 
							fixed.rate=if(use.fixed.rate) inpc$projected.migration.rates[itraj,time] else NULL,
							warn.template=env$warns[["_template_"]])
		    #print(c(list(paste('Country:', country.codes.char[cidx], ', time: ', time, ', traj: ', itraj)), migpred))
		    e$migrm[,cidx] <- migpred$M
		    e$migrf[,cidx] <- migpred$F
		    if(!is.null(migpred$laborM)) {
			    e$migrm.labor[,cidx] <- migpred$laborM
			    e$migrf.labor[,cidx] <- migpred$laborF
		    }
		    env$mig.rate[time+1, cidx, itraj] <- migpred$rate
		    if (!is.null(migpred$warns)) {
			    env$warns[[country.codes.char[cidx]]] <- if(is.null(env$warns[[country.codes.char[cidx]]])) migpred$warns else
															env$warns[[country.codes.char[cidx]]] + migpred$warns			
		    }
	    }
	    e$popm <- totpm
	    e$popf <- totpf
	    # for cases when dimension is dropped (if there is one country)
	    if(is.null(dim(e$popm))) e$popm <- abind(e$popm, along=2)
	    if(is.null(dim(e$popf))) e$popf <- abind(e$popf, along=2)
	
	    if(rebalance) {
		    rebalance.migration2groups(e, pop, itraj)			
		    negatives <- as.character(country.codes[unique(e$negatives)])
		    for(country in negatives) add.pop.warn(country, time, 1, env) # 'Population negative while balancing'
	    }
	    #h5write(e$migrm + e$migrm.labor, env$temp.file.name, "pop/migrationm", index = list(time, 1:lages, NULL, itraj))
	    #h5write(e$migrf + e$migrf.labor, env$temp.file.name, "pop/migrationf", index = list(time, 1:lages, NULL, itraj))
	    for(cidx in 1:nr.countries) {
	        env$migrationm[[country.codes.char[cidx]]][1:lmig.ages,itraj] <- e$migrm[,cidx] + e$migrm.labor[,cidx]
            env$migrationf[[country.codes.char[cidx]]][1:lmig.ages,itraj] <- e$migrf[,cidx] + e$migrf.labor[,cidx]
	    }
	    #env$migrationm[1:lages,,itraj] <- e$migrm + e$migrm.labor
	    #env$migrationf[1:lages,,itraj] <- e$migrf + e$migrf.labor
	}
	rm(list=ls(e), envir=e)
	rm(e)
	return(NULL)
}


do.pop.predict.one.country.no.migration <- function(time, country.name, inpc, kannisto, kantor.pasfr, nr.traj,
												popM.prev=NULL, popF.prev=NULL, ages, keep.vital.events, verbose=FALSE) {
	pop.ini <- list(M=inpc$POPm0, F=inpc$POPf0)
	res.env <- new.env()
	mx.ages <- c(0,1,ages[2:length(ages)])

	with(res.env, {
	    totp <-  rep(NA, nr.traj) 
	    totpm <- totpf <- matrix(NA, nrow=27, ncol=nr.traj, dimnames=list(ages, NULL))
	    if(keep.vital.events) {
	        btm <- btf <- asfert <- pasfert <- matrix(0, nrow=7, ncol=nr.traj)
	        deathsm <- deathsf <- matrix(0, nrow=27, ncol=nr.traj, dimnames=list(ages, NULL))
	        mxm <- mxf <- matrix(0, nrow=28, ncol=nr.traj, dimnames=list(mx.ages, NULL))
	    }
	})
	for(itraj in 1:nr.traj) {
	    # TODO: support for fixed.pasfr & fixed.mx
	    #asfr <- inpc$PASFR[,time,drop=FALSE]/100.
	    pasfr <- kantor.pasfr[[itraj]][,time,drop=FALSE]
	    asfr <- pasfr
	    for(i in 1:nrow(asfr)) asfr[i,] <- inpc$TFRpred[time,itraj] * asfr[i,]
	    LTres <- modifiedLC(1, kannisto, inpc$e0Mpred[time,itraj], 
								inpc$e0Fpred[time,itraj], verbose=verbose)		
	    if(time > 1) { # reset initial population to the one at the previous time step
		    pop.ini$M <- popM.prev[,itraj]
		    pop.ini$F <- popF.prev[,itraj]
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
	res <- as.list(res.env)
	rm(list=ls(res.env), envir=res.env)
	return(res)
}

do.pop.predict.one.country.no.migration.half.child <- function(time, country.name, inpc, kannisto, kantor.pasfr, nr.traj,
												popM.hch.prev=NULL, popF.hch.prev=NULL, 
												ages, nvariants, keep.vital.events, verbose=FALSE) {
	res.env <- new.env()
	mx.ages <- c(0,1,ages[2:length(ages)])
	# TODO: support for fixed.pasfr & fixed.mx
	pop.ini <- list(M=NULL, F=NULL)
	LTres <- modifiedLC(1, kannisto, inpc$e0Mmedian[time], 
								inpc$e0Fmedian[time], verbose=verbose)
    with(res.env, {
		totp.hch <- rep(NA, nvariants)
		totpm.hch <- totpf.hch <- matrix(NA, nrow=27, ncol=nvariants, dimnames=list(ages, NULL))
		if(keep.vital.events) {
			btm.hch <- btf.hch <- asfert.hch <- pasfert.hch <- matrix(0, nrow=7, ncol=nvariants)
			deathsm.hch <- deathsf.hch <- matrix(0, nrow=27, ncol=nvariants, dimnames=list(ages, NULL))
			mxm.hch <- mxf.hch <- matrix(0, nrow=28, ncol=nvariants, dimnames=list(mx.ages, NULL))
		}
	})
	for (variant in 1:nvariants) {
		pasfr <- kantor.pasfr[[nr.traj+variant]][,time,drop=FALSE]
		asfr <- pasfr
		for(i in 1:nrow(asfr)) asfr[i,] <- inpc$TFRhalfchild[variant,time] * asfr[i,]		
		pop.ini$M <- popM.hch.prev[, variant]
		pop.ini$F <- popF.hch.prev[, variant]

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


project.migration.one.country.one.step <- function(mu, phi, sigma, oldRates, country.code, rlim=c(NULL, NULL), relaxed.bounds=FALSE, is.small=FALSE){
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
	fun.max <- paste0("max.multiplicative.pop.change", if(isGCC) "" else ".no.gcc", if(relaxed) ".small" else "")
	fun.min <- paste0("min.multiplicative.pop.change", if(relaxed) ".small" else "")
	xmin <- .get.rate.mult.limit(oldRates, nrates, fun.min, max, nperiods=6)
	xmax <- .get.rate.mult.limit(oldRates, nrates, fun.max, min, nperiods=6)
	if(!is.null(rlim[2])) xmax <- min(xmax, rlim[2])
	if(!is.null(rlim[1])) xmin <- max(xmin, rlim[1])
	if(xmin > xmax) {
		avg <- (xmin + xmax)/2.
		xmin <- avg - 1e-3
		xmax <- avg + 1e-3 
	}
  	determ.part <- mu + phi*(oldRate-mu)
  	newRate <- rtruncnorm(n=1,a=xmin-determ.part, b=xmax-determ.part, mean=0, sd=sigma) + determ.part
	return(newRate)
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
			
migthresh <- function() return(migration.thresholds)	
			
max.multiplicative.pop.change <- function(l) 
	migration.thresholds$cummulative.bounds$upper[l]

max.multiplicative.pop.change.no.gcc <- function(l) 
	migration.thresholds$cummulative.bounds$upper.nogcc[l]
	
max.multiplicative.pop.change.no.gcc.small <- function(l)
	migration.thresholds$cummulative.bounds$upper.small[l]
	
min.multiplicative.pop.change <- function(l) 
	migration.thresholds$cummulative.bounds$lower[l]
	
min.multiplicative.pop.change.small <- function(l)
	migration.thresholds$cummulative.bounds$lower.small[l]

gcc.upper.threshold <- function(country.char) 
	if(country.char %in% colnames(migration.thresholds$absolute.bounds)) migration.thresholds$absolute.bounds[[country.char]] else NA

sample.migration.trajectory.from.model <- function(inpc, itraj=NULL, time=NULL, pop=NULL, 
													popM=NULL, popF=NULL, country.code=NULL, mig.rates=NULL, 
													fixed.rate=NULL, warn.template=NULL) {													
	pars <- inpc$migration.parameters[itraj,]
	land.area <- NA
	if(country.code %in% land_area_wpp2015$country_code)
		land.area <- land_area_wpp2015[land_area_wpp2015$country_code==country.code,'land_area']
	i <- 0
	k <- 1
	zero.constant <- get.zero.constant()
	warns <- NULL
	popM21 <- popM[1:21]
	popF21 <- popF[1:21]
	popMdistr <- popM21/pop
	popFdistr <- popF21/pop
	emigrant.rate.bound <- -0.8
	country.code.char <- as.character(country.code)
	while(i <= 1000) {
		i <- i + 1
		if(is.null(fixed.rate)) {
			if(all(pars == 0)) rate <- 0
			else rate <- project.migration.one.country.one.step(pars$mu, pars$phi, pars$sigma, 
					c(as.numeric(inpc$migration.rates), mig.rates[1:time]), country.code, 
					rlim=c(if(pop>0 && !is.na(land.area)) -(pop - 0.0019*land.area)/pop else NULL, 
						if(pop>0) min(gcc.upper.threshold(country.code.char)/pop, if(!is.na(land.area)) 44*land.area/pop - 1 else NA, na.rm=TRUE) else NULL),
					relaxed.bounds=time < 6, is.small=pop < 200
					# max(colSums(inpc$observed$MIGm + inpc$observed$MIGf)
					)
		} else rate <- fixed.rate
		if(is.na(rate)) stop('Migration rate is NA for country ', country.code, ', time ', time, ', traj ', itraj, 
					'.\npop=', paste(pop, collapse=', '), '\nmig rate=', paste(c(as.numeric(inpc$migration.rates), mig.rates[1:time]), collapse=', '))
		mig.count <- rate * pop
		schedMname <- 'M'
		schedFname <- 'F'
		if(rate < 0 && !is.null(inpc$migration.age.schedule[['Mnegative']])) {
			schedMname <- 'Mnegative'
			schedFname <- 'Fnegative'
		}
		msched <- inpc$migration.age.schedule[[schedMname]][,time]
		fsched <- inpc$migration.age.schedule[[schedFname]][,time]
		 if(is.gcc(country.code)) { 
			modeloutsched <- inpc$migration.age.schedule[['Mnegative']][,time]
			insched <- inpc$migration.age.schedule[['M']][,time] # China
			coefs <- gcc.inrate.coefs()
			Ict <- coefs[1] + coefs[2] * max(rate, 0)
	  		Oct <- Ict - rate
	  		sum.msched <- sum(insched) # proportion of male
	  		insched <- c(insched, fsched)
	  		outmodsched <- c(modeloutsched, fsched)
	  		outsched <- outmodsched * c(popMdistr, popFdistr)
	  		outsched <- outsched/sum(outsched)
	  		inrate <- insched * Ict
	  		outrate <- Oct*outsched
	  		netrate <- inrate - outrate
	  		sched <- netrate/sum(netrate)
			 msched <- sched[1:21]
			 fsched <- sched[22:42]
			 # check depopulation for negative rates
			 isneg <- netrate < 0
			 netmiggcc <- mig.count*sched
			 idepop <- which(isneg & abs(netmiggcc) > abs(emigrant.rate.bound)*c(popM21, popF21))
			 if(length(idepop) > 0) {
			   delta.abs <- abs(netmiggcc) - abs(emigrant.rate.bound)*c(popM21, popF21)
			   delta.abs.sum <- sum(delta.abs[idepop])
			   netmiggcc.mod <- netmiggcc
			   # add delta to depopulated age groups
			   netmiggcc.mod[idepop] <- netmiggcc[idepop] + delta.abs[idepop]
			   # remove the total amount of shifter migration from the remaining age groups
			   netmiggcc.mod[-idepop] <- netmiggcc[-idepop] - delta.abs.sum*abs(sched[-idepop])/sum(abs(sched[-idepop]))
			   sched <- netmiggcc.mod/mig.count # should sum to 1 
			    msched <- sched[1:21]
			    fsched <- sched[22:42]
			    #denom <- sum(msched * popMdistr + fsched * popFdistr)
			 }
			 #msched[isneg[1:21]] <- msched[isneg[1:21]] * popMdistr[isneg[1:21]] / denom
			 #fsched[isneg[22:42]] <- fsched[isneg[22:42]] * popFdistr[isneg[22:42]] / denom
		}
		if(rate < 0 && !is.gcc(country.code)) {
				denom <- sum(msched * popMdistr + fsched * popFdistr)
				denom2 <- c(msched, fsched)/denom
				if(abs(rate) > min((abs(emigrant.rate.bound) / denom2)[denom2 > 0]) && i < 1000) next
				#stop('')
				msched <- msched * popMdistr / denom
				fsched <- fsched * popFdistr / denom
		}
		# age-specific migration counts		
		migM <- mig.count*msched
		migF <- mig.count*fsched
		if(!is.null(fixed.rate) || rate == 0) break
		#if(all(popM21 + migM >= zero.constant) && all(popF21 + migF >= zero.constant))  break # assure positive count
		lower.bounds <- c(popM21 + emigrant.rate.bound * popM21, popF21 + emigrant.rate.bound * popF21)
		if(all(c(popM21 + migM, popF21 + migF) >= lower.bounds))  break
		if(((sum(popM21[abs(msched)>0]) + sum(popF21[abs(fsched)>0]) + mig.count) > sum(lower.bounds[c(abs(msched)>0, abs(fsched)>0)]))
				) { # adjust age schedules
			prev.isneg <- rep(FALSE, 42)
			j <- 1
			sample.new.rate <- FALSE
			#while(any(c(popM21 + migM, popF21 + migF) < zero.constant)) {
			while(any(c(popM21 + migM, popF21 + migF) < lower.bounds)) { 
				isneg <- prev.isneg | (c(popM21 + migM, popF21 + migF) < lower.bounds)
				shifts <- -c(migM + popM21, migF + popF21) + lower.bounds
				#stop('')
				shifts[!isneg] <- 0
				shifts[prev.isneg] <- 0
				if(sum(shifts)==0) {sample.new.rate <- TRUE; break}
				sched.new <- c(msched, fsched) + shifts/mig.count
				delta <- sum(c(msched, fsched) - sched.new)
				sched.new[!isneg] <- sched.new[!isneg] + delta*abs(sched.new[!isneg])/sum(abs(sched.new[!isneg]))
				msched <- sched.new[1:21]
				fsched <- sched.new[22:length(sched.new)]	
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
			if(is.null(warns)) warns <- warn.template
			warns[2,time] <- warns[2,time] + 1 # 'Unable to modify age schedule to get positive population'
			break
		}
	}
	if(i>500) {
		if(is.null(warns)) warns <- warn.template
		warns[3,time] <- warns[3,time] + 1 # 'Migration rate resampled more than 500 times'
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

rebalance.migration <- function(e, pop, what='', check.negatives=FALSE) {
	sumpop <- sum(pop)
	which.negative <- c()
	zero.constant <- get.zero.constant()
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
	rebalance.migration(e, pop) # rest of the world; does not assure positive counts
	
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

rebalance.migration.for.all.trajectories <- function(env) {
    nr.traj <- dim(env$migrationm)[3]
    e <- new.env()
    migbalattr <- list(m = "migrm", f = "migrf")
    migresattr <- list(m = "migrationm", f = "migrationf")
    popbalattr <- list(m = "popm", f = "popf")
    popresattr <- list(m = "totpm", f = "totpf")
    for(itraj in 1:nr.traj) {
        # prepare balancing environment
        for(sex in c("m", "f")) {
            e[[migbalattr[[sex]]]] <- env[[migresattr[[sex]]]][,,itraj]
            e[[popbalattr[[sex]]]] <- env[[popresattr[[sex]]]][,,itraj]
            # for cases when dimension is dropped (if there is one country)
            if(is.null(dim(e[[popbalattr[[sex]]]]))) e[[popbalattr[[sex]]]] <- abind(e[[popbalattr[[sex]]]], along=2)
        }
        pop <- colSums(env$totpm[,,itraj] + env$totpf[,,itraj])
        # rebalance and store results
        rebalance.migration(e, pop, check.negatives=TRUE)
        for(sex in c("m", "f")) env[[migresattr[[sex]]]][,,itraj] <- e[[migbalattr[[sex]]]]
    }
}

labor.countries <- function()
	return(prop.labor.migration()$country_code)
	
prop.labor.migration <- function()
	return(migration.thresholds$labor.proportions)
		
prop.labor.migration.for.country <- function(code){
	props <- prop.labor.migration()
	return(props[props$country_code==code,'proportion'])
}


restructure.pop.data.and.compute.quantiles <- function(source.dir, dest.dir, npred, country.codes, ages, nr.traj, inputs, observed, kannisto, 
									present.and.proj.years, keep.vital.events=FALSE, parallel=FALSE, nr.nodes=NULL, 
									verbose=FALSE, ...){
	
	restructure.pop.data.and.compute.quantiles.one.country <- function(cidx) {		
		country <- country.codes[cidx]
		inpc <- inputs[[country.codes[cidx]]]
		obs <- observed[[country.codes[cidx]]]
		MxKan <- kannisto[[country.codes[cidx]]]
		repi <- rep(1,nr.traj) # index for repeating columns
		res.env <- new.env()
		with(res.env, {
			totp <- matrix(NA, nrow=npredplus1, ncol=nr.traj, dimnames=list(present.and.proj.years.pop, NULL))
			totpm <- totpf <- migm <- migf <- array(NA, dim=c(27, npredplus1, nr.traj), 
										dimnames=list(ages, present.and.proj.years.pop, NULL))
			# values from current year
			totp[1,] <- sum(inpc$POPm0) + sum(inpc$POPf0)
			totpm[1:length(inpc$POPm0),1,] <- as.matrix(inpc$POPm0)[,repi]
			totpf[1:length(inpc$POPf0),1,] <- as.matrix(inpc$POPf0)[,repi]
			migm[1:dim(inpc$observed$MIGm)[1],1,] <- as.matrix(inpc$observed$MIGm[,dim(inpc$observed$MIGm)[2]])[,repi]
			migf[1:dim(inpc$observed$MIGf)[1],1,] <- as.matrix(inpc$observed$MIGf[,dim(inpc$observed$MIGf)[2]])[,repi]
			if(keep.vital.events) {
				btm <- btf <- asfert <- pasfert <- array(0, dim=c(7, npredplus1, nr.traj), 
							dimnames=list(NULL, present.and.proj.years, NULL))
				deathsm <- deathsf <- array(0, dim=c(27, npredplus1, nr.traj), 
							dimnames=list(ages, present.and.proj.years, NULL))
				mxm <- mxf <- array(0, dim=c(28, npredplus1, nr.traj), dimnames=list(mx.ages, present.and.proj.years, NULL))

				# values from current year		
				btm[1:dim(obs$btm)[1],1,] <- obs$btm[,dim(obs$btm)[2],repi]
				btf[1:dim(obs$btf)[1],1,] <- obs$btf[,dim(obs$btf)[2],repi]
				deathsm[1:dim(obs$deathsm)[1],1,] <- obs$deathsm[,dim(obs$deathsm)[2],repi]
				deathsf[1:dim(obs$deathsf)[1],1,] <- obs$deathsf[,dim(obs$deathsf)[2],repi]
				asfert[1:dim(obs$asfert)[1],1,] <- obs$asfert[,dim(obs$asfert)[2],repi]
				pasfert[1:dim(obs$pasfert)[1],1,] <- obs$pasfert[,dim(obs$pasfert)[2],repi]
				mxm[1:dim(MxKan[[1]]$mx)[1],1,] <- as.matrix(MxKan[[1]]$mx[,dim(MxKan[[1]]$mx)[2]])[repi]
				mxf[1:dim(MxKan[[2]]$mx)[1],1,] <- as.matrix(MxKan[[2]]$mx[,dim(MxKan[[2]]$mx)[2]])[repi]
			}
		})
		for(time in 1:npred) {
			for(par in c('totp'))
				res.env[[par]][time + 1,] <- envs[[time]][[country]][[par]][]
			for(par in c('totpm', 'totpf', 'migm', 'migf'))
				res.env[[par]][,time + 1,] <- envs[[time]][[country]][[par]][]
			if(keep.vital.events) {
				for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'pasfert', 'mxm', 'mxf'))
					res.env[[par]][,time+1,] <- envs[[time]][[par]][,cidx,]
			}
		}
		observed <- obs
		file.name <- file.path(dest.dir, paste0('totpop_country', country, '.rda'))
		file.name.ve <- file.path(dest.dir, paste0('vital_events_country', country, '.rda'))
		with(res.env, {
				save(totp, totpm, totpf, migm, migf, file = file.name)
				if(keep.vital.events) 
					save(btm, btf, deathsm, deathsf, asfert, pasfert, mxm, mxf, 
							observed, file=file.name.ve)
		})
		stotpm <- colSums(res.env$totpm, na.rm=TRUE)
		stotpf <- colSums(res.env$totpf, na.rm=TRUE)
	
		quant.env <- new.env()
		with(quant.env, {
			PIs_cqp <- quantM <- quantF <- matrix(NA, nrow=nquant, ncol=npredplus1,
						dimnames=list(quantiles.to.keep, present.and.proj.years.pop))
			quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(nages, nquant, npredplus1),
						dimnames=list(ages, quantiles.to.keep, present.and.proj.years.pop))
			mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(2, npredplus1), 
						dimnames=list(c('mean', 'sd'), present.and.proj.years.pop))

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

	ncountries <- length(country.codes)
	nages <- length(ages)
	envs <- list()
	for(time in 1:npred) {
	    envs[[time]] <- new.env()
	    for(cntry in country.codes) {
	        envs[[time]][[cntry]] <- list()
	        for(par in c('totp'))
	            envs[[time]][[cntry]][[par]] <- matter_vec(paths = file.path(source.dir, paste0(par, '_', time, '_', cntry, '.bin')), 
	                                                       length = nr.traj)
	        for(par in c('totpm', 'totpf', 'migm', 'migf'))
	            envs[[time]][[cntry]][[par]] <- matter_mat(paths = file.path(source.dir, paste0(par, '_', time, '_', cntry, '.bin')), 
	                                                       nrow = nages, ncol = nr.traj)
	        #drop(h5read(source.dir, paste0("pop/", par), index = list(NULL, cidx, NULL)))
	        if(keep.vital.events) {
	            for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'pasfert', 'mxm', 'mxf'))
	                res.env[[par]][,time+1,] <- envs[[time]][[par]][,cidx,]
	        }
	    }
	    #load(file.path(source.dir, paste0('pop_time_', time, '.rda')), envir=envs[[time]])
	#    if(keep.vital.events) 
	#        load(file.path(source.dir, paste0('vital_events_time_', time, '.rda')), envir=envs[[time]])
	}
	#h5f = H5Fopen(source.dir)
	#ages <- as.character(h5f$ages)

	#country.codes <- rownames(envs[[1]]$totp)
	#nr.traj <- h5f$nr.traj
	#ages <- dimnames(envs[[1]]$totpm)[[1]]

	mx.ages <- c(0,1,ages[2:nages])
	npredplus1 <- npred + 1
	present.and.proj.years.pop <- present.and.proj.years + 2
	quantiles.to.keep <- get.quantiles.to.keep()
	nquant <- length(quantiles.to.keep)
	#H5Fclose(h5f)
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
								 "ages", "mx.ages", "npred", "observed", "present.and.proj.years",
								 "keep.vital.events", "verbose"), envir=environment())
		res.list <- parLapplyLB(cl, 1:ncountries, 
							restructure.pop.data.and.compute.quantiles.one.country)
		stopCluster(cl)
	} else { # process sequentially
	    if(verbose) cat('(sequentially) ... \n')
		res.list <- list()				
		for(cidx in 1:ncountries) {
		    if(verbose && interactive()) cat("\r", round(cidx/ncountries * 100), '%')
				res.list[[cidx]] <- restructure.pop.data.and.compute.quantiles.one.country(cidx)
		}
		if(verbose) cat("\n")
	}
		
	for(cidx in 1:ncountries) {
		for(par in c('PIs_cqp', 'mean_sd', 'quantM', 'quantF', 'mean_sdM', 'mean_sdF'))
			quant.env[[par]][cidx,,] <- res.list[[cidx]][[par]]
		for(par in c('quantMage', 'quantFage', 'quantPropMage', 'quantPropFage'))
			quant.env[[par]][cidx,,,] <- res.list[[cidx]][[par]]
	}
	return(quant.env)
}



migration.age.schedule <- function(country, npred, inputs) {
	####################################
	# Set future migration age schedules. 
	# Most countries get the UN future schedule. 
	# In cases in which the UN schedule is messy, use Rogers-Castro (i.e. China) schedule.
	# Special handling of the GCC countries.
	# Original code by Jon Azose
	####################################
	nAgeGroups <- 21
	sched.country <- country
	first.year <- FALSE # indicates if the schedule is taken from one time period only (defined by first.year.period)
	first.year.period <- paste(inputs$proj.years[1]-3, inputs$proj.years[1]+2, sep='-')
	scale.to.totals <- NULL
	mig.settings <- inputs[['MIGtype']][inputs[['MIGtype']]$country_code==country,]
	# Country can take a schedule from a different country, e.g. China
	if(mig.settings[,"MigAgeSchedule"] > 0) {
		   sched.country <- mig.settings[,"MigAgeSchedule"]
		   first.year <- TRUE
	}
	# Should the Male/Female ratio be kept or set equal. E.g. China schedule has larger migration for female, so rescale
	if(mig.settings[,"MigAgeEqualMFratio"] == 1) 
		scale.to.totals <- list(M=0.5, F=0.5) 

	cidxM <- which(inputs$MIGm$country_code==sched.country)
	cidxF <- which(inputs$MIGf$country_code==sched.country)
	col.idx <- which(colnames(inputs$MIGm)==first.year.period):ncol(inputs$MIGm)
	if(is.gcc(country)) {
		cidxM.neg <- which(inputs$MIGm$country_code==country)
		cidxF.neg <- which(inputs$MIGf$country_code==country)
		first.year.neg <- FALSE
	}
	if(!is.null(inputs$migration.year.of.schedule)) { 
		first.year.period <- inputs$migration.year.of.schedule
		first.year <- TRUE
		if(is.gcc(country)) first.year.neg <- TRUE
	}
	get.schedule <- function(fyear, idxM, idxF, scale.to.totals=NULL) {
		maleArr <- matrix(0, nrow=nAgeGroups, ncol=npred)
		femaleArr <- matrix(0, nrow=nAgeGroups, ncol=npred)
		if(fyear) { # take one year as the age-schedule for all future years
			cix <- rep(which(colnames(inputs$MIGm)==first.year.period), length(col.idx))
			if(length(cix)==0) stop("Time period ", first.year.period, " not found in the migration data.")
			maleV <- as.matrix(inputs$MIGm[idxM,cix])
			femaleV <- as.matrix(inputs$MIGf[idxF,cix])
		} else {# take all years starting from present year
		  	maleV <- as.matrix(inputs$MIGm[idxM,col.idx])
	    	femaleV <- as.matrix(inputs$MIGf[idxF,col.idx])
		}
		total.mig.positive <- (colSums(maleV) + colSums(femaleV)) > 0
		if(!is.null(scale.to.totals)) {
			maleV <- t(scale.to.totals$M * apply(maleV, 1, '/', colSums(maleV)))
			femaleV <- t(scale.to.totals$F * apply(femaleV, 1, '/', colSums(femaleV)))
		}
	    colnames(maleV) <- colnames(femaleV) <- colnames(inputs$MIGm)[col.idx]
	    tot <- colSums(maleV+femaleV)
	    if(any(tot == 0)) {
			#Pull a model schedule to use in scenarios where the projection is 0
			#Use China's 2010-2015 data as the model
			modelmaleVec <- inputs$MIGm[inputs$MIGm$country_code==156, first.year.period]
			modelfemaleVec <- inputs$MIGf[inputs$MIGf$country_code==156, first.year.period]
			modelscale <- if(is.null(scale.to.totals)) list(M=0.5, F=0.5) else scale.to.totals	
			modelmaleVec <- modelscale$M[1] * modelmaleVec/sum(modelmaleVec)
			modelfemaleVec <- modelscale$F[1] * modelfemaleVec/sum(modelfemaleVec)
			modeltot <- sum(modelmaleVec+modelfemaleVec)
			modelM <- modelmaleVec/modeltot
			modelF <- modelfemaleVec/modeltot
		}
		non.zero <- tot != 0
		if(any(non.zero)) {
			maleArr[,which(non.zero)] <- t(apply(maleV[,which(non.zero), drop=FALSE], 1, '/', tot[which(non.zero)]))
	    	femaleArr[,which(non.zero)] <- t(apply(femaleV[,which(non.zero), drop=FALSE], 1, '/', tot[which(non.zero)]))
	    }
	    if(any(!non.zero)) {
	    	maleArr[,which(!non.zero)] <- matrix(modelM, nrow=nAgeGroups, ncol=sum(!non.zero))
	    	femaleArr[,which(!non.zero)] <- matrix(modelF, nrow=nAgeGroups, ncol=sum(!non.zero))
	    }
	    return(list(maleArr, femaleArr, total.mig.positive))
	}

	if(is.gcc(country)) { # for GCC countries keep the original ratio of male to female (for positive net migration)
    	unscheds <- get.schedule(first.year.neg, cidxM.neg, cidxF.neg)
    	scale <- c(colSums(unscheds[[1]]), colSums(unscheds[[2]]))
    	scale.to.totals <- list(M=scale[1:npred], F=scale[(npred+1):length(scale)])
    }
	scheds <- get.schedule(first.year, cidxM, cidxF, scale.to.totals=scale.to.totals)
	maleArray <- scheds[[1]]
	femaleArray <- scheds[[2]]
	
    # special handling for negative and positive migration rates
    negM <- negF <- NULL

    if(is.gcc(country)) {
    	negF <- femaleArray # female gets China schedule; should be correctly scaled
    	# male - use model out-migration schedule (derived from SA)
    	negMvec <- gcc.model.outschedule()
    	negMvec <- negMvec/sum(negMvec)
    	negM <- matrix(negMvec, nrow=nAgeGroups, ncol=npred)*matrix(scale[1:npred], ncol=npred, nrow=nAgeGroups, byrow=TRUE) # scale
    }
    # For some counties like Egypt, if positive migration rate, set negative schedules to zero, since they would mean out-migration
    if(mig.settings[,"MigAgeZeroNeg"] == 1) {
    	negM <- maleArray
    	maleArray[maleArray<0] <- 0
    	negF <- femaleArray
    	femaleArray[femaleArray<0] <- 0
    	# rescale the rest
    	tot <- apply(maleArray + femaleArray, 2, sum)
    	maleArray <- t(apply(maleArray, 1, '/', tot))
    	femaleArray <- t(apply(femaleArray, 1, '/', tot))
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
			srb=as.numeric(as.matrix(srb)), 
			mxm=LT$mx[[1]], mxf=LT$mx[[2]],
			returnNothingIfNegative=as.integer(returnIfNegative), 
			popm=popm, popf=popf, totp=totp,
			btagem=as.numeric(btageM), btagef=as.numeric(btageF), 
			deathsm=as.numeric(deathsM), deathsf=as.numeric(deathsF)
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

get.migration.thresholds <- function(wpp.year=2017, nperiods=6) {
	# Setting various thresholds used in the migration model
	
	# Cummulative thresholds
	do.call("data", list(paste0("migration_rates_wpp", wpp.year)))
	rates.all <- get(paste0("migration_rates_wpp", wpp.year))
	rates <- rates.all[,3:ncol(rates.all)]
	rownames(rates) <- rates.all$country_code
	# GCC plus Western Sahara & Djibouti
	gcc.plus <- rates.all$country_code[is.gcc(rates.all$country_code) | rates.all$country_code %in% c(732, 262)] 

	rMat <- 1 + rates
	tu <- apply(rMat, 1, max)
	tl <- apply(rMat, 1, min)
	for (i in 2:nperiods) {
		p <- 0*rates[,1:(ncol(rates)-i+1)] + 1 # init with 1
		for(j in 1:i) 	
			p <- p * rMat[,j:(ncol(rates)-i+j)]
		tu <- cbind(tu, apply(p, 1, max))
		tl <- cbind(tl, apply(p, 1, min))
	}
	upper.bounds <- apply(tu, 2, max)
	upper.bounds.nogcc <- apply(tu[!rownames(tu) %in% gcc.plus,], 2, max)
	lower.bounds <- apply(tl, 2, min)

	df <- data.frame(upper=upper.bounds, upper.nogcc=upper.bounds.nogcc, lower=lower.bounds)
	# need to update this
	df$upper.small <- c(1.59473251994369, 2.12674802921728, 2.5940342528116, 3.27346148755602, 4.38710818083359, 6.07190344426719)
	df$lower.small <- c(0.428948397185301, 0.373042042828842, 0.31526281780426, 0.24164343419914, 0.20427051444469, 0.172963052158649)
	cumbounds <- df
	
	# Absolute upper bounds for GCC
	absbounds <- data.frame(
	        '634'= 1473, # Qatar
			'784'= 5674, # UAE
			'414'= 1649, # Kuwait
			'48'= 290,   # Bahrain
			'512'= 2253, # Oman
			'682'= 2716, # SA
		check.names=FALSE)
		
	labor.props <- data.frame(
		country_code=c(48, 50, 818, 356, 360, 414, 512, 586, 608, 634, 682, 784),
		proportion=c(0.981495, 0.708461, 0.50817, 0.336513, 0.2947, 1.71062, 1.083171, 0.186738, 0.133457, 0.85757, 1.035762, 0.959949)
				)
	labor.props <- merge(labor.props, UNlocations[,c('country_code', 'name')], sort=FALSE)
	return(list(cummulative.bounds=cumbounds, absolute.bounds=absbounds, labor.proportions=labor.props))
}

gcc.model.outschedule <- function()
	c(0.049662, 0.020928, 0.012675, 0.031806, 0.042943, 0.047495, 0.138648, 0.148062, 0.142722, 0.05879, 
		0.024265, 0.016177, 0.010399, 0.006933, 0.004622, 0.002889, 0.002889, 0, 0, 0, 0)
		
gcc.inrate.coefs <- function() c(0.07008, 1.05246)