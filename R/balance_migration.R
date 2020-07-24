if(getRversion() >= "2.15.1") utils::globalVariables(c("land_area_wpp2019", "migration.thresholds"))

do.pop.predict.balance <- function(inp, outdir, nr.traj, ages, pred=NULL, countries=NULL, keep.vital.events=FALSE, 
                                   fixed.mx=FALSE, fixed.pasfr=FALSE, function.inputs=NULL, 
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
	world.pop.distr.ini <- NULL
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
		pop0all <- data.table(inp$POPm0)[country_code %in% country.codes][, 1:3, with = FALSE]
		colnames(pop0all)[3] <- "pop0"
		pop0all$pop0 <- pop0all$pop0 + data.table(inp$POPf0)[country_code %in% country.codes][, 3, with = FALSE]
		pop0all <- pop0all[, .(totpop = sum(pop0)), by = "age"]
		world.pop.distr.ini <- pop0all$totpop/sum(pop0all$totpop)
		# save a Roger's Castro curve (taken from China)
		migRC <- age.specific.migration(wpp.year=inp$wpp.year, countries=156, verbose=FALSE)
		mig.rogers.castro <- migRC$male[,"2020-2025"]/sum(migRC$male[,"2020-2025"])
	} 
	outdir.tmp <- file.path(outdir, '_tmp_')
	if(file.exists(outdir.tmp) && start.time.index==1 && !reformat.only) unlink(outdir.tmp, recursive=TRUE)
	if(start.time.index==1 && !reformat.only) dir.create(outdir.tmp)
	
	if(start.time.index > 1) { # reload last rates and population
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
	}
	countries.input <- new.env()
	# Extract the country-specific stuff from the inputs
	if(verbose) cat('\nLoading inputs for ', ncountries, ' countries ')
	if(parallel) {
		if(verbose) cat('(in parallel on ', nr.nodes, ' nodes).')
		cl <- create.pop.cluster(min(nr.nodes, ncountries), ...)
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
	}
	if(!is.null(mig.rate.prev)) {
		mig.rate.prev <- mig.rate.prev[country.codes.char]
		mig.rate.prev <- matrix(mig.rate.prev, nrow=length(mig.rate.prev), ncol=nr.traj)
		if(is.null(mig.rate)) 
		    mig.rate <- array(NA, c(npred+1, ncountries, nr.traj), dimnames=list(NULL, country.codes, NULL))
        mig.rate[start.time.index,,] <- mig.rate.prev
	}
	kannisto <- surv <- mortcast.args <- list()
	observed <- new.env()
	kantor.pasfr <- new.env()
	# compute Kannisto, Lee-Carter parameters and PASFR
	for(country in country.codes.char) {
	    if(!fixed.mx) {
		    kann <- runKannisto(countries.input[[country]], inp$start.year, lc.for.all = inp$lc.for.all, npred = npred)
		    mortcast.args[[country]] <- .prepare.for.mortality.projection(pattern = countries.input[[country]]$MXpattern, mxKan = kann, 
		                                                       hiv.params = inpc$HIVparams, lc.for.all = inp$lc.for.all)
	    } else {
	        kann.pred <- runKannisto.noLC(countries.input[[country]])
	        surv[[country]] <- survival.fromLT(npred, kann.pred, verbose = verbose)
	    }
		if(keep.vital.events) 
			observed[[country]] <- compute.observedVE(countries.input[[country]], inp$pop.matrix, 
										countries.input[[country]]$MIGtype, kann, 
										as.integer(country), inp$estim.years)
		tfr.med <- apply(countries.input[[country]]$TFRpred, 1, median)[nrow(countries.input[[country]]$TFRpred)]
		kantor.pasfr[[country]] <- list()
		for(itraj in 1:nr.traj)	
			kantor.pasfr[[country]][[itraj]] <- if(fixed.pasfr) countries.input[[country]]$PASFR/100. else
			                                           kantorova.pasfr(c(countries.input[[country]]$observed$TFRpred, countries.input[[country]]$TFRpred[,itraj]), 
										countries.input[[country]], norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)
		
		for (variant in 1:nvariants) 
			kantor.pasfr[[country]][[variant + nr.traj]] <- if(fixed.pasfr) countries.input[[country]]$PASFR/100. else                      
			                                                       kantorova.pasfr(c(countries.input[[country]]$observed$TFRpred, 
			                            countries.input[[country]]$TFRhalfchild[variant,]), 
										countries.input[[country]], norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)					
	}
	
    if(!reformat.only) {
    # prepare environment for storing results
    res.env <- new.env()
    with(res.env, {
        totp <- matrix(0, nrow=ncountries, ncol=nr.traj, dimnames=list(country.codes, NULL))
        totpm <- totpf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
        migm <- migf <- array(0, dim=c(27, ncountries, nr.traj), dimnames=list(ages, country.codes, NULL))
        warns <- list()
    })
    res.env$mig.rate <- mig.rate
    
    # prepare working environment for storing population for one trajectory
    create.work.env <- function() {
        work.env <- new.env()
        mapply(assign, global.objects, mget(global.objects, inherits = TRUE), MoreArgs = list(envir = work.env))
	    with(work.env, {
            totpm <- totpf <- array(0, dim=c(27, ncountries), dimnames=list(ages, country.codes))
            migm <- migf <- array(0, dim=c(27, ncountries), dimnames=list(ages, country.codes))
		    if(keep.vital.events) {
		        btm <- btf <- asfert <- pasfert <- array(0, dim=c(7, ncountries), dimnames=list(NULL, country.codes))
                deathsm <- deathsf <- array(0, dim=c(27, ncountries), dimnames=list(ages, country.codes))
                #btm.hch <- btf.hch <- asfert.hch <- pasfert.hch <- array(0, dim=c(7, ncountries, nvariants), dimnames=list(NULL, country.codes, NULL))
                #deathsm.hch <- deathsf.hch <- array(0, dim=c(27, ncountries, nvariants), dimnames=list(ages, country.codes, NULL))
                mxm <- mxf <- array(0, dim=c(28, ncountries), dimnames=list(mx.ages, country.codes))
                #mxm.hch <- mxf.hch <- array(0, dim=c(28, ncountries, nvariants), dimnames=list(mx.ages, country.codes, NULL))
		    }
		    warns <- list()
		    warns[["_template_"]] <- matrix(0, nrow=get.nr.warns(), ncol=npred)
		    mig.rate <- array(NA, c(npred+1, ncountries), dimnames=list(NULL, country.codes))
	    })
        return(work.env)
    }
    create.mig.env <- function() {
        work.env <- new.env()
        lages <- length(ages21)
        with(work.env, {
            migrm <- migrf <- matrix(NA, ncol=ncountries, nrow=lages, dimnames=list(ages21, country.codes))
            migrm.labor <- migrf.labor <- matrix(0, ncol=ncountries, nrow=lages, dimnames=list(ages21, country.codes))
        })
        return(work.env)
    }
	debug <- FALSE
	
	migration.thresholds <- get.migration.thresholds(inp$wpp.year)
	
	# For parallel processing we need to make these objects global, so that 
	# they do not need to be send back and fort for each trajectory.
	# To have a consistent sequential application, we make these global
	# even if not running in parallel.
	global.objects <- c("nr.traj", "country.codes",  "UNnames", "countries.input", 
	                    "mortcast.args", "surv", "npasfr", "fixed.mx",
	                    "ages", "ages21", "mx.ages", "ncountries",
	                    "nvariants", "keep.vital.events", "verbose",
	                    "npred", "country.codes.char", "kantor.pasfr", 
	                    "rebalance", "use.migration.model", "fixed.mig.rate", "outdir.tmp", 
	                    "migration.thresholds", "adjust.mig", 
	                    "world.pop.distr.ini", "mig.rogers.castro")
	assign("migration.thresholds", migration.thresholds, envir=.GlobalEnv)
	work.env <- create.work.env()
	#mapply(assign, global.objects, mget(global.objects, inherits = TRUE), MoreArgs = list(envir=settings.env))
	if(parallel) {
		nr.nodes.traj <- min(nr.nodes, nr.traj)
		if(verbose) cat(' (in parallel on ', nr.nodes.traj, ' nodes).')
		cl <- create.pop.cluster(nr.nodes.traj, ...)
		#clusterExport(cl, settings.env, envir=environment())
		clusterExport(cl, c(global.objects, "global.objects"), envir=environment())
		clusterExport(cl, c("create.work.env", "create.mig.env",
		                    ".ini.pop.res.env", "cleanup.env"), envir=environment())
		clusterEvalQ(cl, {
		    data(land_area_wpp2019)
		    #library(pryr)
		    work.env <- create.work.env()
		    work.mig.env <- create.mig.env()
		})
	} else {
	    if(verbose) cat(' (sequentially).')
	    work.mig.env <- create.mig.env()
	}
	data(land_area_wpp2019)
	
	predict.1time.period.1trajectory <- function(itraj, time) {
	    #gc()
	    #memch1a <- mem_change({
	    with(work.mig.env, {
	        migrm[] <- migrf[] <- migrm.labor[] <- migrf.labor[] <- 0
	    })
	    .ini.pop.res.env(work.env, keep.vital.events)
	    work.env$mig.rate[] <- mig.rate[,,itraj]
	    #})
	    #memch1d <- mem_change({
	    if(time > 1) {
	        work.env$popM.prev <- popM.prev[,,itraj]
	        work.env$popF.prev <- popF.prev[,,itraj]
	        if(is.null(dim(work.env$popM.prev))) # one country only; dimension dropped
	            work.env$popM.prev <- abind(work.env$popM.prev, along=2)
	        if(is.null(dim(work.env$popF.prev))) 
	            work.env$popF.prev <- abind(work.env$popF.prev, along=2)
	        # recompute world pop age distribution
	        popall <- rowSums(work.env$popM.prev + work.env$popF.prev)
	        popall <- c(popall[1:20], sum(popall[21:length(popall)]))
	        work.env$world.pop.distr <- popall/sum(popall)
	    } else work.env$world.pop.distr <- work.env$world.pop.distr.ini
	    #})
        #memch2 <- mem_change({
	    balanced.migration.1traj(time, itraj, work.env = work.mig.env, res.env = work.env)
        #})
        #memch3 <- mem_change({
	    if(keep.vital.events) {
	        # save vital events by trajectories so that we don't need to send it back to master 
	        #file.name <- file.path(outdir.tmp, paste0('vital_events_time_traj_', time, '_', itraj, '.rda'))
	        file.name <- file.path(outdir.tmp, paste0('vital_events_traj_', itraj, '.rda'))
	        items.to.save <- c("btm", "btf", "deathsm", "deathsf", "asfert", "pasfert", "mxm", "mxf")
	        store.env <- new.env()
	        if(time == 1) { # add time dimension
	            for(item in items.to.save)
	                store.env[[item]] <- abind(work.env[[item]], rev.along = 3)
	        } else {
	            load(file.name, envir = store.env)
	            for(item in items.to.save)
	                store.env[[item]] <- abind(store.env[[item]], work.env[[item]], along = 1)
	        }
	        save(list = items.to.save, envir = store.env, file = file.name)
	        #with(work.env, {
	        #    save(btm, btf, deathsm, deathsf, asfert, pasfert, mxm, mxf,  file=file.name)
	        #})
	        #print(file.name)
	        cleanup.env(store.env)
	    }
	    res <- list()
	    #for(item in c('warns'))
	    #    res[[item]] <- work.env[[item]]
	    for(item in c('totpm', 'totpf', 'migm', 'migf', 'mig.rate', 'warns'))
	        res[[item]] <- work.env[[item]]
	    #})
	    #if(itraj %in% c(1, 2, nr.traj, nr.traj - 1)) {
	    #    cat("\nTime: ", time, " traj: ", itraj, "\n")
	    #    print(c(memch1a, memch1d, memch2, memch3))
	    #    print(mem_used())
	    #}
	    return(res)
	}
	
	for(time in start.time.index:npred) {
	    unblock.gtk.if.needed(paste('finished', time, status.for.gui), gui.options)
	    if(verbose) cat('\nProcessing time period ', time)
		.ini.pop.res.env(res.env, vital.events = FALSE)
		if(parallel) {
		    clusterExport(cl, c("time", "mig.rate"), envir=environment())
		    if(time > 1) clusterExport(cl, c("popM.prev", "popF.prev"), envir=environment())
		    thispred <- clusterApplyLB(cl, 1:nr.traj, predict.1time.period.1trajectory, time = time)
		} else { # process sequentially
		    thispred <- list()
		    for(itraj in 1:nr.traj) {
		        thispred[[itraj]] <- predict.1time.period.1trajectory(itraj, time)
		    }
		}
		#memch1 <- mem_change({
		# collect results
		for(itraj in 1:nr.traj) { 
			for(par in c('totpm', 'totpf', 'migm', 'migf'))
				res.env[[par]][,,itraj] <- thispred[[itraj]][[par]]
			res.env$mig.rate[,,itraj] <- thispred[[itraj]]$mig.rate
			for(country in setdiff(names(thispred[[itraj]]$warns), "_template_")) {
			    if(is.null(res.env$warns[[country]]))
			        res.env$warns[[country]] <- thispred[[itraj]]$warns[["_template_"]]
			    res.env$warns[[country]] <- res.env$warns[[country]] + thispred[[itraj]]$warns[[country]]
			}
		}
		rm(thispred)
		#})
		# Migration adjustments
		if(adjust.mig) {
		    # population without migration
		    res.env$totpm <- res.env$totpm - res.env$migm
		    res.env$totpf <- res.env$totpf - res.env$migf
		    # 1. adjust
		    #mig.before <- list(M = res.env$migm, F = res.env$migf)
		    adjust.migration.if.needed(time, present.and.proj.years.pop[time], country.codes, countries.input, res.env)
		    #mig.after <- list(M = res.env$migrationm, F = res.env$migrationf)
		    # 2. re-balance
		    if(rebalance)
		        rebalance.migration.for.all.trajectories(res.env)
		    # adjust migration rates
		    res.env$mig.rate[time + 1,,] <- apply(res.env$migm + res.env$migf, c(2,3), sum)/apply(res.env$totpm + res.env$totpf, c(2,3), sum)
		    #mig.after.balance <- list(M = res.env$migrationm, F = res.env$migrationf)
		    # New population counts
		    res.env$totpm <- res.env$totpm + res.env$migm
		    res.env$totpf <- res.env$totpf + res.env$migf
		}
		#memch2 <- mem_change({
		spop <- res.env$totpm + res.env$totpf
		res.env$totp <- if(dim(res.env$totp)[1]==1) sum(spop) else apply(spop, c(2,3), sum) # distinction if there is only one country
		popM.prev <- res.env$totpm
		popF.prev <- res.env$totpf
		if(is.null(dim(popM.prev))) # one country only; dimension dropped
			popM.prev <- abind(popM.prev, along=2)
		if(is.null(dim(popF.prev))) 
			popF.prev <- abind(popF.prev, along=2)
        mig.rate <- res.env$mig.rate
		#})
		#memch3 <- mem_change({
		if (any(res.env$totpm < get.zero.constant()) || any(res.env$totpf < get.zero.constant())){
			cntries.m <- which(apply(res.env$totpm, 2, function(x) any(x < get.zero.constant())))
			cntries.f <- which(apply(res.env$totpf, 2, function(x) any(x < get.zero.constant())))
			for(country in unique(c(cntries.m, cntries.f))) {
				neg.times <- unique(which(apply(res.env$totpm[,country,], 2, function(x) any(x<0))),
								which(apply(res.env$totpf[,country,], 2, function(x) any(x<0))))
				add.pop.warn(country.codes.char[country], neg.times, 5, res.env)  #'Final population negative for some age groups'
			}
		}#})
		#memch4 <- mem_change({
		with(res.env, {
			file.name <- file.path(outdir.tmp, paste0('pop_time_', time, '.rda'))
			save(totp, totpm, totpf, mig.rate, migm, migf, file = file.name)
		})#})
		#cat("\nMain loop time: ", time, "\n")
		#print(c(memch1, memch2, memch3, memch4))
		#print(mem_used())
		#gc()
	} # end time
	
	if(parallel) stopCluster(cl)
	else {
	    cleanup.env(work.mig.env)
	}
	rm(migration.thresholds, envir = .GlobalEnv)
    } # end if(!reformat.only)
	
	if(verbose) cat('\nRe-formatting data ')
	quant.env <- restructure.pop.data.and.compute.quantiles(outdir.tmp, outdir, npred, nr.traj, countries.input, observed, 
					present.and.proj.years, keep.vital.events, 
					#parallel=parallel,  # this can cause memory swapping
					parallel=FALSE, 
					nr.nodes=nr.nodes.cntry, 
					chunk.size=chunk.size, verbose=verbose)
	if(verbose) cat(' done.\n')
	unlink(outdir.tmp, recursive=TRUE)
	pop.predict.half.child(res.env, outdir, work.env)
	cleanup.env(res.env)
	cleanup.env(work.env)
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
	# cleanup 
	cleanup.env(countries.input)
	cleanup.env(kantor.pasfr)
	cleanup.env(observed)
	
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

.ini.pop.res.env <- function(e, vital.events=FALSE, debug = FALSE) {
    if(debug) {
        print("before")
        print(c(address(e), refs(e)))
        if(exists('totpm', envir = e)) print(c("totpm: ", inspect("totpm", env = e)))
        if(exists('btm', envir = e)) print(c("btm: ", inspect("btm", env = e)))
    }
    for(item in c('totp', 'totpm', 'totpf', 'migm', 'migf')) 
        if(exists(item, envir = e)) e[[item]][] <- 0
    if(vital.events) {
        for(item in c('btm', 'btf', ' deathsm', ' deathsf', 'mxm', 'mxf'))
            if(exists(item, envir = e)) e[[item]][] <- 0
    }
        if(debug) {
            print("after")
            print(c(address(e), refs(e)))
            if(exists('totpm', envir = e)) print(c("totpm: ", inspect("totpm", env = e)))
            if(exists('btm', envir = e)) print(c("btm: ", inspect("btm", env = e)))
        } 
}

get.zero.constant <- function() -1e-4

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

pop.predict.half.child <- function(res.env, outdir, wenv) {
    # add res matrices
    present.and.proj.years <- wenv$present.and.proj.years
    present.and.proj.years.pop <- present.and.proj.years + 2
    npred <- wenv$npred
    npredplus1 <- npred + 1
    nvariants <- wenv$nvariants
    ages <- wenv$ages
    mx.ages <- wenv$mx.ages
    with(res.env, {
        totpm.hch <- totpf.hch <- array(0, dim=c(27, npredplus1, nvariants), 
                                        dimnames=list(ages, present.and.proj.years.pop, NULL))
        if(keep.vital.events) {
            btm.hch <- btf.hch <- asfert.hch <- pasfert.hch <- array(0, dim=c(7, npredplus1, nvariants),
                                                                     dimnames=list(NULL, present.and.proj.years, NULL))
            deathsm.hch <- deathsf.hch <- array(0, dim=c(27, npredplus1, nvariants), 
                                                dimnames=list(ages, present.and.proj.years, NULL))
            mxm.hch <- mxf.hch <- array(0, dim=c(28, npredplus1, nvariants), dimnames=list(mx.ages, present.and.proj.years, NULL))
        }
    })
    .ini.pop.res.env(res.env, wenv$keep.vital.events)
    for(cidx in 1:wenv$ncountries) {
        country <- wenv$country.codes[cidx]
        file.name <- file.path(outdir, paste0('totpop_country', country, '.rda'))
        computed.env <- new.env()
        load(file.name, envir=computed.env)
        medmigm <- apply(computed.env$migm, c(1,2), "median")
        medmigf <- apply(computed.env$migf, c(1,2), "median")
        medpopm <- apply(computed.env$totpm, c(1,2), "median")
        medpopf <- apply(computed.env$totpf, c(1,2), "median")
        res.env$totpm.hch[,1,] <- computed.env$totpm[,1,1]
        res.env$totpf.hch[,1,] <- computed.env$totpf[,1,1]
        if(wenv$keep.vital.events) {
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
            do.pop.predict.one.country.no.migration.half.child(time, cidx, res.env, 
                                        list(M=popM.hch.prev, F=popF.hch.prev), wenv)
            # add migration
            res.env$totpm.hch[,time+1,] <- res.env$totpm.hch[,time+1,] + medmigm[,time+1]
            res.env$totpf.hch[,time+1,] <- res.env$totpf.hch[,time+1,] + medmigf[,time+1]
            # for ages where both variants are equal, they should be equal to the pop median
            res.env$totpm.hch[,time+1,] <- adjust.half.child(res.env$totpm.hch[,time+1,], medpopm[,time+1])
            res.env$totpf.hch[,time+1,] <- adjust.half.child(res.env$totpf.hch[,time+1,], medpopf[,time+1])
            
            #if(keep.vital.events) {
            #    for(par in c('btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'pasfert.hch', 'mxm.hch','mxf.hch')) 
            #        res.env[[par]][,time+1,] <- nomigpred[[par]]	
            #}
            popM.hch.prev <- res.env$totpm.hch[,time+1,]
            popF.hch.prev <- res.env$totpf.hch[,time+1,]
        }
        spop <- res.env$totpm.hch + res.env$totpf.hch
        res.env$totp.hch <- apply(spop, c(2,3), sum, na.rm = TRUE)
        for(par in c('totp.hch', 'totpm.hch', 'totpf.hch'))
            computed.env[[par]] <- res.env[[par]]
        save(list=ls(computed.env, all.names=TRUE), envir=computed.env, file=file.name)
        if(wenv$keep.vital.events) {
            for(par in c('btm.hch', 'btf.hch', 'deathsm.hch', 'deathsf.hch', 'asfert.hch', 'pasfert.hch', 'mxm.hch','mxf.hch')) 
                computed.env.ve[[par]] <- res.env[[par]]
            save(list=ls(computed.env.ve, all.names=TRUE), envir=computed.env.ve, file=file.name.ve)
        }
    }
}

balanced.migration.1traj <- function(time, itraj, work.env, res.env) {
    lages <- length(res.env$ages21)
    for(cidx in 1:res.env$ncountries) {
        inpc <- res.env$countries.input[[res.env$country.codes.char[cidx]]]
        inpc$world.pop.distr <- res.env$world.pop.distr
        inpc$world.pop.distr.ini <- res.env$world.pop.distr.ini
        inpc$mig.rogers.castro <- res.env$mig.rogers.castro
        pop.ini <- if(time == 1) list(M = inpc$POPm0, F = inpc$POPf0) else
                          list(M = res.env$popM.prev[,cidx,drop=FALSE],
                               F = res.env$popF.prev[,cidx,drop=FALSE])
        do.pop.predict.1country.1traj.no.migration(time, itraj, cidx, res.env, pop.ini)
        pop <- sum(res.env$totpm[,cidx] + res.env$totpf[,cidx])
        migpred <- .get.migration.one.trajectory(res.env$use.migration.model, inpc, itraj, time, pop, 
                                                 popM=res.env$totpm[,cidx], popF=res.env$totpf[,cidx], country.code=res.env$country.codes[cidx], 
                                                 mig.rates=if(!is.null(res.env$mig.rate)) res.env$mig.rate[,cidx] else NULL, 
                                                 fixed.rate=if(res.env$fixed.mig.rate) inpc$projected.migration.rates[itraj,time] else NULL,
                                                 warn.template=res.env$warns[["_template_"]])
        #print(c(list(paste('Country:', country.codes.char[cidx], ', time: ', time, ', traj: ', itraj)), migpred))
        work.env$migrm[,cidx] <- migpred$M
        work.env$migrf[,cidx] <- migpred$F
        if(!is.null(migpred$laborM)) {
            work.env$migrm.labor[,cidx] <- migpred$laborM
            work.env$migrf.labor[,cidx] <- migpred$laborF
        }
        res.env$mig.rate[time+1, cidx] <- migpred$rate
        if (!is.null(migpred$warns)) {
            res.env$warns[[res.env$country.codes.char[cidx]]] <- if(is.null(res.env$warns[[res.env$country.codes.char[cidx]]])) migpred$warns else
                res.env$warns[[res.env$country.codes.char[cidx]]] + migpred$warns			
        }
    }
    work.env$popm <- res.env$totpm
    work.env$popf <- res.env$totpf
    # for cases when dimension is dropped (if there is one country)
    #if(is.null(dim(work.env$popm))) work.env$popm <- abind(work.env$popm, along=2)
    #if(is.null(dim(work.env$popf))) work.env$popf <- abind(work.env$popf, along=2)
    
    if(res.env$rebalance) {
        pop <- colSums(res.env$totpm + res.env$totpf)
        rebalance.migration2groups(work.env, pop, itraj)			
        negatives <- res.env$country.codes.char[unique(work.env$negatives)]
        for(country in negatives) add.pop.warn(country, time, 1, res.env) # 'Population negative while balancing'
    }
    res.env$migm[1:lages,] <- work.env$migrm + work.env$migrm.labor
    res.env$migf[1:lages,] <- work.env$migrf + work.env$migrf.labor
    
    # New population counts
    res.env$totpm <- res.env$totpm + res.env$migm
    res.env$totpf <- res.env$totpf + res.env$migf
    return(NULL)
}


do.pop.predict.1country.1traj.no.migration <- function(time, itraj, cidx, env, pop.ini) {
    country.name <- env$UNnames[cidx]
    inpc <- env$countries.input[[env$country.codes.char[cidx]]] 
    LeeC <- env$surv[[env$country.codes.char[cidx]]]
    pasfr <- env$kantor.pasfr[[env$country.codes.char[cidx]]][[itraj]][,time,drop=FALSE]
    asfr <- pasfr
    for(i in 1:nrow(asfr)) asfr[i,] <- inpc$TFRpred[time,itraj] * asfr[i,]
    LTres <- if(env$fixed.mx) sapply(LeeC, function(x) list(x[[1]][,time,drop=FALSE], 
                                                        x[[2]][,time,drop=FALSE]),
                                 simplify = FALSE, USE.NAMES = TRUE) else
                    project.mortality(inpc$e0Mpred[time,itraj], inpc$e0Fpred[time,itraj], npred = 1, 
                                      mortcast.args = env$mortcast.args[[env$country.codes.char[cidx]]], 
                                      verbose = env$verbose)

    popres <- PopProjNoMigr(1, pop.ini, LTres, asfr, inpc$SRB[time], country.name=country.name,
                            keep.vital.events=env$keep.vital.events)

    env$totpm[,cidx] <- popres$mpop[,2]
    env$totpf[,cidx] <- popres$fpop[,2]
    if(env$keep.vital.events) {
        env$btm[,cidx] <- popres$mbt
        env$btf[,cidx] <- popres$fbt
        env$deathsm[,cidx] <- popres$mdeaths
        env$deathsf[,cidx] <- popres$fdeaths
        env$asfert[,cidx] <- asfr
        env$pasfert[,cidx] <- pasfr*100
        env$mxm[,cidx] <- LTres$mx[[1]]
        env$mxf[,cidx] <- LTres$mx[[2]]
    }
    return(NULL)
}

do.pop.predict.one.country.no.migration.half.child <- function(time, cidx, env, pop.ini, wenv) {
    country.name <- wenv$UNnames[cidx]
    inpc <- wenv$countries.input[[wenv$country.codes.char[cidx]]] 
    LeeC <- wenv$surv[[wenv$country.codes.char[cidx]]]											    

	LTres <- if(wenv$fixed.mx) sapply(LeeC, function(x) list(x[[1]][,time,drop=FALSE], 
	                                                    x[[2]][,time,drop=FALSE]),
	                             simplify = FALSE, USE.NAMES = TRUE) else
	   project.mortality(inpc$e0Mmedian[time], inpc$e0Fmedian[time], npred = 1, 
	                     mortcast.args = wenv$mortcast.args[[wenv$country.codes.char[cidx]]], 
	                     verbose = wenv$verbose)

	for (variant in 1:wenv$nvariants) {
	    pasfr <- wenv$kantor.pasfr[[wenv$country.codes.char[cidx]]][[wenv$nr.traj+variant]][,time,drop=FALSE]
	    asfr <- pasfr
		for(i in 1:nrow(asfr)) asfr[i,] <- inpc$TFRhalfchild[variant,time] * asfr[i,]		
        this.pop.ini <- list(M = pop.ini$M[, variant], F = pop.ini$F[, variant])
		popres <- PopProjNoMigr(1, this.pop.ini, LTres, asfr, inpc$SRB[time], 
								country.name=country.name, keep.vital.events=wenv$keep.vital.events)
		env$totpm.hch[,time+1, variant] <- popres$mpop[,2]
		env$totpf.hch[,time+1, variant] <- popres$fpop[,2]
		if(wenv$keep.vital.events) {
		    env$btm.hch[,time+1,variant] <- popres$mbt
		    env$btf.hch[,time+1,variant] <- popres$fbt
		    env$deathsm.hch[,time+1,variant] <- popres$mdeaths
		    env$deathsf.hch[,time+1,variant] <- popres$fdeaths
		    env$asfert.hch[,time+1,variant] <- asfr
		    env$pasfert.hch[,time+1,variant] <- pasfr*100
		    env$mxm.hch[,time+1,variant] <- LTres$mx[[1]]
		    env$mxf.hch[,time+1,variant] <- LTres$mx[[2]]
		}
	}
	return(NULL)
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


project.migration.one.country.one.step <- function(mu, phi, sigma, oldRates, country.code, rlim=list(NULL, NULL), relaxed.bounds=FALSE, is.small=FALSE){
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
	if(!is.null(rlim[[2]])) xmax <- min(xmax, rlim[[2]])
	if(!is.null(rlim[[1]])) xmin <- max(xmin, rlim[[1]])
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
	if(country.code %in% land_area_wpp2019$country_code)
		land.area <- land_area_wpp2019[land_area_wpp2019$country_code==country.code,'land_area']
	i <- 0
	k <- 1
	zero.constant <- get.zero.constant()
	warns <- NULL
	popM21 <- c(popM[1:20], sum(popM[21:length(popM)]))
	popF21 <- c(popF[1:20], sum(popF[21:length(popF)]))
	popMdistr <- popM21/pop
	popFdistr <- popF21/pop
	emigrant.rate.bound <- -0.3
	country.code.char <- as.character(country.code)
	# adjustment constant for adjusting the rate by the population distribution
	popdistr <- (popM21 + popF21)/sum(popM21 + popF21)
	adj.constant.neg.numer <- sum(inpc$mig.rogers.castro * popdistr)
	pop0distr <- c(inpc$POPm0 + inpc$POPf0)/sum(c(inpc$POPm0 + inpc$POPf0))
	adj.constant.neg.denom <- sum(inpc$mig.rogers.castro * pop0distr)
	adj.constant.neg <- adj.constant.neg.numer/adj.constant.neg.denom
	adj.constant.pos.numer <- sum(inpc$mig.rogers.castro * inpc$world.pop.distr)
	adj.constant.pos.denom <- sum(inpc$mig.rogers.castro * inpc$world.pop.distr.ini)
	while(i <= 1000) {
		i <- i + 1
		if(is.null(fixed.rate)) {
			if(all(pars == 0)) rate <- 0
			else {
			    rlim1 <- if(pop>0 && !is.na(land.area)) -(pop - 0.0019*land.area)/pop else NULL
			    #rlim2a <- c(gcc.upper.threshold(country.code.char)/pop, if(!is.na(land.area)) 44*land.area/pop - 1 else NA)
			    rlim2a <- c(gcc.upper.threshold(country.code.char)/pop, if(!is.na(land.area)) exp(5.118 + 0.771*log(land.area))/pop - 1 else NA)
			    rlim <- list(rlim1, 
			                 if(pop>0 && any(!is.na(rlim2a))) min(rlim2a, na.rm=TRUE) else NULL)
			    rate <- project.migration.one.country.one.step(pars$mu, pars$phi, pars$sigma, 
					c(as.numeric(inpc$migration.rates), mig.rates[1:time]), country.code, rlim = rlim,
					relaxed.bounds=time < 6, is.small=pop < 200
					# max(colSums(inpc$observed$MIGm + inpc$observed$MIGf)
					)
			}
		} else rate <- fixed.rate
		if(is.na(rate)) stop('Migration rate is NA for country ', country.code, ', time ', time, ', traj ', itraj, 
					'.\npop=', paste(pop, collapse=', '), '\nmig rate=', paste(c(as.numeric(inpc$migration.rates), mig.rates[1:time]), collapse=', '))
		# adjustment of the rate by the population distribution
		if(is.null(fixed.rate)) {
		    if(rate < 0) { # adjust using the country pop distr
		        adj.constant <- adj.constant.neg.numer/adj.constant.neg.denom
		    } else { # adjust using the world pop distr
		        adj.constant <- adj.constant.pos.numer/adj.constant.pos.denom
		    }
		    #browser()
		    #print(c(country.code, time, itraj, ":", round(adj.constant,3), round(rate, 3), round(rate * adj.constant,3), 
		    #        round(rate * pop), round(rate * adj.constant * pop)))
		    rate <- rate * adj.constant
		}
		mig.count <- rate * pop
		if(rate < 0 && mig.count < emigrant.rate.bound*pop && i < 1000 && is.null(fixed.rate)) next # resample if the outmigration would be larger than what is allowed
		schedMname <- 'M'
		schedFname <- 'F'
		if(rate < 0 && !is.null(inpc$migration.age.schedule[['Mnegative']])) {
			schedMname <- 'Mnegative'
			schedFname <- 'Fnegative'
		}
		msched <- inpc$migration.age.schedule[[schedMname]][,time]
		fsched <- inpc$migration.age.schedule[[schedFname]][,time]
		 if(FALSE && is.gcc(country.code)) { # This is switched off for 2300 projections
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
		#if(rate < 0 && !is.gcc(country.code)) {
		if(rate < 0) { # this is on for GCC countries for the 2300 projections 
				denom <- sum(msched * popMdistr + fsched * popFdistr)
				if(denom == 0) { # not enough people to migrate out
				    msched[] <- 0
				    fsched[] <- 0
				} else {
				    denom2 <- c(msched, fsched)/denom
				    if(abs(rate) > min((abs(emigrant.rate.bound) / denom2)[denom2 > 0]) && i < 1000 && is.null(fixed.rate)) next
				    msched <- msched * popMdistr / denom
				    fsched <- fsched * popFdistr / denom
				}
		}
		# age-specific migration counts
		migM <- mig.count*msched
		migF <- mig.count*fsched
		if(!is.null(fixed.rate) || rate == 0 || mig.count == 0) break
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
    nr.traj <- dim(env$migm)[3]
    e <- new.env()
    migbalattr <- list(m = "migrm", f = "migrf")
    migresattr <- list(m = "migm", f = "migf")
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
    cleanup.env(e)
}

labor.countries <- function()
	return(prop.labor.migration()$country_code)
	
prop.labor.migration <- function()
	return(migration.thresholds$labor.proportions)
		
prop.labor.migration.for.country <- function(code){
	props <- prop.labor.migration()
	return(props[props$country_code==code,'proportion'])
}


restructure.pop.data.and.compute.quantiles <- function(source.dir, dest.dir, npred, nr.traj, inputs, observed, 
									present.and.proj.years, keep.vital.events=FALSE, parallel=FALSE, nr.nodes=NULL, 
									verbose=FALSE, ...){
	
	restructure.pop.data.and.compute.quantiles.one.country <- function(cidx) {		
		country <- country.codes[cidx]
		inpc <- inputs[[country]]
		obs <- observed[[country]]
		MxKan <- runKannisto.noLC(inputs[[country]], observed = TRUE)
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
				btm <- btf <- asfert <- pasfert <- array(NA, dim=c(7, npredplus1, nr.traj), 
							dimnames=list(NULL, present.and.proj.years, NULL))
				deathsm <- deathsf <- array(NA, dim=c(27, npredplus1, nr.traj), 
							dimnames=list(ages, present.and.proj.years, NULL))
				mxm <- mxf <- array(NA, dim=c(28, npredplus1, nr.traj), dimnames=list(mx.ages, present.and.proj.years, NULL))

				# values from current year		
				btm[1:dim(obs$btm)[1],1,] <- obs$btm[,dim(obs$btm)[2],repi]
				btf[1:dim(obs$btf)[1],1,] <- obs$btf[,dim(obs$btf)[2],repi]
				deathsm[1:dim(obs$deathsm)[1],1,] <- obs$deathsm[,dim(obs$deathsm)[2],repi]
				deathsf[1:dim(obs$deathsf)[1],1,] <- obs$deathsf[,dim(obs$deathsf)[2],repi]
				asfert[1:dim(obs$asfert)[1],1,] <- obs$asfert[,dim(obs$asfert)[2],repi]
				pasfert[1:dim(obs$pasfert)[1],1,] <- obs$pasfert[,dim(obs$pasfert)[2],repi]
				mxm[1:dim(MxKan[[1]]$mx)[1],1,] <- as.matrix(MxKan[[1]]$mx[,dim(MxKan[[1]]$mx)[2], drop = FALSE])[,repi]
				mxf[1:dim(MxKan[[2]]$mx)[1],1,] <- as.matrix(MxKan[[2]]$mx[,dim(MxKan[[2]]$mx)[2], drop = FALSE])[,repi]
			}
		})
		for(time in 1:npred) {
			for(par in c('totp'))
				res.env[[par]][time+1,] <- envs[[time]][[par]][cidx,]
			for(par in c('totpm', 'totpf', 'migm', 'migf'))
				res.env[[par]][,time+1,] <- envs[[time]][[par]][,cidx,]
		}
		if(keep.vital.events) {
		    for(itraj in 1:nr.traj) {
				for(par in c('btm', 'btf', 'deathsm', 'deathsf', 'asfert', 'pasfert', 'mxm', 'mxf'))
				    res.env[[par]][,2:npredplus1,itraj] <- t(envs.ve[[itraj]][[par]][,,cidx])
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
		cleanup.env(res.env)
		return(quant.env)
	}

	envs <- list()
	for(time in 1:npred) {
	    envs[[time]] <- new.env()
	    load(file.path(source.dir, paste0('pop_time_', time, '.rda')), envir=envs[[time]])
	}
    if(keep.vital.events) {
        envs.ve <- list()
	    for(itraj in 1:nr.traj) {
            envs.ve[[itraj]] <- new.env()
	        load(file.path(source.dir, paste0('vital_events_traj_', itraj, '.rda')), envir=envs.ve[[itraj]])
	        #load(file.path(source.dir, paste0('vital_events_time_', time, '.rda')), envir=envs[[time]])
	   }
	}
	ncountries <- nrow(envs[[1]]$totp)
	country.codes <- rownames(envs[[1]]$totp)
	ages <- dimnames(envs[[1]]$totpm)[[1]]
	nages <- length(ages)
	mx.ages <- c(0,1,ages[2:nages])
	npredplus1 <- npred + 1
	present.and.proj.years.pop <- present.and.proj.years + 2
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
	# cleanup
	rm(res.list)
	for(i in 1:length(envs)) cleanup.env(envs[[i]])
	rm(envs)
	if(keep.vital.events) {
	    for(i in 1:length(envs.ve)) cleanup.env(envs.ve[[i]])
	    rm(envs.ve)
	}
	return(quant.env)
}

cleanup.env <- function(env) {
    rm(list = ls(env), envir = env)
    rm(env)
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
	if((schc <- .pattern.value("MigAgeSchedule", mig.settings, 0)) > 0) {
		   sched.country <- schc
		   first.year <- TRUE
	}
	# Should the Male/Female ratio be kept or set equal. E.g. China schedule has larger migration for female, so rescale
	if(.pattern.value("MigAgeEqualMFratio", mig.settings, 0) == 1) 
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
		if(!is.null(scale.to.totals) && !any(colSums(maleV) == 0) && !any(colSums(femaleV) == 0) 
		    ){
			maleV <- t(scale.to.totals$M * apply(maleV, 1, '/', colSums(maleV)))
			femaleV <- t(scale.to.totals$F * apply(femaleV, 1, '/', colSums(femaleV)))
		}
	    colnames(maleV) <- colnames(femaleV) <- colnames(inputs$MIGm)[col.idx]
	    tot <- colSums(maleV+femaleV)
	    if(any(abs(tot) <= 0.0001)) {
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
		non.zero <- abs(tot) > 0.0001
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
    #if(country == 64) stop("")
    if(is.gcc(country)) {
    	negF <- femaleArray # female gets China schedule; should be correctly scaled
    	# male - use model out-migration schedule (derived from SA)
    	negMvec <- gcc.model.outschedule()
    	negMvec <- negMvec/sum(negMvec)
    	negM <- matrix(negMvec, nrow=nAgeGroups, ncol=npred)*matrix(scale[1:npred], ncol=npred, nrow=nAgeGroups, byrow=TRUE) # scale
    }
    # For some countries like Egypt, if positive migration rate, set negative schedules to zero, since they would mean out-migration
    if(.pattern.value("MigAgeZeroNeg", mig.settings, 0) == 1) {
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