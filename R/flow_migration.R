if(getRversion() >= "2.15.1") utils::globalVariables(c("land_area_wpp2019", "migration.thresholds"))

do.pop.predict.flow <- function(inp, outdir, nr.traj, ages, batch=1, batch.size=10, pred=NULL, countries=NULL, keep.vital.events=FALSE, 
                                fixed.mx=FALSE, fixed.pasfr=FALSE, function.inputs=NULL, start.time.index=1,
								verbose=TRUE, compute.summary=FALSE, migration.settings=NULL, ...) {
	# if migration.settings is NULL, don't use the migration model
	if (is.null(migration.settings))
		stop('Argument migration.settings has to be given if use.migration.flow.model is TRUE.')
	countries.idx <- if(is.null(countries)) which(UNlocations$location_type==4) else which(UNlocations$country_code %in% countries)
	country.codes <- UNlocations$country_code[countries.idx]
	country.codes.char <- as.character(country.codes)
	ncountries <- length(country.codes)
	nr_project <- length(inp$proj.years)
	nages <- length(ages)
	mx.ages <- c(0,1,ages[2:nages])
	ages21 <- ages[ages<=100]
	nest <- length(inp$estim.years)
	outdir.batch <- file.path(outdir, paste0("batch_", batch))
	if(!file.exists(outdir.batch)) dir.create(outdir.batch, recursive=TRUE)
	present.and.proj.years <- c(inp$estim.years[nest], inp$proj.years)
	present.and.proj.years.pop <- present.and.proj.years + 2
	inp.to.save <- list()
	# remove big or redundant items from inputs to be saved
	for(item in ls(inp)[!grepl('^migMpred$|^migFpred$|^TFRpred$|^e0Fpred$|^e0Mpred$|^estim.years$|^proj.years$|^wpp.years$', ls(inp))]) 
		inp.to.save[[item]] <- get(item, inp)
	npred <- nr_project

	# set intitial outflow rate values
	mig.out.rates.ini <- migration.settings$ini.out.rates[country.codes.char]
	flow.obs <- migration.settings$flowData

	# save a Roger's Castro curve (taken from China)
	migRC <- age.specific.migration(wpp.year=inp$wpp.year, countries=156, verbose=FALSE)
	mig.rogers.castro <- migRC$male[,"2020-2025"]/sum(migRC$male[,"2020-2025"])

	outdir.raw <- file.path(outdir.batch, 'raw')
	if(file.exists(outdir.raw) && start.time.index==1) unlink(outdir.raw, recursive=TRUE)
	if(start.time.index==1) dir.create(outdir.raw)
    
    fixed.mig.rate <- FALSE    
	UNnames <- UNlocations[countries.idx,'name']
	countries.input <- new.env()
	# Extract the country-specific stuff from the inputs
	if(verbose) cat('\nLoading inputs for ', ncountries, ' countries (sequentially)')
	for(cidx in 1:ncountries) {
		inpc <- get.country.inputs(country.codes[cidx], inp, nr.traj, UNnames[cidx])
		if(is.null(inpc) || length(inpc$POPm0)==0) next
		countries.input[[country.codes.char[cidx]]] <- as.environment(inpc)
	}
	nr.traj <- min(nr.traj, sapply(ls(countries.input), function(ccode) return(ncol(countries.input[[ccode]]$TFRpred))))
	npred <- min(npred, sapply(ls(countries.input), function(ccode) return(nrow(countries.input[[ccode]]$TFRpred))))
	npredplus1 <- npred+1
	npasfr <- nrow(countries.input[[country.codes.char[1]]]$PASFR)
	nvariants <- nrow(countries.input[[country.codes.char[1]]]$TFRhalfchild)
	if(verbose) cat('\nProcessing ', npred, ' time periods for ', ncountries, ' countries for ', nr.traj, ' trajectories (sequentially)')

	# generate batch subset indicies
	b <- min(batch, 10)
	trajSubsetIndex <- (batch.size*(b-1) + 1):(batch.size*b)
	nr.traj.subset <- length( trajSubsetIndex )
	muColNames <- c("country_code", as.character(trajSubsetIndex) )

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
		tfr.med <- apply(countries.input[[country]]$TFRpred[,trajSubsetIndex], 1, median)[nrow(countries.input[[country]]$TFRpred)]
		kantor.pasfr[[country]] <- list()
		for(itraj in 1:nr.traj.subset)	
			kantor.pasfr[[country]][[itraj]] <- if(fixed.pasfr) countries.input[[country]]$PASFR/100. else
			                                           kantorova.pasfr(c(countries.input[[country]]$observed$TFRpred, countries.input[[country]]$TFRpred[,trajSubsetIndex[itraj]]), 
										countries.input[[country]], norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)
		
		for (variant in 1:nvariants) 
			kantor.pasfr[[country]][[variant + nr.traj.subset]] <- if(fixed.pasfr) countries.input[[country]]$PASFR/100. else                      
			                                                       kantorova.pasfr(c(countries.input[[country]]$observed$TFRpred, 
			                            countries.input[[country]]$TFRhalfchild[variant,]), 
										countries.input[[country]], norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)						
	}
	debug <- FALSE

	# set flow model parameter values
	mu.oos <- migration.settings$mu.oos[,muColNames,with=FALSE]
	phi.oos <- migration.settings$phi.oos[trajSubsetIndex]
	sigma.oos <- migration.settings$sigma.oos[trajSubsetIndex]

	mu <- migration.settings$mu[,muColNames,with=FALSE]
	phi <- migration.settings$phi[trajSubsetIndex]
	sigma <- migration.settings$sigma[trajSubsetIndex]

	pi.dir <- migration.settings$pi.dir

	country.codes.char <- country.codes[order(as.numeric(country.codes))] %>% as.character()
    
    # prepare environment for storing results
    res.env <- new.env()
    with(res.env, {
    	totp <- matrix(0, nrow=ncountries, ncol=nr.traj.subset, dimnames=list(country.codes.char, NULL))
    	totpm <- totpf <- array(0, dim=c(27, ncountries, nr.traj.subset), dimnames=list(ages, country.codes.char, NULL))
    	flow.age.female <- flow.age.male <- array(NA, dim=c(ncountries, length(ages21), ncountries, nr.traj.subset), 
        							dimnames=list("origin"=country.codes.char, "age"=ages21, "destination"=country.codes.char, "trajectory"=1:nr.traj.subset))
    	out.mig.rate <- array(NA, dim=c(ncountries, npred+1, nr.traj.subset), dimnames=list(country.codes.char, 1:(npred+1), 1:nr.traj.subset))
    	warns <- list()
    	}
    )
    
    # prepare working environment for storing population for one trajectory
    create.work.env <- function() {
    	work.env <- new.env()
    	mapply(assign, global.objects, mget(global.objects, inherits = TRUE), MoreArgs = list(envir = work.env))
		with(work.env, {
    	    totpm <- totpf <- array(0, dim=c(27, ncountries), dimnames=list(ages, country.codes.char))
    	    flow.age.female <- flow.age.male <- array(NA, dim=c(ncountries, length(ages21), ncountries), 
        							dimnames=list("origin"=country.codes.char, "age"=ages21, "destination"=country.codes.char))
    	    out.mig.rate <- matrix(NA, nrow=ncountries, ncol=npred+1, dimnames=list(country.codes.char, 1:(npred+1)))
    	    names(out.mig.rate) <- country.codes.char
    	    inmigm <- outmigm <- netmigm <- inmigf <- outmigf <- netmigf <- matrix(NA, nrow=27, ncol=ncountries, dimnames=list(ages, country.codes.char))
			if(keep.vital.events) {
			    btm <- btf <- asfert <- pasfert <- array(0, dim=c(7, ncountries), dimnames=list(NULL, country.codes.char))
    	        deathsm <- deathsf <- array(0, dim=c(27, ncountries), dimnames=list(ages, country.codes.char))
               	mxm <- mxf <- array(0, dim=c(28, ncountries), dimnames=list(mx.ages, country.codes.char))
            }
			warns <- list()
			warns[["_template_"]] <- matrix(0, nrow=get.nr.warns(), ncol=npred)
			}
		)
        return(work.env)
    }
    create.mig.env <- function() {
        work.env <- new.env()
        with(work.env, {
        	flow.age.female <- flow.age.male <- array(NA, dim=c(ncountries, length(ages21), ncountries), 
        							dimnames=list("origin"=country.codes.char, "age"=ages21, "destination"=country.codes.char))
        })
        return(work.env)
    }
	# To have a consistent sequential application, we make these global
	global.objects <- c("nr.traj", "nr.traj.subset", "trajSubsetIndex", 
						"country.codes",  "UNnames", "countries.input", 
	                    "mortcast.args", "surv", "npasfr", "fixed.mx",
	                    "ages", "ages21", "mx.ages", "ncountries",
	                    "nvariants", "keep.vital.events", "verbose",
	                    "npred", "country.codes.char", "kantor.pasfr", "outdir.raw", "outdir.batch",
	                    "mu.oos", "phi.oos", "sigma.oos", "mu", "phi", "sigma", "pi.dir", 
	                    "mig.rogers.castro")
	work.env <- create.work.env() 
	work.mig.env <- create.mig.env()
	
	data(land_area_wpp2019)	

	
	predict.1time.period.1flow.trajectory <- function(itraj, time) {
	    with(work.mig.env, {
	        flow.age.male[] <- flow.age.female[] <- 0
	    })
	    .ini.pop.flow.res.env(work.env, keep.vital.events)
	    if(time==1) work.env$out.mig.rate.prev <- mig.out.rates.ini
	    if(time > 1) {
	        work.env$popM.prev <- popM.prev[,,itraj]
	        work.env$popF.prev <- popF.prev[,,itraj]
	        work.env$out.mig.rate.prev <- out.mig.rate.prev[,itraj] 
	    }
	    if(time==2) work.env$out.mig.cum <- out.mig.cum[,itraj]
	    if(time>2) work.env$out.mig.cum <- out.mig.cum[,,itraj]
	    get.flow.1traj(time, itraj, work.env = work.mig.env, res.env = work.env)
	    res <- list()
	    for(item in c('totpm', 'totpf', 'out.mig.rate', 'flow.age.male', 'flow.age.female', 'warns'))
	        res[[item]] <- work.env[[item]]
	    return(res)
	}
	
	for(time in start.time.index:npred) {
		if(verbose) cat("\nProcessing time period", time, "at", format(Sys.time(), "%X") )
			.ini.pop.res.env(res.env, vital.events = FALSE)
		thispred <- list()
		for(itraj in 1:nr.traj.subset) {
		    thispred[[itraj]] <- predict.1time.period.1flow.trajectory(itraj, time)
		}
		# collect results
		for(itraj in 1:nr.traj.subset) { 

			for( par in c('totpm', 'totpf', 'out.mig.rate') )
				res.env[[par]][,,itraj] <- thispred[[itraj]][[par]]
			for( par in c('flow.age.male', 'flow.age.female') )
				res.env[[par]][,,,itraj] <- thispred[[itraj]][[par]]
			res.env$out.mig.rate[,,itraj] <- thispred[[itraj]]$out.mig.rate
			for( country in setdiff(names(thispred[[itraj]]$warns), "_template_") ) {
			    if(is.null(res.env$warns[[country]]))
			        res.env$warns[[country]] <- thispred[[itraj]]$warns[["_template_"]]
				res.env$warns[[country]] <- res.env$warns[[country]] + thispred[[itraj]]$warns[[country]]
			}
		}
		rm(thispred)
		spop <- res.env$totpm + res.env$totpf
		res.env$totp <- if(dim(res.env$totp)[1]==1) sum(spop) else apply(spop, c(2,3), sum) # distinction if there is only one country
		popM.prev <- res.env$totpm
		popF.prev <- res.env$totpf
		out.mig.cum <- res.env$out.mig.rate[,1:time,] 
		out.mig.rate.prev <- res.env$out.mig.rate[,time,]
		if (any(res.env$totpm < get.zero.constant()) || any(res.env$totpf < get.zero.constant())){
			cntries.m <- which(apply(res.env$totpm, 2, function(x) any(x < get.zero.constant())))
			cntries.f <- which(apply(res.env$totpf, 2, function(x) any(x < get.zero.constant())))
			for(country in unique(c(cntries.m, cntries.f))) {
				neg.times <- unique(which(apply(res.env$totpm[,country,], 2, function(x) any(x<0))),
				which(apply(res.env$totpf[,country,], 2, function(x) any(x<0))))
				add.pop.warn(country.codes.char[country], neg.times, 5, res.env)  #'Final population negative for some age groups'
			}
		}
		with(res.env, {
			file.name <- file.path(outdir.raw, paste0('pop_time_', time, '.rda'))
			colnames(totp) <- as.character( 1:nr.traj.subset )
			save(totp, totpm, totpf, file = file.name)

			female.flow.file.name <- file.path(outdir.raw, paste0('flow_female_time_', time, '.rda'))
			save(flow.age.female, file = female.flow.file.name)

			male.flow.file.name <- file.path(outdir.raw, paste0('flow_male_time_', time, '.rda'))
			save(flow.age.male, file = male.flow.file.name)
		})
	} # end time
	cleanup.env(work.mig.env)

	# update restructure.pop.data.and.compute.quantiles and bayesPop.prediction quantity as needed
	if(compute.summary){
		if(verbose) cat('\nRe-formatting data (sequentially)\n')
		quant.env <- restructure.flow.data.and.compute.quantiles(source.dir=outdir, nr.traj=nr.traj, 
			proj.years=inp$proj.years, country.codes.char=country.codes.char, verbose=verbose, ...)
		if(verbose) cat(' done.\n')
		cleanup.env(res.env)
		cleanup.env(work.env)
		#save meta file
		bayesPop.flow.prediction <- structure(list(
								nr.traj = nr.traj,
								popPI = quant.env$popPIarray,
								flowMarginPI = quant.env$flowMarginPIarray, 
								flowPI = quant.env$flowPIarray,
								globeMigratingPI = quant.env$globeMigPIarray,
    	           				proj.years = present.and.proj.years[-1], # excludes present period (middle of periods)
    	           				wpp.year = inp$wpp.year,
				   				countries = country.codes.char,
				   				warnings = res.env$warns), class='bayesPop.flow.prediction')
		prediction.file <- file.path(outdir, 'flow_forecast_summary.rda')
		save(bayesPop.flow.prediction, file=prediction.file)
		cat('\nPrediction stored into', outdir, '\n')
		print.pop.warnings(bayesPop.flow.prediction, which.warns=c(2,3,5))
		# cleanup 
		cleanup.env(countries.input)
		cleanup.env(kantor.pasfr)
		cleanup.env(observed)
		cleanup.env(quant.env)
	
		return(bayesPop.flow.prediction)
		} else {
		cat('\nProcessing complete.\n')
		cleanup.env(countries.input)
		cleanup.env(kantor.pasfr)
		cleanup.env(observed)
		return(NULL)
	}

}



get.flow.1traj <- function(time, itraj, work.env, res.env) {
    lages <- length(res.env$ages21)
    for(cidx in 1:res.env$ncountries) {
    	ccc <- res.env$country.codes.char[cidx]
        inpc <- res.env$countries.input[[ccc]]
        inpc$country.codes.char <- res.env$country.codes.char
        inpc$origin.code.char <- ccc
        inpc$pi.dir <- res.env$pi.dir
        inpc$mig.rogers.castro <- res.env$mig.rogers.castro
        inpc$out.mig.rate.prev <- res.env$out.mig.rate.prev[ccc]
        if(time==2) inpc$out.mig.cum <- res.env$out.mig.cum[ccc]
        if(time>2) inpc$out.mig.cum <- res.env$out.mig.cum[ccc,]
        inpc$ages21 <- res.env$ages21
        if(is.incomplete.flow(ccc)){
        	inpc$mu <- as.numeric( res.env$mu.oos[country_code==ccc,-"country_code"] )[itraj]
        	inpc$sigma <- as.numeric( res.env$sigma.oos )[itraj]
        	inpc$phi <- as.numeric( res.env$phi.oos )[itraj]
        } else {
        	inpc$mu <- as.numeric( res.env$mu[country_code==ccc,-"country_code"] )[itraj]
        	inpc$sigma <- as.numeric( res.env$sigma )[itraj]
        	inpc$phi <- as.numeric( res.env$phi )[itraj]
        }
        pop.ini <- if(time == 1) list(M = inpc$POPm0, F = inpc$POPf0) else
                          list(M = res.env$popM.prev[,cidx,drop=FALSE],
                               F = res.env$popF.prev[,cidx,drop=FALSE])
        # run itraj index forwrad for all parameters other than migration so that we can execute a poor man parallel job
        do.pop.predict.1country.1traj.no.migration(time, itraj, cidx, res.env, pop.ini)
        pop <- sum(res.env$totpm[,cidx] + res.env$totpf[,cidx])
        outmigpred <- sample.migration.trajectory(inpc, itraj, time, pop, 
    									popM=res.env$totpm[,cidx], popF=res.env$totpf[,cidx], 
                                    	country.code=res.env$country.codes[cidx],
                                    	warn.template=res.env$warns[["_template_"]])
        res.env$flow.age.male[ccc,,] <- outmigpred$flow.age.m
        res.env$flow.age.female[ccc,,] <- outmigpred$flow.age.f
		#if(time==1) res.env$out.mig.rate[ccc, 1] <- res.env$out.mig.rate.prev[ccc]
        res.env$out.mig.rate[ccc, time] <- outmigpred$out.mig.rate
    }

    work.env$popm <- res.env$totpm
    work.env$popf <- res.env$totpf

    # net flows by age
    for(cidx in 1:res.env$ncountries) {
    	ccc <- res.env$country.codes.char[cidx]
    	
    	res.env$outmigm[1:lages,ccc] <- rowSums( res.env$flow.age.male[ccc,,] )
    	res.env$inmigm[1:lages,ccc]  <- colSums( res.env$flow.age.male[,,ccc] )
    	res.env$netmigm[1:lages,ccc] <- res.env$inmigm[1:lages,ccc] - res.env$outmigm[1:lages,ccc]

    	res.env$outmigf[1:lages,ccc] <- rowSums( res.env$flow.age.female[ccc,,] )
    	res.env$inmigf[1:lages,ccc]  <- colSums( res.env$flow.age.female[,,ccc] )
    	res.env$netmigf[1:lages,ccc] <- res.env$inmigf[1:lages,ccc] - res.env$outmigf[1:lages,ccc]
    }
    
    # New population counts
    res.env$totpm <- res.env$totpm + ( res.env$netmigm/1e3 )
    res.env$totpf <- res.env$totpf + ( res.env$netmigf/1e3 )

    return(NULL)
}





sample.migration.trajectory <- function(inpc, itraj=NULL, time=NULL, pop=NULL, 
									popM=NULL, popF=NULL, country.code=NULL, warn.template=NULL) {
	# TO DO: 
    # - Switch Rogers-Castro with country-specific age schedules or model schedules

	zero.constant <- get.zero.constant()
	warns <- NULL
	popM21 <- popM[1:21]
	popF21 <- popF[1:21]
	popMdistr <- popM21/pop
	popFdistr <- popF21/pop
	age21 <- inpc$ages21
	mig.age.schedule.f <- inpc$mig.rogers.castro
	mig.age.schedule.m <- inpc$mig.rogers.castro

	origin <- inpc$origin.code.char
	piDestTrajectory <- get.pi.1trajectory(country.code=origin, itraj=itraj, dir=inpc$pi.dir)
	mu <- inpc$mu
	phi <- inpc$phi
	sigma <- inpc$sigma
	out.mig.rate.prev <- inpc$out.mig.rate.prev

	rho <- mu * (1 - phi) + phi * out.mig.rate.prev

    maxIter <- iter <- 0
    while(iter<100){
		out.mig.rate.prop <- rgamma(n=1, shape=rho*sigma, rate=sigma)

		if(time==1){
			popdrop <- 1 - out.mig.rate.prop/200
    	} else {
    		out.mig.cum <- inpc$out.mig.cum 
    	    popdrop <- tail( cumprod(1 - c( out.mig.cum, out.mig.rate.prop)/200 ), 1 )
    	}


    	if( !is.special.case(origin) & !is.gcc(origin) & !is.labor(origin) ){
        	popdrop.lower.bound <- general.outflow.bounds(time)
        	if( popdrop >= popdrop.lower.bound ) break
        } 
        
        if( is.special.case(origin) ){
        	popdrop.lower.bound <- special.case.outflow.bounds(time)
        	if( popdrop >= popdrop.lower.bound ) break
        } 
        
        if( is.gcc(origin) & !is.special.case(origin) ){
        	popdrop.lower.bound <- gcc.outflow.bounds(time)
        	if( popdrop >= popdrop.lower.bound ) break
        } 
        
        if( is.labor(origin) & !is.special.case(origin) ){
        	popdrop.lower.bound <- labor.outflow.bounds(time)
        	if( popdrop >= popdrop.lower.bound ) break
        } 
        
        iter <- iter+1
    }

    out.mig.rate.raw <- out.mig.rate.prop 
    out.mig.count.raw <- round( out.mig.rate.raw * 5 * pop ) 
    # rate = migrants / (5 * Pop / 1000) = migrants * (1 / (5*Pop/1000)) = migrants * 1 / 5*pop000
    # 	=> migrants = rate * 5 * pop000 with pop000 = Pop/1000

    # age, sex, and migration schedule 
    outflow.sched.m <- popMdistr * ( mig.age.schedule.m )
    outflow.sched.f <- popFdistr * ( mig.age.schedule.f )

    outflow.sched <- c( outflow.sched.f, outflow.sched.m ) / sum( outflow.sched.f + outflow.sched.m )

    outmigF <- out.mig.count.raw * outflow.sched[1:21]
    outmigM <- out.mig.count.raw * outflow.sched[22:42]

    # check that no age-sex group drops by more than 40%
    outflow.sched.supremum <- 0.4 * 1e3 * c(popF21*((mig.age.schedule.f>0)*1), popM21*((mig.age.schedule.m>0)*1))
    if( any( c(outmigF, outmigM) > outflow.sched.supremum) ){ # adjust age schedules
		prev.isneg <- rep(FALSE, 42)
		j <- 1
		while( any( c( outmigF, outmigM ) > outflow.sched.supremum) ) { 
			isneg <- prev.isneg | c(outmigF, outmigM) > outflow.sched.supremum
			shifts <- - c(outmigF, outmigM) + outflow.sched.supremum
			shifts[!isneg] <- 0
			shifts[prev.isneg] <- 0
			sched.new <- outflow.sched + shifts/(out.mig.count.raw)
			delta <- sum(outflow.sched - sched.new)
			sched.new[!isneg] <- sched.new[!isneg] + delta*abs(sched.new[!isneg])/sum(abs(sched.new[!isneg]))
			sched.f <- sched.new[1:21]
			sched.m <- sched.new[22:42]	
			outmigF <- out.mig.count.raw * sched.f
			outmigM <- out.mig.count.raw * sched.m
			prev.isneg <- isneg
			outflow.sched <- sched.new
			j <- j+1
			if(j>100) break
			}

		# if the simulation broke before an acceptable shift, then set the remaining outflows to the maximum
		if( any( outmigF > outflow.sched.supremum[1:21] ) ){
			getsMaxOutflow <- which( outmigF > outflow.sched.supremum[1:21] )
			outmigF[ getsMaxOutflow ] <- outflow.sched.supremum[1:21][ getsMaxOutflow ]
		}

		if( any( outmigM > outflow.sched.supremum[22:42] ) ){
			getsMaxOutflow <- which( outmigM > outflow.sched.supremum[22:42] )
			outmigM[ getsMaxOutflow ] <- outflow.sched.supremum[22:42][ getsMaxOutflow ]
		}

		out.mig.count.raw <- sum(outmigF + outmigM)
		outflow.sched <- c(outmigF, outmigM) / out.mig.count.raw

    }


    # once the exact outflow is set, then determine where migrants will go
    piDest.sample <- piDestTrajectory$pi.traj
    piDest.codes <- piDestTrajectory$dest.codes
    flow.vector.raw <- rmultinom( n=1, size=out.mig.count.raw, prob=piDest.sample ) %>% as.numeric()
    flow.matrix.raw <- matrix(flow.vector.raw, nrow=21, ncol=199, byrow=TRUE, dimnames=list(age21, piDest.codes))

   	flow.age.schedule.f <- matrix( outflow.sched[1:21], nrow=21, ncol=199, dimnames=list(age21, piDest.codes))
    flow.age.schedule.m <- matrix( outflow.sched[22:42], nrow=21, ncol=199, dimnames=list(age21, piDest.codes))

    flow.age.f.raw <- round( flow.matrix.raw * flow.age.schedule.f ) %>% cbind(0)
    flow.age.m.raw <- round( flow.matrix.raw * flow.age.schedule.m ) %>% cbind(0)

    # add origin column and name back in to simplify post-processing
    colnames(flow.age.f.raw)[200] <- colnames(flow.age.m.raw)[200] <- origin
    flow.age.f <- flow.age.f.raw[,inpc$country.codes.char]
    flow.age.m <- flow.age.m.raw[,inpc$country.codes.char]

    out.mig.count <- sum(flow.age.f + flow.age.m)
    out.mig.rate <- out.mig.count / ( 5 * pop )

    return(list(out.mig.rate=out.mig.rate, out.mig.count=out.mig.count, flow.age.f=flow.age.f, flow.age.m=flow.age.m) ) 
}


restructure.flow.data.and.compute.quantiles <- function(source.dir, nr.traj, proj.years, country.codes.char, verbose=FALSE, ...){

	probs <- c(0.05, 0.1, 0.5, 0.9, 0.95)
	nprobs <- length(probs)
	ncountries <- length(country.codes.char)
	npred <- length( proj.years )

	quant.env <- new.env()
	with(quant.env, {
		flowPIarray <- 
		array(NA, 
			dim=c(nprobs, ncountries, ncountries, npred), 
			dimnames=list(quantile=probs, orig=country.codes.char, dest=country.codes.char, year=proj.years)
			)
		flowMarginPIarray <- 
			array(NA, 
				dim=c(nprobs, ncountries, 3, npred),
				dimnames=list(quantile=probs, country_code=country.codes.char, margin=c("outflow", "inflow", "netflow"), year=proj.years))
		popPIarray <- 
			array(NA, 
				dim=c(nprobs, ncountries, npred), 
				dimnames=list(quantile=probs, country_code=country.codes.char, year=proj.years))
		globeMigPIarray <- array(NA, dim=c(nprobs, npred), dimnames=list(quantile=probs, year=proj.years))
	})

	work.env <- new.env()

	for(t in 1:npred){

		cat("Processing time period", t, "at", format(Sys.time(), "%X"), "\n")
  		load(paste0(source.dir, "/raw/flow_female_time_", t, ".rda"), envir=work.env)
  		load(paste0(source.dir, "/raw/flow_male_time_", t, ".rda"), envir=work.env)
  		load(paste0(source.dir, "/raw/pop_time_", t, ".rda"), envir=work.env)
  
 		work.env$od.female = 
    		apply(work.env$flow.age.female, MARGIN=c("origin", "destination", "trajectory"), sum)
  		work.env$od.male = 
  		  apply(work.env$flow.age.male, MARGIN=c("origin", "destination", "trajectory"), sum)
  		work.env$od = work.env$od.female + work.env$od.male
  	
  		work.env$outflows = apply( work.env$od, MARGIN=c("origin", "trajectory"), sum)
  		work.env$inflows = apply( work.env$od, MARGIN=c("destination", "trajectory"), sum)
  		work.env$netflows = work.env$inflows - work.env$outflows
  
  		work.env$outflows.quantile = 
  		 	apply(work.env$outflows, MARGIN="origin", quantile, probs=probs, names=FALSE)
  		work.env$inflows.quantile = 
  			apply(work.env$inflows, MARGIN="destination", quantile, probs=probs, names=FALSE)
  		work.env$netflows.quantile = 
  			apply(work.env$netflows, MARGIN=1, quantile, probs=probs, names=FALSE)

  		work.env$globeMig = apply( work.env$od, MARGIN="trajectory", sum)
  		work.env$globeMig.quantile = quantile( work.env$globeMig, probs=probs, names=FALSE)
    
  		work.env$od.quantile = 
  			apply(work.env$od, MARGIN=c("origin", "destination"), quantile, probs=probs, names=FALSE)
  
  		work.env$tpop.quantile = 
  			apply(work.env$totp, MARGIN=1, quantile, probs=probs, names=FALSE)
  
  		quant.env$flowMarginPIarray[probs, country.codes.char, "outflow", t] = work.env$outflows.quantile[probs, country.codes.char] 
  		quant.env$flowMarginPIarray[probs, country.codes.char, "inflow", t] = work.env$inflows.quantile[probs, country.codes.char] 
  		quant.env$flowMarginPIarray[probs, country.codes.char, "netflow", t] = work.env$netflows.quantile[probs, country.codes.char] 
  
  		quant.env$flowPIarray[probs, country.codes.char, country.codes.char, t] = 
  			work.env$od.quantile[probs, country.codes.char, country.codes.char]

  		quant.env$popPIarray[probs, country.codes.char, t] = work.env$tpop.quantile[probs, country.codes.char]

  		quant.env$globeMigPIarray[,t] <- work.env$globeMig.quantile
  
  		rm(list=ls(work.env), envir=work.env)
	}
	return(quant.env)
}




.ini.pop.flow.res.env <- function(e, vital.events=FALSE, debug=FALSE) {
    if(debug) {
        print("before")
        print(c(address(e), refs(e)))
        if(exists('totpm', envir = e)) print(c("totpm: ", inspect("totpm", env = e)))
        if(exists('btm', envir = e)) print(c("btm: ", inspect("btm", env = e)))
    }
    for(item in c('totp', 'totpm', 'totpf', 'outmigm', 'outmigf', 'inmigm', 'inmigf', 'netmigm', 'netmigf')) 
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




get.pi.1trajectory <- function(country.code, itraj, dir){
	trajDestCodes = colnames( fread(file=paste0(dir, "/pi_", country.code, ".csv"), header=TRUE, nrows=1) )
	trajRow = fread(file=paste0(dir, "/pi_", country.code, ".csv"), header=FALSE, nrows=1, skip=itraj )
	return(list(pi.traj=trajRow, dest.codes=trajDestCodes))
}



general.outflow.bounds <- function(time) 
	c(0.8525257, 0.7541725, 0.6692184, 0.6216760, 0.5678551, 0.5216186)[min(time, 6)]

gcc.outflow.bounds <- function(time) 
	# "BHR" "KWT" "OMN" "QAT" "SAU" "ARE"
	c(0.6611600, 0.5974223, 0.5120738, 0.4597284, 0.4022156, 0.3522138)[min(time, 6)]

labor.outflow.bounds <- function(time) 
	# "BGD" "EGY" "IND" "IDN" "PAK" "PHL"
	c(0.9707848, 0.9500301, 0.9297360, 0.9129944, 0.8999324, 0.8854139)[min(time, 6)]

special.case.outflow.bounds <- function(time) 
	# "CHN" "IND" "IDN" "BRA"
	c(0.9923201, 0.9851209, 0.9801190, 0.9750746, 0.9716744, 0.9683950)[min(time, 6)]


is.gcc <- function(country)
	return(country %in% c(784, 48, 414, 512, 634, 682)) # c("BHR", "KWT", "OMN", "QAT", "SAU", "ARE")


is.labor <- function(country)
	return(country %in% c(50, 818, 360, 356, 586, 608)) # c("BGD", "EGY", "IND", "IDN", "PAK", "PHL")


is.special.case <- function(country)
	return(country %in% c(76, 156, 360, 356)) # c("CHN","IND", "IDN", "BRA")


is.incomplete.flow <- function(country){
	return(country %in% c(531, 499, 729, 688, 728)) # c("CUW", "MNE", "SDN", "SRB", "SSD")
}
