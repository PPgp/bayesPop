if(getRversion() >= "2.15.1") utils::globalVariables("UNlocations")

pop.aggregate <- function(pop.pred, regions, input.type=c('country', 'region'),
						name = input.type,
						inputs=list(e0F.sim.dir=NULL, e0M.sim.dir='joint_', tfr.sim.dir=NULL),
						my.location.file=NULL,
						verbose=FALSE, ...) {
	if(is.null(my.location.file))
		bayesTFR:::load.bdem.dataset('UNlocations', pop.pred$wpp.year, envir=globalenv(), verbose=verbose)
	else {
		env <- globalenv()
		if(is.character(my.location.file))
		    locs <- read.delim(my.location.file, comment.char='#')
		else locs <- my.location.file # data.frame
		assign("UNlocations", locs, envir=env)
	}
	regions <- unique(regions)
	method <- match.arg(input.type)
	if(missing(name)) name <- method
	if(method == 'country') 
		aggr.pred <- pop.aggregate.countries(pop.pred, regions, name, verbose=verbose, ...)
	if(method == 'region')
		aggr.pred <- pop.aggregate.regional(pop.pred, regions, name, inputs=inputs, verbose=verbose)
	invisible(get.pop.aggregation(pop.pred=pop.pred, name=name))
}

get.countries.for.region <- function(region, pop.pred) {
	reg.idx <- which(UNlocations[,'country_code'] == region)
	if(length(reg.idx) <= 0) {
		warning('No region code ', region, ' available.')
		return(c())
	}
	all.countries <- UNlocations[UNlocations[,'location_type'] == 4,]
	location.type <- UNlocations[reg.idx,'location_type']
	if(location.type == 0)  # corresponds to world, i.e. all countries
		countries <- all.countries[,'country_code']
	else {
		code.column <- switch(as.character(location.type), '2'='area_code', '3'='reg_code', paste0('agcode_', location.type))
		if(!is.element(code.column, colnames(UNlocations))) {
			warning('Invalid location type ', location.type, ' for region ', region,'. Location file must contain column ', code.column, '.')
			return(c())
		}	
		countries <- all.countries[is.element(all.countries[,code.column], region), 'country_code']
	}
	return(countries[is.element(countries, pop.pred$countries[,'code'])])
}

pop.aggregate.regional <- function(pop.pred, regions, name,
						inputs=list(e0F.sim.dir=NULL, e0M.sim.dir='joint_', tfr.sim.dir=NULL), 
						verbose=FALSE) {
	inp <- load.inputs(pop.pred$function.inputs, pop.pred$inputs$start.year, pop.pred$inputs$present.year, pop.pred$inputs$end.year, 
								pop.pred$wpp.year, verbose=verbose)
	for (item in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MIGm', 'MIGf', 'SRB', 'PASFR', 'MIGtype'))
		inp[[item]] <- NULL
	aggregated.countries <- list()
	aggr.obs.data <- list(male=NULL, female=NULL)
	if(verbose) cat('\nAggregating inputs using regional TFR and e0.')
	status.for.gui <- paste('out of', length(regions), 'regions.')
	gui.options <- list()
	obs.data <- inp$pop.matrix
	region.counter <- 0
	for(id in regions) {
		if(getOption('bDem.PopAgpred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			region.counter <- region.counter + 1
			gui.options$bDem.PopAgpred.status <- paste('finished', region.counter, status.for.gui)
			unblock.gtk('bDem.PopAgpred', gui.options)
		}

		countries <- get.countries.for.region(id, pop.pred)
		if(length(countries)==0) next
		inipop <- .aggregate.initial.pop(pop.pred, countries, id)
		for(item in c('POPm0', 'POPf0')) inp[[item]] <- rbind(inp[[item]], inipop[[item]])
		mort <- .aggregate.mortality.rates(pop.pred, countries, id)
		for(item in c('MXm', 'MXf')) inp[[item]] <- rbind(inp[[item]], mort[[item]])
		mig <- .aggregate.migration(pop.pred, countries, id)
		for(item in c('MIGm', 'MIGf')) inp[[item]] <- rbind(inp[[item]], mig[[item]])
		countries.index <- which(is.element(pop.pred$countries[,'code'], countries))
		inp$SRB <- rbind(inp$SRB, .aggregate.srb(pop.pred, countries, countries.index, id))
		if(inp$fixed.pasfr)
		    inp$PASFR <- rbind(inp$PASFR, .aggregate.pasfr(pop.pred, countries, countries.index, id))
		inp$MIGtype <- rbind(inp$MIGtype, .aggregate.migtype(pop.pred, countries, countries.index, id))
		aggregated.countries[[as.character(id)]] <- countries
		# aggregate pop.matrix
		for(cidx in 1:length(countries.index)) {
			country.obs.idx <- grep(paste('^', countries[cidx], '_', sep=''), rownames(obs.data[['male']]), value=FALSE)
			if(cidx == 1) {
				aggr.obs.dataM <- obs.data[['male']][country.obs.idx,]
				aggr.obs.dataF <- obs.data[['female']][country.obs.idx,]
				rownames(aggr.obs.dataM) <- rownames(aggr.obs.dataF) <- sub(paste(countries[cidx], '_', sep=''), 
													paste(id, '_', sep=''), rownames(obs.data[['male']][country.obs.idx,]))
			} else {
				aggr.obs.dataM <- aggr.obs.dataM + obs.data[['male']][country.obs.idx,]
				aggr.obs.dataF <- aggr.obs.dataF + obs.data[['female']][country.obs.idx,]
			}
		}
		aggr.obs.data[['male']] <- rbind(aggr.obs.data[['male']], aggr.obs.dataM)
		aggr.obs.data[['female']] <- rbind(aggr.obs.data[['female']], aggr.obs.dataF)
	}
	dir <- if(is.null(inputs$tfr.sim.dir)) pop.pred$function.inputs$tfr.sim.dir else inputs$tfr.sim.dir
	if(!is.null(dir)) inp$TFRpred <- get.tfr.prediction(dir)
	dir <- if(is.null(inputs$e0F.sim.dir)) pop.pred$function.inputs$e0F.sim.dir else inputs$e0F.sim.dir
	if(!is.null(dir)) inp$e0Fpred <- get.e0.prediction(dir)
	if((is.null(inputs$e0M.sim.dir) || inputs$e0M.sim.dir == 'joint_') && has.e0.jmale.prediction(inp$e0Fpred)) {
		inp$e0Mpred <- get.e0.jmale.prediction(inp$e0Fpred)
	} else inp$e0Mpred <- get.e0.prediction(if(is.null(inputs$e0M.sim.dir)) pop.pred$function.inputs$e0M.sim.dir else inputs$e0M.sim.dir)
	outdir <- gsub('predictions', paste('aggregations', name, sep='_'), pop.output.directory(pop.pred))
	if(file.exists(outdir)) unlink(outdir, recursive=TRUE)
	dir.create(outdir, recursive=TRUE)
	aggr.pred <- do.pop.predict(regions, inp=inp, outdir=outdir, nr.traj=pop.pred$nr.traj, ages=pop.pred$ages, verbose=verbose)
	aggr.pred <- .cleanup.pop.before.save(aggr.pred, remove.cache=TRUE)
	aggr.pred$aggregation.method <- 'region'
	aggr.pred$aggregated.countries <- aggregated.countries
	aggr.pred$inputs$pop.matrix <- aggr.obs.data
	bayesPop.prediction <- aggr.pred
	save(bayesPop.prediction, file=file.path(outdir, 'prediction.rda'))
	return(bayesPop.prediction)
}

.aggregate.by.sum <- function(pop.pred, countries, what.names, id) {
	.sum.fcn <- function(age, index, what) {
		age.idx <- index & (pop.pred$inputs[[what]][,'age'] == age)
		return(colSums(pop.pred$inputs[[what]][age.idx,3:ncol.inp, drop=FALSE]))
	} 
	l <- 21	
	all.ages <- seq(0, by=5, length=l)
	age.labels <- c(paste(all.ages[1:(l-1)], '-', all.ages[2:l]-1, sep=''), paste(all.ages[l],'+',sep=''))
	res <- list()
	for(item in what.names) {
		ncol.inp <- ncol(pop.pred$inputs[[item]])
		idx <- is.element(pop.pred$inputs[[item]][,'country_code'], countries)
		values <- lapply(age.labels, .sum.fcn, index=idx, what=item)
		values <- matrix(unlist(values), nrow=length(values), byrow=TRUE)
		colnames(values) <- colnames(pop.pred$inputs[[item]])[3:ncol.inp]
		res[[item]] <- data.frame(country_code=rep(id,l), age=age.labels, values, check.names=FALSE)
	}
	return(res)
}

.aggregate.initial.pop <- function(pop.pred, countries, aggr.id) 
	return(.aggregate.by.sum(pop.pred, countries, c('POPm0', 'POPf0'), aggr.id))

.aggregate.migration <- function(pop.pred, countries, aggr.id) 
	return(.aggregate.by.sum(pop.pred, countries, c('MIGm', 'MIGf'), aggr.id))

.aggregate.mortality.rates <- function(pop.pred, countries, aggr.id) {
	popmatrix.all <- list(MXm=pop.pred$inputs$pop.matrix$male, MXf=pop.pred$inputs$pop.matrix$female)
	res <- list()
	for(item in c('MXm', 'MXf')) {
		mort.idx <- is.element(pop.pred$inputs[[item]][,'country_code'], countries)
		popmatrix <- popmatrix.all[[item]]
		pop.countries.ages <- matrix(unlist(strsplit(rownames(popmatrix), '_')), ncol=2, byrow=TRUE)
		pop.idx <- is.element(pop.countries.ages[,1], countries)
		years.mx <- unlist(strsplit(colnames(pop.pred$inputs[[item]])[3:ncol(pop.pred$inputs[[item]])], '-'))
		pop.colidx <- which(is.element(as.integer(colnames(popmatrix)), years.mx[seq(2,length(years.mx),by=2)]))
		res[[item]] <- as.data.frame(matrix(NA, nrow=22, ncol=ncol(pop.pred$inputs[[item]]), 
						dimnames=list(c(0, 1, seq(5, 95, by=5), '100+'), colnames(pop.pred$inputs[[item]]))))
		trim.age <- gsub(' ', '', pop.pred$inputs[[item]][,'age'])
		trim.age.unique <- unique(trim.age)
		for(age in rownames(res[[item]])) {
			mort.age <- if(age == '100+') '100' else age # in wpp2012 there is no '100+' (in the new one there is) [&& !(age %in%  trim.age.unique)]
			mort.age.idx <- mort.idx & trim.age == mort.age
			pop.age.idx <- rep(0, nrow(pop.countries.ages))
			if(age == "0" || age =="1") pattern <- '^0-4'
			else {
				if(age == "100+") pattern <- age
				else pattern <- paste('^', age, '-', sep='')
			}
			pop.age.idx[grep(pattern, pop.countries.ages[,2], value=FALSE)] <- TRUE
			pop.age.idx <- pop.idx & pop.age.idx
			#if(sum(pop.age.idx) == 0 || sum(mort.age.idx) == 0) next
			person.years <- popmatrix[pop.age.idx, pop.colidx] # doesn't need to be multiplied by years because it cancels out in the ratio below
			deaths <- colSums(pop.pred$inputs[[item]][mort.age.idx,3:ncol(pop.pred$inputs[[item]])] * person.years)
			res[[item]][age,3:ncol(pop.pred$inputs[[item]])] <- deaths/colSums(person.years)
		}
		res[[item]][,"age"] <- rownames(res[[item]])
		res[[item]][,"country_code"] <- aggr.id
	}
	return(res)
}

.aggregate.srb <- function(pop.pred, countries, countries.index, aggr.id) {
	pred.meansM <- pop.pred$quantilesMage[countries.index,1, "0.5", 2:dim(pop.pred$quantilesMage)[4]]
	pred.meansF <- pop.pred$quantilesFage[countries.index,1, "0.5", 2:dim(pop.pred$quantilesFage)[4]]
	aggr.srb <- data.frame(country_code=aggr.id, matrix(colSums(pred.meansM)/colSums(pred.meansF), nrow=1))
	colnames(aggr.srb) <- colnames(pop.pred$inputs[["SRB"]])
	return(aggr.srb)
}

.aggregate.pasfr <- function(pop.pred, countries, countries.index, aggr.id) {
	pred.meansF <- pop.pred$quantilesFage[countries.index,4:11, "0.5", 2:dim(pop.pred$quantilesFage)[4]]
	idx <- is.element(pop.pred$inputs[["PASFR"]][,'country_code'], countries)
	res.ncol <- ncol(pop.pred$inputs[["PASFR"]])
	res <- as.data.frame(matrix(NA, nrow=7, ncol=res.ncol, 
					dimnames=list(NULL, colnames(pop.pred$inputs[["PASFR"]]))))
	ages <- pop.pred$inputs[["PASFR"]][1:7,"age"]
	res[,1] <- aggr.id
	res[,2] <- ages
	for (iage in 1:7) {
		age.idx <- idx & pop.pred$inputs[["PASFR"]][,'age'] == ages[iage]
		age.pasfr <- pop.pred$inputs[["PASFR"]][age.idx,3:res.ncol]
		res[iage,3:res.ncol] <- colSums(age.pasfr * pred.meansF[,iage,])/colSums(pred.meansF[,iage,])
	}
	return(res)
}

.aggregate.migtype <- function(pop.pred, countries, countries.index, aggr.id) {
	idx <- is.element(pop.pred$inputs[["MIGtype"]][,'country_code'], countries)
	obsN <- round(pop.pred$quantiles[countries.index, "0.5", 1],0) # observed in present year
	res <- matrix(NA, nrow=1, ncol=3, dimnames=list(NULL, colnames(pop.pred$inputs[["MIGtype"]])))
	res[,1] <- aggr.id
	res[,2] <- max(pop.pred$inputs[["MIGtype"]][idx,"ProjFirstYear"])
	f <- factor(rep(pop.pred$inputs[["MIGtype"]][idx,"MigCode"], obsN), levels=c(0,9))
	res[,3] <- as.integer(f[which.max(tapply(rep(1,sum(obsN)), f, sum))])
	return(res)
}

split.pop05 <- function(dat) {
    # split first pop age group into 0-1 & 1-4
    res <- abind(array(NA, dim(dat)[2:3]), dat, along = 1) # add row for the first age group (0-1)
    # using the sprague formula for P_0
    res[1,,] <- apply(dat[1,,, drop = FALSE]*0.3616 + dat[2,,, drop = FALSE]*-0.2768 + dat[3,,, drop = FALSE]*0.1488 + dat[4,,, drop = FALSE]*-0.0336, c(2,3), sum)
    res[2,,] <- dat[1,,] - res[1,,]
    if(any(res[1:2,,] < 0)) { # negative values, use spline method
        age5 <- seq(0, length = 7, by = 5)
        age1 <- seq(0, max(age5)+10, 1)
        for(time in 1:dim(res)[2]){
            if(any(res[1:2,time,] < 0)) {
                idx <- which(res[1,time,] < 0 | res[2,time,] < 0)
                for(i in idx){
                    pop5cum <- cumsum(c(0, dat[1:length(age5),time,i])) 
                    ## compute a monotone Hermite spline according to the method of Fritsch and Carlson implemented in splinefun {stats} in method "monoH.FC" 
                    yinterpol.fun <- splinefun(c(age5, max(age5)+10), pop5cum, method = "monoH.FC")
                    ## predict interpolated values by single age on cumulative function
                    pop1cum <- yinterpol.fun(age1)
                    ## decumulate
                    pop1 <- diff(pop1cum) 
                    res[1,time,i] <- pop1[1]
                }
                #res[1,time,idx] <- dat[1,time,idx]/5.
                res[2,time,idx] <- dat[1,time,idx] - res[1,time,idx]
            }
        }
    }
    return(res)
}

pop.aggregate.countries <- function(pop.pred, regions, name, 
                                    use.kannisto = TRUE, keep.vital.events = NULL,
                                    verbose=verbose, adjust=FALSE) {
    
	if(verbose) cat('\nAggregating using countries as inputs.')
	nreg <- length(regions)
	quantiles.to.keep <- as.numeric(dimnames(pop.pred$quantiles)[[2]])
	quant <- quantM <- quantF <- array(NA, c(nreg, dim(pop.pred$quantiles)[2:3]), dimnames=c(list(regions), dimnames(pop.pred$quantiles)[2:3]))
	quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(nreg, dim(pop.pred$quantilesMage)[2:4]),
						dimnames=c(list(regions), dimnames(pop.pred$quantilesMage)[2:4]))
	mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(nreg,dim(pop.pred$traj.mean.sd)[2:3]))
	obs.data <- pop.pred$inputs$pop.matrix
	# convert to arrays for performance reasons
	for(sex in names(obs.data)) 
	    obs.data[[sex]] <- data.matrix(obs.data[[sex]])

	aggr.obs.data <- aggrobs <- abrpop <- abrpopobs <- aggrabrpop <- aggrabrpopobs <- abrd <- aggrabrdeaths <- abrdobs <- aggrabrdeathsobs <- abrpop.hch <- aggrabrpop.hch <- abrd.hch <- aggrabrdeaths.hch <- list(male=NULL, female=NULL)
	outdir <- gsub('predictions', paste('aggregations', name, sep='_'), pop.output.directory(pop.pred))
	if(file.exists(outdir)) unlink(outdir, recursive=TRUE)
	dir.create(outdir, recursive=TRUE)
	aggregated.countries <- list()
	id.idx <- 0
	valid.regions <- rep(FALSE, length(regions))
	status.for.gui <- paste('out of', nreg, 'regions.')
	gui.options <- list()
	no.vital.events <- !is.null(keep.vital.events) && keep.vital.events == FALSE
	aggr.quantities <- c('totp', 'totpm', 'totpf', 'totp.hch', 'totpm.hch', 'totpf.hch')
	aggr.quantities.ve <- c('btm', 'btf', 'btm.hch', 'btf.hch', 
								'deathsm', 'deathsf', 'deathsm.hch', 'deathsf.hch',
								'migm', 'migf')
	ind.names <- list(male = list(pop = "totpm", mx = "mxm"), 
	                 female = list(pop = "totpf", mx = "mxf"))
	ind.names.hch <- list(male = list(pop = "totpm.hch", mx = "mxm.hch"), 
	                  female = list(pop = "totpf.hch", mx = "mxf.hch"))
	# The next two lines need to be there for the R checker to get to know these objects
	totp <- totpm <- totpf <- totp.hch <- totpm.hch <- totpf.hch <- NULL
	btm <- btf <- deathsm <- deathsf <- migm <- migf <- btm.hch <- btf.hch <- deathsm.hch <- deathsf.hch <- NULL
	aggr.quantities.all <- aggr.quantities
	max.lage <- dim(pop.pred$quantilesMage)[2]
	max.lage.fert <- age.length.fert(pop.pred$annual)
	fert.age.start <- age.index.fert(pop.pred$annual)[1]
	#mort.ages <- lt.ages(pop.pred$annual)
	if(pop.pred$annual) {
	    time.step <- 1
	    mort.ages <- seq(0, length = max.lage)
	    kannisto.proj.ages <- 100:130
	    kannisto.est.ages <- 90:99
	} else { # 5x5
	    time.step <- 5
	    mort.ages <- c(0,1,seq(5, length = max.lage-1, by = 5))
	    kannisto.proj.ages <- seq(100, 130, by = 5)
	    kannisto.est.ages <- seq(80, 95, by = 5)
	}
	fert.age.idx <- seq(fert.age.start, length = max.lage.fert)
	obs.cols.ve <- colnames(obs.data[["male"]])
	if(!pop.pred$annual) obs.cols.ve <- paste(as.integer(colnames(obs.data[["male"]]))-5, obs.cols.ve, sep = "-")
	prev.year <- as.character(pop.pred$proj.years.pop[1]-time.step)
	
	for(reg.idx in 1:length(regions)) {
		if(getOption('bDem.PopAgpred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			gui.options$bDem.PopAgpred.status <- paste('finished', reg.idx, status.for.gui)
			unblock.gtk('bDem.PopAgpred', gui.options)
		}
		id <- regions[reg.idx]
		if(verbose) cat('\nAggregating region ', id)
		countries <- get.countries.for.region(id, pop.pred)
		if(length(countries)==0) next
		id.idx <- id.idx + 1
		valid.regions[reg.idx] <- TRUE
		countries.index <- which(is.element(pop.pred$countries[,'code'], countries))
		e <- new.env()
		if(adjust && is.null(pop.pred$adjust.env)) pop.pred$adjust.env <- new.env()
		
		#debug.df <- NULL
		for(cidx in 1:length(countries.index)) {
			country.obs.idx <- grep(paste0('^', countries[cidx], '_'), rownames(obs.data[['male']]), value=FALSE)
			traj.file <- file.path(pop.output.directory(pop.pred), paste('totpop_country', countries[cidx], '.rda', sep=''))
			load(traj.file, envir=e)
			if(adjust) adjust.trajectories(countries[cidx], e, pop.pred, pop.pred$adjust.env)
			if(!no.vital.events) {
				ve.file <- file.path(pop.output.directory(pop.pred), 
									paste('vital_events_country', countries[cidx], '.rda', sep=''))
				if(file.exists(ve.file)) {
					if(cidx == 1) {
					    #requireNamespace("DemoTools")
						aggr.quantities.all <- c(aggr.quantities, aggr.quantities.ve)
						observed <- list()
					}
					load(ve.file, envir=e)
					use.obs.cols <- intersect(obs.cols.ve, colnames(e$observed$mxm))
					use.obs.cols.pop <- if(pop.pred$annual) use.obs.cols else substr(use.obs.cols, 6, 9)
					
					# compute abridged population and deaths
					for(sex in c("male", "female")) {
					   # Attach the second last observed year of population
					    popavg <- array(0, c(max.lage, dim(e[[ind.names[[sex]]$pop]])[2]+1, dim(e[[ind.names[[sex]]$pop]])[3]))
					    popavg[1:length(country.obs.idx),1,] <- obs.data[[sex]][country.obs.idx, prev.year]
					    popavg[,-1,] <- e[[ind.names[[sex]]$pop]]
					    popavg <- (popavg[,-1, ,drop = FALSE] + popavg[,-dim(popavg)[2],, drop = FALSE])/2.
					    popavg.hch <- array(0, c(max.lage, dim(e[[ind.names[[sex]]$pop]])[2]+1, 2))
					    popavg.hch[1:length(country.obs.idx),1,] <- obs.data[[sex]][country.obs.idx, prev.year]
					    popavg.hch[,-1,] <- e[[ind.names.hch[[sex]]$pop]]
					    popavg.hch <- (popavg.hch[,-1, ,drop = FALSE] + popavg.hch[,-dim(popavg)[2],, drop = FALSE])/2.
					    popavgobs <- (obs.data[[sex]][country.obs.idx, -1 , drop=FALSE] + 
					                      obs.data[[sex]][country.obs.idx, -dim(obs.data[[sex]])[2] , drop=FALSE])/2.
					    if(!pop.pred$annual){
					        # projected abridged pop and deaths
					        abrpop[[sex]] <- split.pop05(popavg)
					        abrpop.hch[[sex]] <- split.pop05(popavg.hch)
					        abrpopobs[[sex]] <- split.pop05(abind(popavgobs[,use.obs.cols.pop, drop = FALSE], NULL, along = 3))
                        } else {
                            abrpop[[sex]] <- popavg
                            abrpop.hch[[sex]] <- popavg.hch
                            abrpopobs[[sex]] <- abind(popavgobs[,use.obs.cols.pop, drop = FALSE], NULL, along = 3)
                        }
					    abrd[[sex]] <- abrpop[[sex]] * e[[ind.names[[sex]]$mx]] # should be multiplied by five but that would get canceled out later
					    abrd.hch[[sex]] <- abrpop.hch[[sex]] * e[[ind.names.hch[[sex]]$mx]]
					    abrdobs[[sex]] <- abrpopobs[[sex]]*e$observed[[ind.names[[sex]]$mx]][,use.obs.cols,, drop = FALSE]
					}
					# debug.df <- rbind(debug.df, data.frame(
					#                     country = countries[cidx],
					#                     age = c(0,1, seq(5, 130, by = 5)),
					#                     #popMabr_2015 = c(spltmpM[,"2010-2015",], rep(NA, 6)),
					#                   #mxM_2010_15 = c(e$observed$mxm[,"2010-2015",], rep(NA, 6)),
					#                   #deathsM_2010_15 = c(abrdobs[["male"]][,"2010-2015"], rep(NA, 6)),
					#                   popM_2015 = c(NA, obs.data[["male"]][country.obs.idx,"2015",], rep(NA, 6)),
					#                   popM_2020 = c(NA, e$totpm[,"2020",]),
					#                   popMabr_2020 = split.pop05(e$totpm)[,"2020",], 
					#                   mxM_2015_20 = e$mxm[,"2018",],
					#                   deathsM_2015_20 = abrd[["male"]][,"2020",],
					#                   popM_2025 = c(NA, e$totpm[,"2025",]),
					#                   popMabr_2025 = split.pop05(e$totpm)[,"2025",],
					#                   mxM_2020_25 = e$mxm[,"2023",],
					#                   deathsM_2020_25 = abrd[["male"]][,"2025",])
					# )
				} else {
				    no.vital.events <- TRUE
				}
			}
			#debug.dt <- data.table(debug.df)
			#debug.dt[]
			if(cidx == 1) { # first country to be processed
				for(par in aggr.quantities.all)
					assign(par, e[[par]])
			    for(sex in names(aggrobs)) 
				    aggrobs[[sex]] <- obs.data[[sex]][country.obs.idx,, drop=FALSE]
				rownames(aggrobs[["male"]]) <- rownames(aggrobs[["female"]]) <- sub(paste0(countries[cidx], '_'), 
																paste0(id, '_'), rownames(obs.data[['male']][country.obs.idx,, drop=FALSE]))
				if(!no.vital.events) {
					for(par in aggr.quantities.ve)
						if(!is.null(e$observed[[par]]))
							observed[[par]] <- e$observed[[par]]
					for(sex in names(aggrabrdeaths)) {
					    aggrabrdeaths[[sex]] <- abrd[[sex]]
					    aggrabrdeathsobs[[sex]] <- abrdobs[[sex]]
					    aggrabrdeaths.hch[[sex]] <- abrd.hch[[sex]]
					    aggrabrpop[[sex]] <- abrpop[[sex]]
					    aggrabrpop.hch[[sex]] <- abrpop.hch[[sex]]
					    aggrabrpopobs[[sex]] <- abrpopobs[[sex]]
					}
				}
				trajectory.indices <- e$trajectory.indices											
			} else { # at least one country was already processed
				for(par in aggr.quantities.all) 
					assign(par, get(par) + e[[par]])
			    for(sex in names(aggrobs)) 
			        aggrobs[[sex]] <- aggrobs[[sex]] + obs.data[[sex]][country.obs.idx,]
				if(!no.vital.events) {
					for(par in names(observed))
						observed[[par]] <- observed[[par]] + e$observed[[par]]
					for(sex in names(aggrabrdeaths)) {
					    aggrabrpop[[sex]] <- aggrabrpop[[sex]] + abrpop[[sex]]
					    aggrabrpopobs[[sex]] <- aggrabrpopobs[[sex]] + abrpopobs[[sex]]
					    aggrabrpop.hch[[sex]] <- aggrabrpop.hch[[sex]] + abrpop.hch[[sex]]
					    aggrabrdeaths[[sex]] <- aggrabrdeaths[[sex]] + abrd[[sex]]
					    aggrabrdeaths.hch[[sex]] <- aggrabrdeaths.hch[[sex]] + abrd.hch[[sex]]
					    aggrabrdeathsobs[[sex]] <- aggrabrdeathsobs[[sex]] + abrdobs[[sex]]
					}
				}
			}
		}
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch, trajectory.indices,
			 file = file.path(outdir, paste0('totpop_country', id, '.rda')))
		if(!no.vital.events) {
		    popfwprev <- array(0, c(max.lage.fert, dim(totpf)[2]+1, dim(totpf)[3]))
		    popfwprev[,1,] <- aggrobs[["female"]][fert.age.idx, prev.year]
		    popfwprev[,-1,] <- totpf[fert.age.idx,,]
		    # asfert
		    asfert <- 2*(btm + btf)/(popfwprev[,-dim(popfwprev)[2],,drop=FALSE]+
		                                 popfwprev[,-1,,drop=FALSE])
		    popfwprev.hch <- array(0, c(max.lage.fert, dim(totpf)[2]+1, dim(totpf.hch)[3]))
		    popfwprev.hch[,1,] <- aggrobs[["female"]][fert.age.idx, prev.year]
		    popfwprev.hch[,-1,] <- totpf.hch[fert.age.idx,,]
		    asfert.hch <- 2*(btm.hch + btf.hch)/(popfwprev.hch[,-dim(popfwprev.hch)[2],,drop=FALSE]+
		                                             popfwprev.hch[,-1,,drop=FALSE])
		    # pasfert
		    tfr <- apply(asfert, c(2,3), sum)
		    pasfert <- asfert/abind(tfr, NULL, along=0)[rep(1,dim(asfert)[1]),,,drop=FALSE]*100
		    tfr.hch <- apply(asfert.hch, c(2,3), sum)
		    pasfert.hch <- asfert.hch/abind(tfr.hch, NULL, along=0)[rep(1,dim(asfert.hch)[1]),,,drop=FALSE]*100
		    
		    # compute mortality rates from aggregated abridged deaths and population
		    mxm <- aggrabrdeaths[["male"]] / aggrabrpop[["male"]] 
		    mxf <- aggrabrdeaths[["female"]] / aggrabrpop[["female"]]
		    mxm.hch <- aggrabrdeaths.hch[["male"]] / aggrabrpop.hch[["male"]] 
		    mxf.hch <- aggrabrdeaths.hch[["female"]] / aggrabrpop.hch[["female"]]
		    dimnames(mxm)[1:2] <- dimnames(mxf)[1:2] <- dimnames(mxm.hch)[1:2] <- dimnames(mxf.hch)[1:2] <- list(
		        mort.ages, dimnames(btm)[[2]])
		    
			# apply kannisto for old age groups
			if(use.kannisto) {
			    for(itraj in 1:dim(mxm)[3]) {
			        kan <- do.call("cokannisto", list(mxm[1:max.lage,,itraj], mxf[1:max.lage,,itraj], 
			                       est.ages = kannisto.est.ages, proj.ages = kannisto.proj.ages))
			        mxm[,,itraj] <- kan$male
			        mxf[,,itraj] <- kan$female
			    }
			    for(ivar in 1:dim(mxm.hch)[3]) {
			        kan.hch <- do.call("cokannisto", list(mxm.hch[1:max.lage,,ivar], mxf.hch[1:max.lage,,ivar],
			                                              est.ages = kannisto.est.ages, proj.ages = kannisto.proj.ages))
			        mxm.hch[,,ivar] <- kan.hch$male
			        mxf.hch[,,ivar] <- kan.hch$female
			    }
			}
			# asfert, pasfert for observed data
			observed <- within(observed, {
			    tmp <- abind(aggrobs[["female"]][fert.age.idx, , drop=FALSE], NULL, along=3)
				if(dim(tmp)[2] > dim(btf)[2]+1) # if dimension of births doesn't match population
					tmp <- tmp[,-(1:(dim(tmp)[2]-dim(btf)[2]-1)),, drop=FALSE]
				asfert <- 2*(btm + btf)/(tmp[,-dim(tmp)[2],,drop=FALSE] + tmp[,-1,,drop=FALSE])
				tfr <- apply(asfert, c(2,3), sum)
				pasfert <- asfert/abind(tfr, NULL, along=0)[rep(1,dim(asfert)[1]),,,drop=FALSE]*100
				rm(tmp, tfr)
				# mortality
				mxm <- aggrabrdeathsobs[["male"]] / aggrabrpopobs[["male"]] 
				mxf <- aggrabrdeathsobs[["female"]] / aggrabrpopobs[["female"]]
				dimnames(mxm)[1:2] <- dimnames(mxf)[1:2] <- list(mort.ages[1:dim(mxm)[1]], dimnames(btm)[[2]])
			})
			save(btm, btf, deathsm, deathsf, mxm, mxf, migm, migf, asfert, pasfert,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, pasfert.hch, mxm.hch, mxf.hch,
				observed, file=file.path(outdir, paste0('vital_events_country', id, '.rda')))
		}
		#stop("")
		quant[id.idx,,] = apply(totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sd[id.idx,1,] <- apply(totp, 1, mean, na.rm = TRUE)
		mean_sd[id.idx,2,] <- apply(totp, 1, sd, na.rm = TRUE)
		
		for (i in 1:dim(totpm)[1]) {
			if(dim(totpm)[3] == 1) { # 1 trajectory
				quantMage[id.idx,i,,] <- matrix(totpm[i,,1], nrow=dim(quantMage)[3], ncol=dim(quantMage)[4], byrow=TRUE)
				quantFage[id.idx,i,,] <- matrix(totpf[i,,1], nrow=dim(quantFage)[3], ncol=dim(quantFage)[4], byrow=TRUE)
				quantPropMage[id.idx,i,,] <- matrix(totpm[i,,1]/totp[,1], nrow=dim(quantPropMage)[3], 
													ncol=dim(quantPropMage)[4], byrow=TRUE)
				
				quantPropFage[id.idx,i,,] <- matrix(totpf[i,,1]/totp[,1], nrow=dim(quantPropFage)[3], 
													ncol=dim(quantPropFage)[4], byrow=TRUE)
			} else { # multiple trajectories
				quantMage[id.idx,i,,] <- apply(totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantFage[id.idx,i,,] <- apply(totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropMage[id.idx,i,,] <- apply(totpm[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropFage[id.idx,i,,] <- apply(totpf[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			}
		}
		sstotpm <- colSums(totpm)
		quantM[id.idx,,] = apply(sstotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdM[id.idx,1,] <- apply(sstotpm, 1, mean, na.rm = TRUE)
		mean_sdM[id.idx,2,] = apply(sstotpm, 1, sd, na.rm = TRUE)
		sstotpf <- colSums(totpf)
		quantF[id.idx,,] = apply(sstotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdF[id.idx,1,] <- apply(sstotpf, 1, mean, na.rm = TRUE)
		mean_sdF[id.idx,2,] = apply(sstotpf, 1, sd, na.rm = TRUE)
		aggregated.countries[[as.character(id)]] <- countries
		for(sex in names(aggrobs)) 
		    aggr.obs.data[[sex]] <- rbind(aggr.obs.data[[sex]], aggrobs[[sex]])
	}
	aggr.pred <- pop.pred
	which.reg.index <- function(x, set) return(which(set == x))
	reg.idx <- sapply(regions[valid.regions], which.reg.index, set=UNlocations[,'country_code']) 
	aggr.pred$countries <- data.frame(code=UNlocations[reg.idx, 'country_code'], name=UNlocations[reg.idx, 'name'])
	aggr.pred$quantiles <- quant
	aggr.pred$quantilesM <- quantM
	aggr.pred$quantilesF <- quantF
	aggr.pred$quantilesMage <- quantMage
	aggr.pred$quantilesFage <- quantFage
	aggr.pred$quantilesPropMage <- quantPropMage
	aggr.pred$quantilesPropFage <- quantPropFage
	aggr.pred$traj.mean.sd <- mean_sd
	aggr.pred$traj.mean.sdM <- mean_sdM
	aggr.pred$traj.mean.sdF <- mean_sdF
	aggr.pred$aggregation.method <- 'country'
	aggr.pred$aggregated.countries <- aggregated.countries
	aggr.pred$inputs <- as.environment(as.list(pop.pred$inputs, all.names=TRUE)) # clone environment
	aggr.pred$inputs$pop.matrix <- aggr.obs.data
	bayesPop.prediction <- .cleanup.pop.before.save(aggr.pred, remove.cache=TRUE)
	save(bayesPop.prediction, file=file.path(outdir, 'prediction.rda'))
	cat('\nAggregations stored into', outdir, '\n')
	return(bayesPop.prediction)
}