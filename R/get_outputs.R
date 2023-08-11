
get.expression.indicators <- function() {
		return(list(D='deaths', B='births', S='survival', F='fertility', Q='qx', M='mx', 
		            G='migration', R='pasfr', E='ex', A='ax'))
}

age.index.all <- function(...) 
    return(1:age.length.all(...))

age.length.all <- function(annual = FALSE, observed = FALSE){
    if(annual) {
        if(observed) return(101) 
        return(131)
    }
    if(observed) return(21) 
    return(27)
}


ages.all <- function(annual = FALSE, observed = FALSE) {
    l <- age.length.all(annual, observed)
    if(annual) 
        return(seq(0, length = l))
    return(seq(0, by = 5, length = l))
}

age.index.fert <- function(annual = FALSE) {
    if(annual) return(11:55)
    return(4:10)
}

age.length.fert <- function(...)
    return(length(age.index.fert(...)))

ages.fert <- function(...)
    return(ages.all(...)[age.index.fert(...)])

lt.age.length <- function(annual = FALSE, ...) {
    l <- age.length.all(annual = annual, ...)
    if(!annual) l <- l + 1
    return(l)
}

lt.ages <- function(annual = FALSE, observed = FALSE) {
    l <- lt.age.length(annual, observed)
    if(annual) 
        return(seq(0, length = l))
    return(c(0, 1, seq(5, by = 5, length = l-2)))
}

lt.age.index <- function(...) 
    return(1:lt.age.length(...))


has.pop.prediction <- function(sim.dir) {
	if(file.exists(file.path(sim.dir, 'predictions', 'prediction.rda'))) return(TRUE)
	return(FALSE)
}

pop.output.directory <- function(pop.pred) return(file.path(pop.pred$base.directory, pop.pred$output.directory))

get.pop.prediction <- function(sim.dir, aggregation=NULL, write.to.cache=TRUE) {
	############
	# Returns an object of class bayesPop.prediction
	############
	if(!is.null(aggregation)) return(get.pop.aggregation(sim.dir, name=aggregation))
	output.dir <- file.path(sim.dir, 'predictions')
	pop.pred <- .get.prediction.object(output.dir, 'predictions')
	pop.pred$base.directory <- normalizePath(sim.dir)
	pop.pred$cache <- .load.cache(output.dir)
	pop.pred$write.to.cache <- write.to.cache
	pop.pred$is.aggregation <- FALSE
	if(is.null(pop.pred$annual)) pop.pred$annual <- FALSE
	return(pop.pred)
}
.cleanup.pop.before.save <- function(pop.pred, remove.cache=FALSE) {
	# remove items that do not have to be saved because they are added using get.pop.prediction
	names.to.remove <- c('output.directory', 'base.directory', 'write.to.cache', 'is.aggregation', 
							if(remove.cache) 'cache' else NULL)
	names.to.remove <- names.to.remove[names.to.remove %in% names(pop.pred)]
	if(length(names.to.remove) > 0) {
		pred <- pop.pred[-which(names(pop.pred) %in% names.to.remove)] # this removes the class attribute
		class(pred) <- class(pop.pred)
		return(pred)
	}
	return(pop.pred)
}
.load.cache <- function(sim.dir) {
	if(!file.exists(file.path(sim.dir, 'cache.rda'))) return(new.env())
	cache <- local({load(file.path(sim.dir, 'cache.rda'))
					cache})
	return(as.environment(cache))	
}

.save.cache <- function(pop.pred) {
	if(is.null(pop.pred$cache) || (!is.null(pop.pred$write.to.cache) && !pop.pred$write.to.cache)) return()
	cache <- pop.pred$cache
	save(cache, file=file.path(pop.output.directory(pop.pred), 'cache.rda'))
}

.remove.cache.file <- function(dir) {
	file.name <-file.path(dir, 'cache.rda')
	if(file.exists(file.name)) unlink(file.name)
}

pop.cleanup.cache <- function(pop.pred) {
	if(!is.null(pop.pred$write.to.cache) && !pop.pred$write.to.cache) {
		if(!is.null(pop.pred$cache))
			warning('No cache manipulation allowed for this prediction object.')
		return()
	}
	.remove.cache.file(pop.output.directory(pop.pred))
	if(!is.null(pop.pred$cache))
		rm(list=ls(pop.pred$cache), envir=pop.pred$cache)
	gc()
	return()
}

has.pop.aggregation <- function(sim.dir=NULL, pop.pred=NULL, return.dirs=FALSE) {
	if(is.null(sim.dir)) sim.dir <- pop.pred$base.directory
	dirs <- list.files(sim.dir, pattern='^aggregations_', full.names=FALSE)
	if(return.dirs) return(dirs)
	return (length(dirs) > 0)
}

available.pop.aggregations <- function(pop.pred){
	dirs <- has.pop.aggregation(pop.pred=pop.pred, return.dirs=TRUE)
	if(length(dirs)<=0) return(c())
	return(substr(dirs, 14, nchar(dirs)))
}

get.pop.aggregation <- function(sim.dir=NULL, pop.pred=NULL, name=NULL, write.to.cache=TRUE) {
	############
	# Returns an object of class bayesPop.prediction created by aggregation
	############
	if(is.null(sim.dir)) sim.dir <- pop.pred$base.directory
	dirs <- has.pop.aggregation(sim.dir=sim.dir, return.dirs=TRUE)
	if(length(dirs) == 0) {
		warning('No aggregation available in', sim.dir)
		return(NULL)
	}
	output.dir <- file.path(sim.dir, dirs)
	names <- substr(dirs, 14, nchar(dirs))
	if(length(names) == 1){
		if(!is.null(name) && name != names) 
			warning('Mismatch in aggregation names. Available aggregation is called ', names)
		pop.aggr <- .get.prediction.object(output.dir, paste('aggregations', names, sep='_'))
	} else {
		idx <- which(names == name)
		if (length(idx) == 0) idx <- menu(names, title='Available aggregations:')	
		pop.aggr <- .get.prediction.object(output.dir[idx], paste('aggregations', names[idx], sep='_'))
	}
	pop.aggr$base.directory <- normalizePath(sim.dir)
	pop.aggr$is.aggregation <- TRUE
	pop.aggr$cache <- .load.cache(pop.output.directory(pop.aggr))
	pop.aggr$write.to.cache <- write.to.cache
	if(is.null(pop.aggr$annual)) pop.aggr$annual <- FALSE

	return(pop.aggr)	
}

.get.prediction.object <- function(directory, name=directory) {
	pred.file <- file.path(directory, 'prediction.rda')
	if(!file.exists(pred.file)) {
		warning('File ', pred.file, ' does not exist.')
		return(NULL)
	}
	load(file=pred.file)
	bayesPop.prediction$output.directory <- name
	# convert inputs to environment
	if(!is.environment(bayesPop.prediction$inputs))
		bayesPop.prediction$inputs <- list2env(bayesPop.prediction$inputs)
	return(bayesPop.prediction)
}

summary.bayesPop.prediction <- function(object, country=NULL, sex=c('both', 'male', 'female'), compact=TRUE, ...) {
	if(!is.null(country)) country <- get.country.object(country, country.table=object$countries)
	sex <- match.arg(sex)
	if (sex == 'male') {
		object$quantiles <- object$quantilesM
		object$traj.mean.sd <- object$traj.mean.sdM
	} else {
		if (sex == 'female') {
			object$quantiles <- object$quantilesF
			object$traj.mean.sd <- object$traj.mean.sdF
		} else sex <- 'both'
	}
	res <- bayesTFR:::get.prediction.summary.data(object, 
				unchanged.pars=c('nr.traj', 'estim.years'), 
				country=country, compact=compact)
	res$nr.countries <- nrow(object$countries)
	res$sex <- sex
	class(res) <- 'summary.bayesPop.prediction'
	return(res)
}

print.summary.bayesPop.prediction <- function(x, digits = 5, ...) {
	cat('\nProjections:', length(x$projection.years), '(', x$projection.years[1], '-', 
					x$projection.years[length(x$projection.years)], ')')
	cat('\nInitial time point:', x$estim.years[length(x$estim.years)])
	cat('\nObserved time points:', length(x$estim.years), '(', x$estim.years[1], '-', 
					x$estim.years[length(x$estim.years)], ')')
	cat('\nTrajectories:', x$nr.traj)
	cat('\nNumber of countries:', x$nr.countries)
	if(!is.null(x$country.name)) {
		cat('\n\nCountry:', x$country.name, '\n')
		cat('\nProjected Population')
		if (x$sex != 'both') cat(' for', x$sex)
		cat(':\n')
		print(x$projections, digits=digits, ...)
	}
}

get.pop.observed.with.age <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all', data=NULL, annual = FALSE) {
	# Results are not sorted in the same order as values in "country". The caller should take care of it. 
	sex <- match.arg(sex)
	if(is.null(data)) data <- pop.pred$inputs$pop.matrix
	if(sex == 'both') {
		data <- data[['male']][,colnames(data[['male']]),drop=FALSE] + data[['female']][,colnames(data[['male']]),drop=FALSE]
	} else data <- data[[sex]]
	country.idx <- grep(paste('^', country, '_', sep='', collapse='|'), rownames(data), value=FALSE)
	data <- data[country.idx,, drop=FALSE]
	if(is.null(pop.pred$proj.years.pop)) {
		coln <- as.integer(colnames(data))
		if(coln[1] %% 5 != 0 && !annual) # column names should be the end of 5-year interval (not the middle)
			colnames(data) <- as.integer(colnames(data)) + 2
	}
	max.age <- as.integer(round(nrow(data)/length(country),0))
	age.idx <- if(age[1]=='all') 1:max.age else age
	age.idx <- age.idx[age.idx <= max.age]
	return(list(data=data, age.idx=age.idx, max.age=max.age))
}


get.pop.observed <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all', sum.over.ages=TRUE) {
	data.age <- get.pop.observed.with.age(pop.pred, country, sex, age)
	data <- data.age$data
	if(nrow(data) == 0) {
		d <- rep(NA, ncol(data))
		names(d) <- names(data)
		return(d)
	}
	age.idx <- data.age$age.idx
	if(sum.over.ages) return(colSums(data[age.idx,,drop=FALSE]))
	return(data[age.idx,,drop=FALSE])
}

get.pop.observed.multiple.countries <- function(pop.pred, countries, sex=c('both', 'male', 'female'), age='all', sum.over.ages=TRUE) {
	# the function assumes that population data are sorted be age
	data.age <- get.pop.observed.with.age(pop.pred, countries, sex, age)
	data <- data.age$data
	age.idx <- data.age$age.idx
	ncountries <- length(countries)
	max.age <- data.age$max.age
	cindex <- lapply(countries, function(x) grep(paste0("^", x, "_"), rownames(data), value=FALSE)[1:max.age])
	names(cindex) <- countries
	sum.over.countries <- function(country.idx) return(colSums(data[country.idx,][age.idx,]))
	if(sum.over.ages) return(list(data=t(sapply(cindex, sum.over.countries)), age.idx=age.idx))
	res <- array(NA, c(ncountries, length(age.idx), ncol(data)), dimnames=list(countries, age.idx, colnames(data)))
	countries.char <- names(cindex)
	for(i in 1:ncountries) res[i,,] <- as.matrix(data[cindex[[countries.char[i]]],])
	return(list(data=res, age.idx=age.idx))
}

.get.pop.quantiles <- function(pop.pred, what='', adjust=FALSE, allow.negative.adj = TRUE) {
	quant <- pop.pred[[paste0('quantiles', what)]]
	if(!adjust) return(quant)
	return(adjust.quantiles(quant, what, wpp.year=pop.pred$wpp.year, annual = pop.pred$annual, 
	                        env=pop.pred$adjust.env, allow.negatives = allow.negative.adj))
}

.load.traj.file <- function(dir, country, e) {
	traj.file <- file.path(dir, paste('totpop_country', country, '.rda', sep=''))
	if (!file.exists(traj.file)) return(0)
	load(traj.file, envir=e)
	return(1)
}

pop.trajectories <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all', ...) {
	country.object <- get.country.object(country, country.table=pop.pred$countries)
	res <- get.pop.trajectories(pop.pred, country.object$code, sex=sex, age=age, ...)
	return(res$trajectories[,res$index])
}

get.pop.trajectories <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all',
 									nr.traj=NULL, typical.trajectory=FALSE, adjust=FALSE, allow.negative.adj = TRUE) {
	
	quant <- hch <- age.idx <- traj <- traj.idx <-  NULL
	load.traj <- is.null(nr.traj) || nr.traj > 0 || typical.trajectory || adjust
	e <- new.env()
	if (!.load.traj.file(pop.output.directory(pop.pred), country, e))
		return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch))
	if(adjust) {
		if(is.null(pop.pred$adjust.env)) pop.pred$adjust.env <- new.env()
		adjust.trajectories(country, e, pop.pred, pop.pred$adjust.env, allow.negatives = allow.negative.adj)
	}
	sex <- match.arg(sex)
	max.age <- dim(e$totpf)[1] # should be 27 or 131 if annual = TRUE
	age.idx <- if(age[1]=='all') 1:max.age else age
	annual <- pop.pred$annual
	max.age.allowed <- age.length.all(annual)
	if(max(age.idx) > max.age.allowed || min(age.idx) < 1) {
	    if(annual) stop(paste('Age index must be between 0 and ', max.age.allowed - 1, '.')) 
	    stop(paste('Age index must be between 1 and ', max.age.allowed, '(age 130+).'))
	}
	if(sex == 'both' && all((1:max.age) %in% age.idx)) { # for both sexes and all ages
		if(load.traj) traj <- e$totp
		quant <- .get.pop.quantiles(pop.pred, adjust=adjust, allow.negative.adj = allow.negative.adj)
		hch <- e$totp.hch
	} else {
	    if(sex == 'both') {
	        if(load.traj) traj <- colSums(e$totpm[age.idx,,,drop=FALSE]) + colSums(e$totpf[age.idx,,,drop=FALSE])
	        hch <- colSums(e$totpm.hch[age.idx,,,drop=FALSE]) + colSums(e$totpf.hch[age.idx,,,drop=FALSE])
	    } else {
	        if(sex=='male') {
	            if(load.traj) traj <- colSums(e$totpm[age.idx,,,drop=FALSE])
	            hch <- colSums(e$totpm.hch[age.idx,,,drop=FALSE])
	            if (length(age.idx) == max.age) quant <- .get.pop.quantiles(pop.pred, what='M', adjust=adjust, allow.negative.adj = allow.negative.adj)
	            else {if (length(age.idx) == 1) quant <- .get.pop.quantiles(pop.pred, what='Mage', adjust=adjust, allow.negative.adj = allow.negative.adj)[,age.idx,,]}
	        } else { # female
	            if(load.traj) traj <- colSums(e$totpf[age.idx,,,drop=FALSE])
	            hch <- colSums(e$totpf.hch[age.idx,,,drop=FALSE])
	            if (length(age.idx) == max.age) quant <- .get.pop.quantiles(pop.pred, what='F', adjust=adjust, allow.negative.adj = allow.negative.adj)
	            else {if (length(age.idx) == 1) quant <- .get.pop.quantiles(pop.pred, what='Fage', adjust=adjust, allow.negative.adj = allow.negative.adj)[,age.idx,,]}
	        }
	    }
	}
	if(load.traj) {
		if(typical.trajectory) {
			traj.idx <- bayesTFR:::get.typical.trajectory.index(traj)
		} else {
			thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(traj)[2])
			if (thintraj$nr.points > 0) 
		 		traj.idx <- thintraj$index
		}
	}
	if(!is.null(traj)) 
	 	rownames(traj) <- litem('proj.years.pop', pop.pred, pop.pred$proj.years+2)
	return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch))
}


get.pop.trajectories.multiple.age <- function(pop.pred, country, sex=c('both', 'male', 'female'), 
												age='all', nr.traj=NULL, proportion=FALSE, typical.trajectory=FALSE,
												adjust=FALSE, ...) {
	# Like get.pop.trajectories() but it doesn't sum up over ages.
	# Called when creating pop pyramid and pop.byage.*. Doesn't handle potential support ratio.
	age.idx <- traj.idx <- traj <- quant <- hch <- NULL
	e <- new.env()
	if (.load.traj.file(pop.output.directory(pop.pred), country, e)) {
		sex <- match.arg(sex)
		max.age <- dim(e$totpm)[1] # should be 27
		age.idx <- if(age[1]=='all') 1:max.age else age
		if(sex == 'both') {
			traj <- e$totpm[age.idx,,,drop=FALSE] + e$totpf[age.idx,,,drop=FALSE]
			hch <- e$totpm.hch[age.idx,,,drop=FALSE] + e$totpf.hch[age.idx,,,drop=FALSE]
		} else {
			if(sex=='male') {
				traj <- e$totpm[age.idx,,,drop=FALSE] 
				quant <- .get.pop.quantiles(pop.pred, what='Mage', adjust=adjust, ...)[,age.idx,,]
				hch <- e$totpm.hch[age.idx,,,drop=FALSE]
			} else {
				traj <- e$totpf[age.idx,,,drop=FALSE]
				quant <- .get.pop.quantiles(pop.pred, what='Fage', adjust=adjust, ...)[,age.idx,,]
				hch <- e$totpf.hch[age.idx,,,drop=FALSE]
			}
			if(proportion) {
				totpop <- (apply(e$totpm[,,,drop=FALSE], c(2,3), sum) + apply(e$totpf[,,,drop=FALSE], c(2,3), sum))
				for(iage in 1:dim(traj)[1])
					traj[iage,,] <- traj[iage,,]/totpop
			}
		}
		if(typical.trajectory) 
			traj.idx <- bayesTFR:::get.typical.trajectory.index(traj)
		else {
			thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(traj)[3])
			if (thintraj$nr.points > 0)
				traj.idx <- thintraj$index
		}
	} 
	if(!is.null(traj)) 
	 	dimnames(traj)[[2]] <- litem('proj.years.pop', pop.pred, pop.pred$proj.years+2) # pop.pred$proj.years
	return(list(trajectories=traj, index=traj.idx, age.idx=age.idx, quantiles=quant, half.child=hch))
}

is.saved.pi <- function(pop.pred, pi, warning=TRUE) {
	if(length(pi) == 0) return(NULL)
	is.valid.pi <- rep(NA, length(pi))
	quantile.values <- as.numeric(dimnames(pop.pred$quantiles)[[2]])
	for (i in 1:length(pi)) {
		al <- (1+pi[i]/100)/2		
		is.valid.pi[i] <- any(round(quantile.values,6)==round(al,6))
		if(!is.valid.pi[i] && warning)
			warning(pi[i], '% interval not available.')
	}
	return(is.valid.pi)
}

get.pop.traj.quantiles.byage <- function(quantile.array, pop.pred, country.index, country.code, year.index,
									trajectories=NULL, pi=80, q=NULL, ...) {
	# quantile.array should be 4d-array (country x age x quantiles x time)
	al <- if(!is.null(q)) q else c((1-pi/100)/2, (1+pi/100)/2)
	found <- FALSE
	if(!is.null(quantile.array)) {
		quantile.values <- as.numeric(dimnames(quantile.array)[[3]])
		alidx<-round(quantile.values,6)==round(al[1],6)
		cqp <- NULL
		if (any(alidx)) { # pre-saved quantiles
			cqp <- quantile.array[country.index,,alidx,year.index]
			if(length(al) > 1) {
				alidx2 <- round(quantile.values,6)==round(al[2],6)
				cqp <- rbind(cqp, quantile.array[country.index, ,alidx2,year.index])
			}
			found <- TRUE
		} 
	}
	if(!found) { # non-standard quantiles
		if(is.null(trajectories)) {
			warning('Quantiles not found')
			return(NULL)	
		}
		cqp <- apply(trajectories[,year.index,,drop=FALSE], 1, 
						quantile, al, na.rm = TRUE)
	}
	return(cqp)
}


get.pop.traj.quantiles <- function(quantile.array, pop.pred, country.index=NULL, country.code=NULL, 
									trajectories=NULL, pi=80, q=NULL, reload=TRUE, ...) {
	# quantile.array should be 3d-array (country x quantiles x time). 
	# If country.index is NULL or there is just one country in the prediciton object, 
    # the country dimension can be omitted 
	al <- if(!is.null(q)) q else c((1-pi/100)/2, (1+pi/100)/2)
	found <- FALSE
	if(!is.null(quantile.array)) {
		if((is.null(country.index) || nrow(pop.pred$countries) == 1) && length(dim(quantile.array))<3) {
			quantile.array <- abind(quantile.array, along=0)
			country.index <- 1
		}
		quantile.values <- as.numeric(dimnames(quantile.array)[[2]])
		alidx<-round(quantile.values,6)==round(al[1],6)
		cqp <- NULL
		if (any(alidx)) { # pre-saved quantiles
			cqp <- quantile.array[country.index, alidx,]
			if(length(al) > 1) {
				alidx2 <- round(quantile.values,6)==round(al[2],6)
				cqp <- rbind(cqp, quantile.array[country.index, alidx2,])
			}
			found <- TRUE
		} 
	}
	if(!found) { # non-standard quantiles
		if(is.null(trajectories) && !reload) {
			warning('Quantiles not found')
			return(NULL)	
		}
		do.reload <- FALSE
		if (is.null(trajectories)) {
			if(pop.pred$nr.traj > 0) do.reload <- TRUE
		} else { 
			if (dim(trajectories)[2] < 2000 && pop.pred$nr.traj > dim(trajectories)[2] && reload) do.reload <- TRUE
		}
		if(do.reload) {
			#load 2000 trajectories maximum for computing quantiles
			traj.reload <- get.pop.trajectories(pop.pred, country.code, nr.traj=2000, ...)
			trajectories <- traj.reload$trajectories
		}
		if (!is.null(trajectories)) {
			if(is.null(dim(trajectories))) trajectories <- abind(trajectories, along=2) # only one trajectory
			cqp <- apply(trajectories, 1, 
						quantile, al, na.rm = TRUE)
		} else {
			warning('Quantiles not found')
			return(NULL)
		}
	}
	return(cqp)
}

get.migration <- function(pop.pred, country, sex, is.observed=FALSE, VEenv=NULL) {
    par1 <- paste0('mig',tolower(sex))
    par2 <- paste0(tolower(sex), 'mig')
    res <- NULL	
    if(!is.null(VEenv)) {
        for(par in c(par1, par2)) {
            if(!is.null(VEenv[[par]])) {
                res <- VEenv[[par]]
                break
            }
        }
    }
    if(is.null(res)) {
        par3 <- paste0('MIG',tolower(sex))
        inputs <- if(is.observed) pop.pred$inputs$observed else pop.pred$inputs
        res <- .get.par.from.inputs(par3, inputs, country)
        if(!is.observed) {
            #add present year 
            obs <- .get.par.from.inputs(par3, pop.pred$inputs$observed, country)
            res <- cbind(obs[,ncol(obs)], res)
            colnames(res) <- pop.pred$proj.years
        }
        res <- abind(res, along=3)
    }
	return(res)
}

mid.period3d <- function(dat)
    (dat[,-1, ,drop = FALSE] + dat[,-dim(dat)[2],, drop = FALSE])/2.

mx.aggregate <- function(mx, pop, abridged = TRUE) {
    # Aggregate mx over sexes
    # mx and pop are lists with elements for male and female
    abr.deaths <- abr.pop <- list()
    for(s in c("male", "female")) {
        if(length(dim(pop[[s]]))<3) pop[[s]] <- abind(pop[[s]], along=3)
        # abridged average population split to 0-1 and 1-4
        abr.pop[[s]] <- mx[[s]]
        abr.pop[[s]][] <- NA
        apop <- mid.period3d(pop[[s]])
        if(abridged) 
            apop <- split.pop05(apop)
        if(dim(apop)[2] > dim(abr.pop[[s]])[2]) # remove time periods from apop to align with mx
            apop <- apop[,-(1:(dim(apop)[2] - dim(abr.pop[[s]])[2])),,drop = FALSE]
        itime <- (dim(abr.pop[[s]])[2] - dim(apop)[2] + 1):dim(abr.pop[[s]])[2]
        abr.pop[[s]][,itime,] <- pmax(apop, 1e-4)
        # abridged deaths
        abr.deaths[[s]] <- abr.pop[[s]] * mx[[s]]
    }
    # combine to mx
    denom <- abr.pop$male + abr.pop$female
    aggr.mx <- (abr.deaths$male + abr.deaths$female)/denom
    #aggr.mx[is.na(aggr.mx) & denom == 0] <- 1 # it is NA where pop is zero; set the mx to 1
    return(aggr.mx)
}

get.mx <- function(mxm, sex, age05=c(FALSE, FALSE, TRUE), abridged = TRUE, pop = NULL) {
    if(sex == "T")  # mx for both sexes
        mxm <- mx.aggregate(mxm, pop, abridged = abridged)
    if(length(dim(mxm))<3) mxm <- abind(mxm, along=3)
    if(age05[3]) {
        res1 <- LifeTableMxCol(mxm[,, 1], colname='mx', sex=sex, age05=age05, abridged = abridged)
        if(is.null(dim(res1))) res1 <- abind(res1, along=2)
        res <- array(0, dim=c(dim(res1)[1], dim(res1)[2], dim(mxm)[3]))
        res[,,1] <- res1
        if(dim(mxm)[3]> 1) { 
            for (itraj in 2:dim(mxm)[3]) {
                res[,, itraj] <- LifeTableMxCol(mxm[,, itraj], colname='mx', sex=sex, age05=age05, abridged = abridged)
            }
        }
        dimnames(res)[[2]] <- dimnames(mxm)[[2]]
        return(res)
    }
    if(abridged) {
        if(!age05[2]) mxm <- mxm[-2,,,drop=FALSE]
        if(!age05[1]) mxm <- mxm[-1,,,drop=FALSE]
    }
	return (mxm)
}

.mx.replace.na.for.old.ages <- function(mx) {
    # replace NA at old ages with 1
    mask <- is.na(mx)
    mask[1:17,,] <- mask[1:17,,] & FALSE
    mx[mask] <- 1
    return(mx)
}

.get.lt.col <- function(ltcol, mxm, sex, age05=c(FALSE, FALSE, TRUE), 
                        abridged = TRUE, replace.na = TRUE, pop = NULL) {
    if(sex == "T") mxm <- mx.aggregate(mxm, pop, abridged = abridged) # aggregate over sexes
    if(length(dim(mxm))<3) mxm <- abind(mxm, along=3)
    if(abridged && replace.na && any(is.na(mxm)) && any(!is.na(mxm)))
        mxm <- .mx.replace.na.for.old.ages(mxm)
    val1 <- LifeTableMxCol(mxm[,, 1], colname=ltcol, sex=sex, age05=age05, abridged = abridged)
    if(is.null(dim(val1))) val1 <- abind(val1, along=2)
    val <- array(0, dim=c(dim(val1)[1], dim(val1)[2], dim(mxm)[3]))
    val[,,1] <- val1
    dimnames(val)[[2]] <- dimnames(mxm)[[2]]
    if(dim(mxm)[3] <= 1) return(val)
    for (itraj in 2:dim(mxm)[3]) {
        val[,, itraj] <- LifeTableMxCol(mxm[,, itraj], colname=ltcol, sex=sex, age05=age05, abridged = abridged)
    }
    return (val)
}

get.qx <- function(...) {
    return(.get.lt.col('qx', ...))
}

get.ex <- function(...) {
    return(.get.lt.col('ex', ...))
}

get.ax <- function(...) {
    return(.get.lt.col('ax', ...))
}

get.survival <- function(mxm, sex, age05=c(FALSE, FALSE, TRUE), abridged = TRUE,
                         replace.na = TRUE, pop = NULL) {
    if(sex == "T") mxm <- mx.aggregate(mxm, pop, abridged = abridged) # aggregate over sexes
	if(length(dim(mxm))<3) mxm <- abind(mxm, along=3)
	# sx21 <- dim(mxm)[1] < 27
	# for (itraj in 1:dim(mxm)[3]) {
		# LLm <- LifeTableMxCol(mxm[,, itraj, drop=FALSE], colname='Lx', sex=sex, age05=age05)
		# if(itraj == 1) {
			# if(is.null(dim(LLm))) LLm <- abind(LLm, along=2)
			# sx <- array(0, dim=c(dim(LLm)[1], dim(LLm)[2], dim(mxm)[3]), dimnames=vector("list", 3))
			# sr <- rep(0, dim(LLm)[1])
		# }
		# sx[,, itraj] <- apply(LLm, 2, function(x, sr) { if(!any(is.na(x))) {
															# res.sr <- if(sx21) .C("get_sx21_21", as.numeric(x), sx=sr)
																		# else .C("get_sx27", as.numeric(x), sx=sr)
															# return(res.sr$sx)
														# } 
														# return(rep(NA, length(x)))
													# }, sr)
	# }
	if(replace.na && any(is.na(mxm)) && any(!is.na(mxm)))
	    mxm <- .mx.replace.na.for.old.ages(mxm)
	sx1 <- LifeTableMxCol(mxm[,, 1], colname='sx', sex=sex, age05=age05, abridged = abridged)
	if(is.null(dim(sx1))) sx1 <- abind(sx1, along=2)
	sx <- array(0, dim=c(dim(sx1)[1], dim(sx1)[2], dim(mxm)[3]),
	            dimnames = c(dimnames(sx1)[1], dimnames(mxm)[2:3]))
	sx[,,1] <- sx1
	if(dim(mxm)[3] <= 1) return(sx) # one trajectory
	for (itraj in 2:dim(mxm)[3]) 
		sx[,, itraj] <- LifeTableMxCol(mxm[,, itraj], colname='sx', sex=sex, age05=age05, abridged = abridged)
	return (sx)
}

get.popVE.trajectories.and.quantiles <- function(pop.pred, country, 
									event=c('births', 'deaths', 'survival', 'fertility', 'qx', 'mx', 'migration', 'pasfr', 'ex', 'ax'), 
									sex=c('both', 'male', 'female'), age='all', sum.over.ages=TRUE,
 									nr.traj=NULL, q=NULL, typical.trajectory=FALSE, is.observed=FALSE,
									allow.higher.ages = FALSE, ...) {
 	# get trajectories and quantiles for vital events and other indicators
 	input.indicators <- c('migration')
 	#input.indicators <- c()
 	life.table.indicators <- c('survival', 'qx', 'mx', 'ex', 'ax')
 	quant <- hch <- age.idx <- traj <- traj.idx <-  NULL
 	event <- match.arg(event)
 	sex <- match.arg(sex)
 	time.labels <- colnames(pop.pred$inputs$pop.matrix$male)
 	
 	#if (!is.element(event, input.indicators)) {
		traj.file <- file.path(pop.output.directory(pop.pred), paste('vital_events_country', country, '.rda', sep=''))
		if (!file.exists(traj.file) && !is.element(event, input.indicators)) 
			return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch, event=event))
		myenv <- new.env()
		if (file.exists(traj.file)) load(traj.file, envir=myenv)
		if(is.observed) myenv <- myenv$observed
	#}
	max.age.index.allowed <- age.length.all(pop.pred$annual)
	min.age.index.allowed <- 1
	if(is.observed) {
		nperiods <- length(get.pop.observed.periods(pop.pred))
		if(!allow.higher.ages) max.age.index.allowed <- age.length.all(pop.pred$annual, observed = TRUE)
	}
	age.normal <- TRUE
	subtract.max.age <- 0
	if(is.element(event, life.table.indicators)) {
	    if(!pop.pred$annual) min.age.index.allowed <- -1
		mx <- list(male=myenv$mxm, female=myenv$mxf, male.hch=myenv$mxm.hch, female.hch=myenv$mxf.hch)
		mx$both <- list(male = mx$male, female = mx$female)
		mx$both.hch <- list(male = mx$male.hch, female = mx$female.hch)
		for(s in c("male", "female", "both")) {
		    if(sex == s) {
		        sex.names <- s
		        if(!is.observed) sex.names <- c(sex.names, paste0(s, ".hch"))
		        break
		    }
		}
		poplist <- NULL
		if(sex=='both') {
		    if(!is.observed) {
		        pop <- list(male = get.pop.trajectories.multiple.age(pop.pred, country, sex = "male"),
		                    female = get.pop.trajectories.multiple.age(pop.pred, country, sex = "female"))
		        poplist <- list(both = list(male = pop$male$trajectories, female = pop$female$trajectories),
		                        both.hch = list(male = pop$male$half.child, female = pop$female$half.child))
		    } else {
		        poplist <- list(both = list(male = get.pop.observed(pop.pred, country, sex = "male", sum.over.ages=FALSE), 
		                                    female = get.pop.observed(pop.pred, country, sex = "female", sum.over.ages=FALSE)
		                                    )
		                    )
		    }
		}
		sexarg <- list(male = "M", female = "F", male.hch = "M", female.hch = "F", both = "T", both.hch =  "T")
		alltraj <- list(male=NULL, female=NULL, both = NULL, male.hch=NULL, female.hch=NULL, both.hch = NULL)
		age05 <- c(FALSE, FALSE, TRUE)
		if(! pop.pred$annual) {
		    if(age[1]!='all' && (any(age < 1))) {
			    age05 <- rep(TRUE, 3)
			    age.normal <- FALSE
			    subtract.max.age <- 2
		    }
		    abridged <- TRUE
	    } else abridged <- FALSE
		for(sn in sex.names) {
			if(!is.null(mx[[sn]])) {
				if(is.observed) {
				    if(sn == "both") {
				        dimnames(mx[[sn]]$female)[[2]] <- dimnames(mx[[sn]]$male)[[2]] <- time.labels[(length(time.labels)-dim(mx[[sn]]$male)[2]+1):length(time.labels)]
				    } else dimnames(mx[[sn]])[[2]] <- time.labels[(length(time.labels)-dim(mx[[sn]])[2]+1):length(time.labels)]
				}
				alltraj[[sn]] <- do.call(paste0('get.', event), 
									list(mx[[sn]], sex=sexarg[[sn]], age05=age05, pop = poplist[[sn]], abridged = abridged))
			}
		}
	} else { # no life table events
		if (is.element(event, input.indicators)) { # migration
 			alltraj <- list(male=NULL, female=NULL, male.hch=NULL, female.hch=NULL)
 			sex.index <- 1:2
			if(sex=='male') sex.index <- sex.index[1]
			if(sex=='female') sex.index <- sex.index[2]
			for(is in sex.index) {
				alltraj[[names(alltraj)[is]]] <- do.call(paste0('get.', event), 
									list(pop.pred, country, sex=c("M","F","M","F")[is], is.observed=is.observed, VEenv=myenv))
			}
		} else
		    alltraj <- switch(event,
		                      births = list(male=myenv$btm, female=myenv$btf, male.hch=myenv$btm.hch, female.hch=myenv$btf.hch),
		                      deaths = list(male=myenv$deathsm, female=myenv$deathsf, male.hch=myenv$deathsm.hch, female.hch=myenv$deathsf.hch),
		                      fertility = list(female=myenv$asfert, female.hch=myenv$asfert.hch),
		                      pasfr = list(female=myenv$pasfert, female.hch=myenv$pasfert.hch)
		                )
	}
	has.hch <- !is.observed && (!is.null(alltraj$male.hch) || !is.null(alltraj$female.hch) || !is.null(alltraj$both.hch))
	max.age <- NULL
	for(s in c("male", "female", "both")) {# max.age should be 7, 21, 27, 28 (28 only if zero is explicitely included in 'age', so never when age=='all')
        if(is.null(alltraj[[s]])) next
        max.age <- dim(alltraj[[s]])[1] - subtract.max.age
        trajdimnames <- dimnames(alltraj[[s]])
        break
	}
	if(is.null(max.age)) # no trajectories available
	    return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch, event=event))
	max.age.shift <- 0 # shift for fertility indicators
	age.obs.length <- age.length.all(pop.pred$annual, observed = TRUE)
	first.fert.idx <- age.index.fert(pop.pred$annual)[1]
	if(max.age < age.obs.length) max.age.shift <- first.fert.idx - 1
	
	age.idx <- age.idx.raw  <- if(age[1]=='all') 1:max.age else age[age <= (max.age + max.age.shift)] 
	quantiles <- if(is.null(q)) get.quantiles.to.keep() else q

	if(event %in% c('births', 'fertility', 'pasfr')) {
		if(age[1] != 'all') {
			age.idx <- age.idx - max.age.shift # translate age index into mother's child-bearing age index
			if(length(age.idx)==0 || max(age.idx) > max.age || min(age.idx) < 1){
			    if(pop.pred$annual) stop('Age index for ', event, ' must be between ', max.age.shift, ' and ', max.age + max.age.shift - 1, '.')
			    allowed.ages <- ages.fert(pop.pred$annual)
			    stop('Age index for ', event, ' must be between ', max.age.shift + 1, ' (age ', allowed.ages[1], '-', allowed.ages[1]+4, ') and ', 
				     max.age + max.age.shift, ' (age ', allowed.ages[length(allowed.ages)], '-', allowed.ages[length(allowed.ages)]+4, ').')
			}
		} else age.idx.raw <- age.idx + max.age.shift
	} 
	if(length(age.idx)==0 || max(age.idx) > max.age.index.allowed || (min(age.idx) < min.age.index.allowed)) {
	    age.ranges <- c(min.age.index.allowed, max.age.index.allowed, max.age)
	    if(pop.pred$annual) age.ranges <- age.ranges - 1 # age argument is assumed to be already shifted for 1-year age groups, i.e. starting with 1 instead of allowed 0
		stop('Age index must be between ', age.ranges[1], ' (first age category) and ', min(age.ranges[3], age.ranges[2]),  
						' (open-ended age category).')
	}
	if(!age.normal) {
		age.idx.raw <- age.idx
		age.idx <- age.idx+2
	}
	#if(!is.observed) stop('')
	hch <- NULL
	if(event  %in% c('fertility', 'pasfr')) sex <- 'female'
	if(sex == 'both' && !is.element(event, life.table.indicators)) { # summing over sexes
		if(sum.over.ages) {
			traj <- colSums(alltraj$male[age.idx,,,drop=FALSE]) + colSums(alltraj$female[age.idx,,,drop=FALSE])
			if(has.hch) hch <- colSums(alltraj$male.hch[age.idx,,,drop=FALSE]) + colSums(alltraj$female.hch[age.idx,,,drop=FALSE])
		} else {
			traj <- alltraj$male[age.idx,,,drop=FALSE] + alltraj$female[age.idx,,,drop=FALSE]
			if(has.hch) hch <- alltraj$male.hch[age.idx,,,drop=FALSE] + alltraj$female.hch[age.idx,,,drop=FALSE]
		}
	} else { # one sex
		if(sum.over.ages) {
			traj <- colSums(alltraj[[sex]][age.idx,,,drop=FALSE])
			if(has.hch) hch <- colSums(alltraj[[paste(sex,'hch', sep='.')]][age.idx,,,drop=FALSE])
		} else {
			traj <- alltraj[[sex]][age.idx,,,drop=FALSE]
			if(has.hch) hch <- alltraj[[paste(sex,'hch', sep='.')]][age.idx,,,drop=FALSE]
		}
	}
	if(is.observed) {
		if(length(dim(traj)) < 3) # age dimension is missing
			traj <- abind(traj, NULL, along=0)
		 if(dim(traj)[[2]] < nperiods) {		
			traj <- abind(array(NA, dim=c(dim(traj)[[1]], nperiods-dim(traj)[[2]], dim(traj)[[3]]), 
						dimnames=list(NULL, colnames(pop.pred$inputs$pop.matrix$male)[1:(nperiods-dim(traj)[[2]])], NULL)),
						traj, along=2)
		}
	}
	quant <- NULL
	if(!is.observed) {
		if(sum.over.ages) { # quantiles are 2-d arrays
			quant <- apply(traj, 1, quantile, quantiles, na.rm = TRUE)
			dimnames(quant) <- list(quantiles, trajdimnames[[2]])
			traj.for.thinning <- traj
			year.dim <- 1
		} else { # quantiles are 3-d arrays: age x quantiles x period
			quant <- aperm(apply(traj, c(1,2), quantile, quantiles, na.rm = TRUE), c(2,1,3))
			dimnames(quant) <- list(#pop.pred$ages[age.idx.raw], 
								trajdimnames[[1]][age.idx], quantiles, trajdimnames[[2]])
			traj.for.thinning <- traj[1,,]
			if(is.vector(traj.for.thinning)) # in case of one trajectory
				traj.for.thinning <- abind(traj.for.thinning, NULL, along=2)
			year.dim <- 2
		}
	}
	if((is.null(nr.traj) || nr.traj > 0) && !is.observed) {
		if(typical.trajectory) {
			traj.idx <- bayesTFR:::get.typical.trajectory.index(traj.for.thinning)
			if(!sum.over.ages && length(age.idx) > 1) {
				for(i in 2:length(age.idx))
					traj.idx <- rbind(traj.idx, bayesTFR:::get.typical.trajectory.index(traj[i,,]))
			}
		} else {
			thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(traj.for.thinning)[2])
			if (thintraj$nr.points > 0) 
		 		traj.idx <- thintraj$index
		}
		dimnames(traj)[[year.dim]] <- trajdimnames[[2]]
	} else if(!is.observed) traj <- NULL
 	
	return(list(trajectories=traj, index=traj.idx, quantiles=quant, 
				age.idx=age.idx, age.idx.raw=age.idx.raw, half.child=hch, event=event))
}


get.age.labels <- function(ages, collapsed=FALSE, age.is.index=FALSE, last.open=FALSE, single.year = FALSE) {
	all.age <- ages.all(annual = single.year, observed = FALSE)
	ages.idx <- if(age.is.index) ages else which(is.element(all.age, ages))
	if(is.element(0, ages.idx) && !single.year) { # age 1-4 is included
		all.age <- c(-1, all.age)
		idx14 <- which(is.element(ages.idx,0))
		ages.idx[idx14] <- 1
		idx01 <- which(is.element(ages.idx,-1))
		ages.idx[-c(idx14,idx01)] <- ages.idx[-c(idx14,idx01)] + 1	
	}
	if(is.element(-1, ages.idx) && !single.year) { # age 0-1 is included
		all.age <- c(-2, all.age)
		idx01 <- which(is.element(ages.idx,-1))
		ages.idx[idx01] <- 1
		ages.idx[-idx01] <- ages.idx[-idx01] + 1
	}
	ages.idx.shift <- ages.idx+1
	if(collapsed) {
		ages.idx.dif <- which(!is.element(ages.idx, ages.idx.shift))
		ages.idx.shift <- ages.idx.shift[!is.element(ages.idx.shift, ages.idx)]
		ages.idx <- ages.idx[ages.idx.dif]
	}
	lages <- all.age[ages.idx]
	uages <- all.age[ages.idx.shift]
	l <- length(lages)
	from <- all.age[ages.idx[1:(l-1)]]
	to <- all.age[ages.idx.shift[1:(l-1)]]-1
	if(!single.year) {
	    which.01 <- which(from < -1)
	    which.14 <- which(from < 0 & from > -2)
	    if(length(which.01)>0) { # includes age 0-1
		    from[which.01] <- 0
		    to[which.01][!is.na(to[which.01])] <- 1
	    }
	    if(length(which.14)>0) { # includes age 1-4
		    from[which.14] <- 1
		    to[which.14][!is.na(to[which.14])] <- 4
	    }
	}
	result <- if(l==1 && is.na(to)) paste0(from, '+') # open-ended
			  else ifelse(from == to, from, paste0(from, '-', to))
	if (l > 1) result <- c(result, if(is.na(all.age[ages.idx.shift[l]]) || last.open) paste0(all.age[ages.idx[l]], '+')
			else ifelse(all.age[ages.idx[l]] == all.age[ages.idx.shift[l]]-1, all.age[ages.idx[l]], 
			            paste0(all.age[ages.idx[l]], '-', all.age[ages.idx.shift[l]]-1)))
	return(as.character(result))
}	

.get.year.index <- function(year, years, annual = FALSE) {
    if(annual) return(if(year %in% years) which(years == year) else NA)
	lyears <- length(years)
	res <- as.integer(cut(year, labels=1:lyears, breaks=c(years-3, years[lyears]+2)))
	return(res)
	
	#breaks <- c(years-3, years[lyears]+2)
	#h <- try(hist(year, breaks=breaks, plot=FALSE)$count, silent=TRUE)
	#return(if(inherits(h, "try-error")) NULL else which(h > 0)[1])
}
get.pop.prediction.periods <- function(pop.pred, end.time.only=FALSE) {
    if(pop.pred$annual) return(pop.pred$proj.years)
	if(end.time.only) return(litem('proj.years.pop', pop.pred, pop.pred$proj.years+2))
	return(sapply(lapply(pop.pred$proj.years, '+', c(-3, 2)), paste, collapse='-'))
}
get.prediction.year.index <- function(pop.pred, year) {
	years <- pop.pred$proj.years
	return(.get.year.index(year, years, annual = pop.pred$annual))
}

get.observed.year.index <- function(pop.pred, year) 
	return(.get.year.index(year, as.integer(colnames(pop.pred$inputs$pop.matrix$male)), annual = pop.pred$annual))

get.pop.observed.periods <- function(pop.pred, end.time.only=FALSE) {
	years <- as.integer(colnames(pop.pred$inputs$pop.matrix$male))
	if(pop.pred$annual) return(years)
	if(is.null(pop.pred$proj.years.pop)) { # assuring compatibility before version 5.0-1 where years are the middle time points
		if(end.time.only) return(years + 2)
	} else {
		if(end.time.only) return(years)
		years <- years - 2
	}
	return(sapply(lapply(years, '+', c(-3, 2)), paste, collapse='-'))
}

get.predORobs.year.index <- function (pred, year) 
{
    projection.index <- get.prediction.year.index(pred, year)
    projection <- TRUE
    if (is.null(projection.index) || is.na(projection.index)) {
        projection <- FALSE
        projection.index <- get.observed.year.index(pred, year)
    }
    return(c(index = projection.index, is.projection = projection))
}

get.quantiles.to.keep <- function() {
	#return(c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1))
	return(c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975))
}

get.countries.table.bayesPop.prediction <- function(object, ...) 
	return(object$countries)
	
get.pop <- function(object, pop.pred, aggregation=NULL, observed=FALSE, ...) {
	split.object <- strsplit(object, '_', fixed=TRUE)[[1]]
	what <- substr(split.object[1],1,1)
	# Is it a vital event
	has.ve <- is.element(what, names(get.expression.indicators()))
	split.object[1] <- sub(what, '', split.object[1], fixed=TRUE)
	
	# Parse country
	country.string <- regmatches(split.object[1], 
						regexpr('^[[:alnum:]]*\\[|^[[:alnum:]]*\\{?|^XXX\\[|^XXX\\{?', split.object[1])) 
	country.code.char <- gsub('\\[|\\{', '', country.string)
	if(nchar(country.code.char) == 0) stop('No country specified.')
	if(country.code.char != 'XXX') {
		e <- new.env()
		data('iso3166', package='bayesTFR', envir=e)
		charcodename <- if(nchar(country.code.char)==2) 'charcode' else 'charcode3'
		country.code.iso <- e$iso3166[e$iso3166[,charcodename] == country.code.char,'uncode']
		if(length(country.code.iso)>0)
			country.code <- country.code.iso
		else {
			country.code <- as.integer(country.code.char)
			if(is.na(country.code)) stop('Unknown country ', country.code.char)
		}
	} else {
		if(!observed) stop('Country must be specified. No XXX allowed.')
		country.code <- country.code.char
	}
	# remove country code from the string
	split.object[1] <- sub(paste(country.code.char, sep=''), '', split.object[1], fixed=TRUE)
	if(nchar(split.object[1])<=0) split.object <- split.object[-1]
	
	# Parse sex
	sex <- 'both'
	if(length(split.object) > 0) {
		for(sx in c('F', 'M')) {
			sx.idx <- grep(sx, split.object)
			if(length(sx.idx) > 0) {
				sex <- list(F='female', M='male')[[sx]]
				split.object[sx.idx] <- sub(sx, '', split.object[sx.idx])
				break
			}
		}
	}
	# Parse age
	sum.over.ages <- TRUE
	age <- 'all'
	age.part.idx <- grep("\\[|\\{", split.object)
	if(length(age.part.idx) > 0) {
		if(length(age.part.idx) > 1) stop('Only one age vector is allowed.')
		age <- eval(parse(text=gsub('\\[|\\]|\\{|\\}', '', split.object[age.part.idx])))
		if(is.null(age)) age <- 'all' else if(pop.pred$annual) age <- age + 1
		if(grepl("{", split.object, fixed=TRUE)) sum.over.ages <- FALSE
	}
	# find country (search aggregations if not found, or if it is an aggregation search base prediction object)
	if(country.code != 'XXX') {
		country.object <- get.country.object(country.code, country.table=pop.pred$countries)
		if(is.null(country.object$code)) {
			if(pop.pred$is.aggregation) {
				base.pred <- get.pop.prediction(sim.dir=pop.pred$base.directory)
				country.object <- get.country.object(country.code, country.table=base.pred$countries)
				if(is.null(country.object$code)) stop('Invalid country code used.')
				pop.pred <- base.pred
			} else {
				av.aggrs <- available.pop.aggregations(pop.pred)
				indep.idx <- which(is.element('country', av.aggrs))
				if(length(indep.idx) > 0)  # put country-type aggregation first
					av.aggrs <- c('country', av.aggrs[-indep.idx])
				if(is.null(aggregation)) aggregation <- av.aggrs
				for(aggr in aggregation) {
					if(!is.element(aggr, av.aggrs)) {warning('Aggregation', aggr, 'not available.'); next}
					aggr.obj <- get.pop.aggregation(pop.pred=pop.pred, name=aggr)
					country.object <- get.country.object(country.code, country.table=aggr.obj$countries)
					if(!is.null(country.object$code)) {pop.pred <- aggr.obj; break}
				}
				if(is.null(country.object$code)) stop('Invalid country code used.')
			}
		}
	}
	if(observed) {
		if(country.code != 'XXX') {
			if(!has.ve) {
				traj <- get.pop.observed.with.age(pop.pred, country=country.object$code, sex=sex, age=age)
				d <- traj$data[traj$age.idx,,drop = FALSE]
			} else {
				traj <- get.popVE.trajectories.and.quantiles(pop.pred, country.object$code, 
											event=get.expression.indicators()[[what]], sex=sex, age=age, 
											sum.over.ages=FALSE, is.observed=TRUE, ...)
				traj$age.idx <- traj$age.idx.raw
				d <- traj$trajectories
			}
		    if(is.null(d)) return(NULL)
			if(sum.over.ages) {
				if(!is.null(nrow(d)) && nrow(d) == 0) # country not found in the observed data
					d <- rep(NA, ncol(d))
				else d <- colSums(d)
				data <- as.matrix(d) # adds trajectory dimension if missing
				dim(data) <- c(1, dim(data)) # adding age dimension
				dimnames(data) <- list(NULL, colnames(traj$data), NULL)
			} else {# only if it was not summed up, because then the as.matrix command adds a dimension
				data <- if(is.null(dim(d)) || !is.array(d)) as.matrix(d) else d
				#data <- as.matrix(d)
				#if(age[1]=='all') { # extend to 27 age categories
				#	if(dim(data)[1] < 27) data <- rbind(data, matrix(0, nrow=27-nrow(data), ncol=ncol(data)))
				#	traj$age.idx <- c(traj$age.idx, (nrow(d)+1):27)
				#}
				if(length(dim(data)) < 3) {
					dim(data) <- c(dim(data), 1)
					dimnames(data)[2:length(dim(d))] <- dimnames(d)[2:length(dim(d))]
				}
			}
			dimnam <- dimnames(data)
			dim(data) <- c(1, dim(data)) # adding country dimension
			dimnames(data)[2:length(dim(data))] <- dimnam
		} else { # multiple countries
		    if(has.ve) stop("XXX is only allowed for population.")
			traj <- get.pop.observed.multiple.countries(pop.pred, countries=pop.pred$countries$code, sex=sex, 
														age=age, sum.over.ages=sum.over.ages)
			data <- traj$data
			if(!sum.over.ages) {
				if(age[1]=='all') { # extend to 27 age categories
				    lages <- age.length.all(pop.pred$annual)
					traj$age.idx <- c(traj$age.idx, (dim(data)[[2]]+1):lages)
					data <- abind(data, array(0, c(dim(data)[[1]], lages-dim(data)[[2]], dim(data)[[3]])), along=2)
				}
				dim(data) <- c(dim(data), 1) # adding trajectory dimension
				dimnames(data)[[3]] <- dimnames(traj$data)[[3]]
			} else {
				dim(data) <- c(dim(data)[[1]], 1, dim(data)[[2]], 1) # adding age and trajectory dimension
				dimnames(data)[[3]] <- dimnames(traj$data)[[2]]
			}
			if(!is.null(dimnames(traj$data[[1]]))) dimnames(data)[[1]] <- dimnames(traj$data)[[1]]
		}
	} else { # projections
		if(what == 'P')
			traj <- .get.trajectories(sum.over.ages=sum.over.ages, pop.pred, country=country.object$code, sex=sex, age=age, ...)
		else {
			if(has.ve){
				traj <- get.popVE.trajectories.and.quantiles(pop.pred, country.object$code, 
											event=get.expression.indicators()[[what]], sex=sex, age=age, 
											sum.over.ages=sum.over.ages, ...)
				time.dim <- if(sum.over.ages) 1 else 2
				if(dim(traj$trajectories)[time.dim] < length(pop.pred$proj.years)) { # add current year
					traj$trajectories <- abind(array(NA, dim=dim(traj$trajectories)[-time.dim]), traj$trajectories, along=time.dim)
					dimnames(traj$trajectories)[[time.dim]] <- pop.pred$proj.years
				}
				traj$age.idx <- traj$age.idx.raw
			} else stop('Indicator ', what, 'not implemented.')
		}
		data <- traj$trajectories
		dim(data) <- c(1,dim(data)) # adding country  dimension
		if(sum.over.ages) {
			dim(data) <- c(1,dim(data)) # adding age dimension
			dimnames(data) <- list(NULL, NULL, dimnames(traj$trajectories)[[1]], NULL)
		} else dimnames(data) <- list(NULL, dimnames(traj$trajectories)[[1]], dimnames(traj$trajectories)[[2]], NULL)
	}
	if(length(traj$age.idx) == dim(data)[[2]]) dimnames(data)[[2]] <- traj$age.idx
	    #dimnames(data)[[2]] <- if(pop.pred$annual) traj$age.idx - 1 else traj$age.idx
	    
	return(data)
}

.get.trajectories <- function(sum.over.ages=TRUE, ...){
	traj <- if(sum.over.ages) get.pop.trajectories(...) else get.pop.trajectories.multiple.age(...)
	return(traj)
}

.parse.pop.expression <- function(expression, args='...') {
	# Add spaces around binary operators so that expression components can be identified
	#expression <- gsub('(\\*|\\/|\\+|\\-|\\^|\\%\\%|\\%\\/\\%)', ' \\1 ', expression)
	expression <- gsub('(\\*|\\/|\\+|\\^|\\%\\%|\\%\\/\\%)', ' \\1 ', expression) # removed minus 
	indicators <- c('P', names(get.expression.indicators()))
	indOR <- paste(indicators, collapse='|')
	# Replace expression components by 'get.pop' calls
	return (gsub(paste('((', indOR, ')[[:graph:]]*[[:alnum:]]|(', indOR, ')[[:graph:]]*\\]|(', indOR, ')[[:graph:]]*\\})', sep=''),
			paste("get.pop('\\1', pop.pred,", args, ")"), expression))
}

get.pop.ex <- function(expression, pop.pred, observed = FALSE, as.dt = FALSE, ...) {
	# Return trajectories or observed values for an expression (defined by time)
    if(grepl("XXX", expression, fixed = TRUE))
        return(get.pop.ex.all(expression, pop.pred, observed = observed, as.dt = as.dt, ...))
	if(observed) 
		return(get.pop.observed.from.expression(expression, pop.pred, as.dt = as.dt, ...))
	result <- get.pop.trajectories.from.expression(expression, pop.pred, ...)
	res <- result$trajectories[,result$index]
	if(as.dt){
	    trajectory <- NULL
	    res <- data.table::data.table(result$trajectories[,result$index, drop = FALSE], keep.rownames = TRUE)
	    res <- data.table::melt(res, id.vars = "rn", value.name = "indicator", variable.name = "trajectory",
	                  variable.factor = FALSE)
	    res[, trajectory := as.integer(gsub("V", "", trajectory))]
	    data.table::setnames(res, "rn", "year")
	}
	return(res)
}

get.pop.exba <- function(expression, pop.pred, observed = FALSE, as.dt = FALSE, ...) {
	# Return trajectories or observed values for an expression (defined by age, i.e. includes {})
    if(grepl("XXX", expression, fixed = TRUE))
        return(get.pop.exba.all(expression, pop.pred, observed = observed, as.dt = as.dt, ...))
	if(observed) 
		return(get.pop.observed.from.expression.multiple.age(expression, pop.pred, as.dt = as.dt, ...))
	result <- get.pop.trajectories.from.expression.multiple.age(expression, pop.pred, ...)
	res <- result$trajectories[,,result$index]
	if(as.dt){
	    age <- NULL
	    res <- data.table::as.data.table(result$trajectories[,,result$index, drop = FALSE], keep.rownames = TRUE, sorted = FALSE)
	    colnames(res) <- c("age", "year", "trajectory", "indicator")
        res <- res[, `:=`(age = ages.all(pop.pred$annual, observed = observed)[as.integer(age)], 
                          year = as.integer(year))]
        data.table::setcolorder(res, "year", "age")
	}
	return(res)
}

get.pop.ex.all <- function(expression, pop.pred, observed=FALSE, as.dt = TRUE, parallel = FALSE, ...){
    # Return trajectories or observed values for an expression for all countries
    if(observed) 
        return(get.pop.observed.from.expression.all.countries(expression, pop.pred, as.dt = as.dt, ...))

    countries.idx <- seq_along(pop.pred$countries$code)
    ncores <- 0
    if(parallel)
        ncores <- getOption("cl.cores", detectCores(logical = FALSE))
    if(ncores > 1) {
        # This can take lots of time. Run it in parallel
        cat('Evaluating expression for all countries in parallel on', ncores, 'cores.\n',
            'Use options(cl.cores = {number_of_cores}) to change the number of cores used.\n')
        cl <- create.pop.cluster(ncores)
        clusterExport(cl, c("pop.pred", "expression"), envir=environment())
        result <- parLapplyLB(cl, countries.idx, function(i) .extract.trajectories.from.expression(i, pop.pred, expression, as.dt = as.dt, ...))
        stopCluster(cl)
    } else { # sequential run
        result <- list()
        for(icntry in countries.idx){
            result[[icntry]] <- .extract.trajectories.from.expression(icntry, pop.pred, expression, as.dt = as.dt, ...)
        }
    }
    if(as.dt) {
        trajectory <- NULL
        resdt <- data.table::rbindlist(result)
        resdt <- data.table::melt(resdt, id.vars = c("country_code", "rn"), value.name = "indicator", variable.name = "trajectory",
                    variable.factor = FALSE)
        data.table::setnames(resdt, "rn", "year")
        resdt[, `:=`(year = as.integer(year), trajectory = as.integer(gsub("V", "", trajectory)))]
    } else resdt <- drop(do.call("abind", c(result, list(along = 0))))
    return(resdt)
}

.extract.trajectories.from.expression <- function(icntry, pop.pred, expression, as.dt = TRUE, ...){
    cntry <- pop.pred$countries$code[icntry]
    expr <- gsub('XXX', as.character(cntry), expression, fixed=TRUE)
    thisres <- get.pop.trajectories.from.expression(expr, pop.pred, ...)
    res <- thisres$trajectories[, thisres$index, drop = FALSE]
    if(as.dt)
        res <- data.table::data.table(res, keep.rownames = TRUE)[, `:=`(country_code = cntry)]
    return(res)
}

get.pop.exba.all <- function(expression, pop.pred, observed=FALSE, as.dt = TRUE, parallel = FALSE, ...){
    # Return trajectories or observed values for an expression (defined by age, i.e. includes {}) for all countries
    if(observed) 
        return(get.pop.observed.from.expression.multiple.age.all.countries(expression, pop.pred, as.dt = as.dt))
    
    countries.idx <- seq_along(pop.pred$countries$code)
    ncores <- 0
    if(parallel)
        ncores <- getOption("cl.cores", detectCores(logical = FALSE))
    if(ncores > 1) {
        # This can take lots of time. Run it in parallel
        cat('Evaluating expression for all countries in parallel on', ncores, 'cores.\n',
            'Use options(cl.cores = {number_of_cores}) to change the number of cores used.\n')
        cl <- create.pop.cluster(ncores)
        clusterExport(cl, c("pop.pred", "expression"), envir=environment())
        result <- parLapplyLB(cl, countries.idx, function(i) .extract.trajectories.from.expression.byage(i, pop.pred, expression, as.dt = as.dt, ...))
        stopCluster(cl)
    } else { # sequential run
        result <- list()
        for(icntry in countries.idx){
            result[[icntry]] <- .extract.trajectories.from.expression.byage(icntry, pop.pred, expression, as.dt = as.dt, ...)
        }
    }
    if(as.dt) {
        age <- NULL
        resdt <- data.table::rbindlist(result)
        colnames(resdt)[1:3] <- c("age", "year", "trajectory")
        resdt[, `:=`(age = ages.all(pop.pred$annual, observed = observed)[as.integer(age)], year = as.integer(year))]
        data.table::setcolorder(resdt, c("country_code", "year"))
    } else resdt <- drop(do.call("abind", c(result, list(along = 0))))
    return(resdt)
}

.extract.trajectories.from.expression.byage <- function(icntry, pop.pred, expression, as.dt = TRUE, ...){
    cntry <- pop.pred$countries$code[icntry]
    expr <- gsub('XXX', as.character(cntry), expression, fixed=TRUE)
    thisres <- get.pop.trajectories.from.expression.multiple.age(expr, pop.pred, ...)
    res <- thisres$trajectories[, , thisres$index, drop = FALSE]
    if(as.dt)
        res <- data.table::as.data.table(res, keep.rownames = TRUE, 
                                        sorted = FALSE, value.name = "indicator")[, `:=`(country_code = cntry)]
    return(res)
}
    
get.pop.trajectories.from.expression <- function(expression, pop.pred, nr.traj=NULL, typical.trajectory=FALSE, 
													adj.to.file=NULL, adj.country=NULL, allow.negative.adj = TRUE, ...) {
	result <- eval(parse(text=.parse.pop.expression(expression, args='...')))
	odim <- length(dim(result))
	ntraj <- dim(result)[odim]
	traj.idx <- NULL
	if(odim == 4 && dim(result)[[1]] == 1) 
		result <- adrop(result, drop=1) # remove country dimension
	if(length(dim(result)) == 3 && dim(result)[[1]] == 1) 
		result <- adrop(result, drop=1) # remove age dimension
	if(ntraj > 1) {
		if(typical.trajectory) {
			traj.idx <- bayesTFR:::get.typical.trajectory.index(result)
		} else {
			thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(result)[2])
			if (thintraj$nr.points > 0) 
		 		traj.idx <- thintraj$index
		}
	} else  {# only 1 trajectory
		traj.idx <- 1
		l<-length(dim(result))
		if(is.null(dim(result)) || (l==2 && dim(result)[2]==1)) along <- 2
		else along <- if(odim > l) l+1 else l 
		result <- abind(result, NULL, along=along)
	}
	if(all(is.na(result[1,]))) { # first time period is NA probably due to using "mid.period"; replace with observed values
	    observed <- get.pop.observed.from.expression(expression, pop.pred, allow.higher.ages = TRUE, ...)
	    result[1,] <- observed[length(observed)]
	}
	if(!is.null(adj.to.file))
		result <- adjust.to.dataset(adj.country, result, adj.file=adj.to.file, use='trajectories',
		                            allow.negatives = allow.negative.adj)
	return(list(trajectories=result, index=traj.idx))
}

get.pop.trajectories.from.expression.multiple.age <- function(expression, pop.pred, nr.traj=NULL, typical.trajectory=FALSE, ...) {
	result <- eval(parse(text=.parse.pop.expression(expression, args='...')))
	odim <- length(dim(result))
	ntraj <- dim(result)[odim]
	traj.idx <- NULL
	country.dropped <- FALSE
	if(odim == 4 && dim(result)[[1]] == 1) {
		result <- adrop(result, drop=1) # remove country dimension
		country.dropped <- TRUE
	}
	if(ntraj > 1) {
		if(typical.trajectory) {
			traj.idx <- bayesTFR:::get.typical.trajectory.index(result[,1,])
		} else {
			thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(result)[3])
			if (thintraj$nr.points > 0) 
		 		traj.idx <- thintraj$index
		}
	} else  {# only 1 trajectory
		traj.idx <- 1
		if(is.null(dim(result))) along <- 3
		else {
			l<-length(dim(result))
			along <- if(odim > l && !country.dropped) l+1 else l 
		}
		result <- abind(result, NULL, along=along)
	}
	if(all(is.na(result[,1,]))) { # first time period is NA probably due to using "mid.period"; replace with observed values
	    observed <- get.pop.observed.from.expression.multiple.age(expression, pop.pred, allow.higher.ages = TRUE, ...)
	    result[1:nrow(observed),1,] <- observed[,ncol(observed)]
	}
	return(list(trajectories=result, index=traj.idx))
}

.get.compatible.pop.matrix.middle.years <- function(pop.pred) {
	years <- colnames(pop.pred$inputs$pop.matrix[['male']])
	if(pop.pred$annual) return(years)
	return(if(is.null(pop.pred$proj.years.pop)) years else as.character(as.integer(years)-2))
}

.get.compatible.pop.matrix.end.years <- function(pop.pred) {
	years <- colnames(pop.pred$inputs$pop.matrix[['male']])
	if(pop.pred$annual) return(years)
	return(if(is.null(pop.pred$proj.years.pop)) as.character(as.integer(years)+2) else years)
}

get.pop.observed.from.expression <- function(expression, pop.pred, as.vector=TRUE, as.dt = FALSE, ...) {
	result <- eval(parse(text=.parse.pop.expression(expression, args='observed=TRUE, ...')))
	if(as.vector && !is.null(result)) {
		result <- drop(result)
		if(is.null(names(result))) {
			l <- length(result)
			end <- ncol(pop.pred$inputs$pop.matrix[['male']])		
			names(result) <- .get.compatible.pop.matrix.middle.years(pop.pred)[(end-l+1):end]
		}
	}
	if(as.dt) {
	    result <- data.table::as.data.table(result, keep.rownames = TRUE, sorted = FALSE)
	    colnames(result) <- c("year", "indicator")   
	}
	return(result)
}

get.pop.observed.from.expression.multiple.age <- function(expression, pop.pred, as.dt = FALSE, ...) {
	result <- adrop(eval(parse(text=.parse.pop.expression(expression, args='observed=TRUE, ...'))), drop=c(1,4)) # remove country and trajectory dimension
	years <- .get.compatible.pop.matrix.end.years(pop.pred)
	if(is.null(dimnames(result)[[2]]) || !all(dimnames(result)[[2]] %in% years)){
	    end.time <- ncol(pop.pred$inputs$pop.matrix[['male']])
	    years <- .get.compatible.pop.matrix.middle.years(pop.pred)
	    dimnames(result)[[2]] <- years[(end.time-dim(result)[2]+1):end.time]
	}
	if(as.dt) {
	    age <- NULL
	    result <- data.table::as.data.table(result, keep.rownames = TRUE, sorted = FALSE)
	    result <- data.table::melt(result, id.vars = "rn", value.name = "indicator", 
	                               variable.name = "year", variable.factor = FALSE)
	    data.table::setnames(result, "rn", "age")
	    result[, `:=`(age = ages.all(pop.pred$annual, observed = TRUE)[as.integer(age)], year = as.integer(year))]
	    data.table::setcolorder(result, "year")
	}
	return(result) 
}

collect.all.countries.observed <- function(expression, pop.pred, time.index = NULL, as.dt = FALSE, ...) {
    years <- .get.compatible.pop.matrix.middle.years(pop.pred)
    years.pop <- .get.compatible.pop.matrix.end.years(pop.pred)
    set.year.dimnames <- FALSE
    result <- list()
    for(icntry in seq_along(pop.pred$countries$code)){
        cntry <- pop.pred$countries$code[icntry]
        expr <- gsub('XXX', as.character(cntry), expression, fixed=TRUE)
        thisres <- get.pop.observed.from.expression(expr, pop.pred, as.vector=FALSE, ...)
        if(icntry == 1){
            if(is.null(dimnames(thisres)[[3]]) || !all(dimnames(thisres)[[3]] %in% years.pop)){
                end.time <- ncol(pop.pred$inputs$pop.matrix[['male']])
                years.to.use <- years[(end.time-dim(thisres)[3]+1):end.time]
                set.year.dimnames <- TRUE
            } else years.to.use <- years.pop[1:dim(thisres)[3]]
            if(is.null(time.index)) time.index <- seq_along(years.to.use)
        }
        if(set.year.dimnames) dimnames(thisres)[[3]] <- years.to.use
        result[[icntry]] <- thisres[1, , years.to.use[time.index],]
        if(as.dt){
            result[[icntry]] <- data.table::as.data.table(result[[icntry]], keep.rownames = TRUE, sorted = FALSE)[, `:=`(country_code = cntry)]
        }
    }
    return(result)
}

get.pop.observed.from.expression.all.countries <- function(expression, pop.pred, time.index = NULL, as.dt = FALSE, ...) {
    result <- collect.all.countries.observed(expression, pop.pred, time.index = time.index, as.dt = as.dt, ...)
    if(as.dt){
        resdt <- data.table::rbindlist(result)
        data.table::setnames(resdt, c("V1", "V2"), c("year", "indicator"))
        resdt[, year := as.integer(year)]
        data.table::setcolorder(resdt, "country_code")
    } else {
        resdt <- do.call("abind", c(result, list(along = 0)))
        dimnames(resdt)[[1]] <- pop.pred$countries$code
    }
	return(resdt)
}

get.pop.observed.from.expression.multiple.age.all.countries <- function(expression, pop.pred, 
                                                                        time.index = NULL, as.dt = FALSE, ...) {
    result <- collect.all.countries.observed(expression, pop.pred, time.index = time.index, as.dt = as.dt, ...)
    
    if(as.dt){
        rn <- NULL
        resdt <- data.table::rbindlist(result)
        resdt <- data.table::melt(resdt, id.vars = c("country_code", "rn"), value.name = "indicator", variable.name = "year",
                      variable.factor = FALSE)
        resdt[, `:=`(age = ages.all(pop.pred$annual, observed = TRUE)[as.integer(rn)],
                     year = as.integer(year))][, rn := NULL]
        data.table::setcolorder(resdt, c("country_code", "year", "age"))
    } else {
        resdt <- do.call("abind", c(result, list(along = 0)))
        dimnames(resdt)[[1]] <- pop.pred$countries$code
    }
    return(resdt)
}


pop.combine <- function(data1, data2, fun, ..., split.along=c('age', 'traj', 'country')) {
	if(length(dim(data1))==length(dim(data2)) && all(dim(data1)==dim(data2))) 
		return(do.call(.remove.trailing.spaces(fun), list(data1, data2)))
	split.along <- match.arg(split.along)
	if(dim(data1)[3] != dim(data2)[3])
		stop('Mismatch in time dimension.', dim(data1)[3], ' vs. ', dim(data2)[3])
	all.splits.mrg <- c(2,4,1) 
	all.splits <- c('age', 'traj', 'country')
	margin <- all.splits.mrg[which(split.along == all.splits)]
	for(imar in 1:length(all.splits.mrg)) {
		if(split.along != all.splits[imar] && dim(data1)[all.splits.mrg[imar]] != dim(data2)[all.splits.mrg[imar]])
			stop('Mismatch in ', all.splits[imar], ' dimension.', 
				dim(data1)[all.splits.mrg[imar]], ' vs. ', dim(data2)[all.splits.mrg[imar]])
	}	
	data <- data2
	unchanged <- data1
	dropd <- c()
	if(split.along != 'country' && dim(data)[1]==1) dropd <- c(dropd, 1)
	if(dim(data)[margin]==1) dropd <- c(dropd, margin)
	data <- adrop(data, drop=dropd)
	return(pop.apply(unchanged, fun, split.along=split.along, data, ...))
	
}

.do.mid.period <- function(data, annual = FALSE) {
    periods <- as.integer(dimnames(data)[[3]])
    if(!annual && periods[1] %% 5 != 0) return(data) # only apply if the time is not the middle of 5-year periods
    newdata <- data[,,-dim(data)[3],, drop = FALSE] + (data[,,-1,, drop = FALSE] - data[,,-dim(data)[3],, drop = FALSE])/2.
    data[] <- NA
    data[,,-1,] <- newdata
    # adjust period names
    if(!annual) dimnames(data)[[3]] <- periods - 2
    data
}

mid.period <- function(data) return(.do.mid.period(data))
mid.period1 <- function(data) return(.do.mid.period(data, annual = TRUE))

period.diff <- function(data) {
    newdata <- data[,,-1,, drop = FALSE] - data[,,-dim(data)[3],, drop = FALSE]
    data[] <- NA
    data[,,-1,] <- newdata
    data
}

period.diff1 <- function(data) return(period.diff(data))

.do.period.ratio <- function(data, annual = FALSE) {
    newdata <- data[,,-1,, drop = FALSE] / data[,,-dim(data)[3],, drop = FALSE]
    data[] <- NA
    data[,,-1,] <- newdata
    # adjust period names
    periods <- as.integer(dimnames(data)[[3]])
    if(!annual && periods[1] %% 5 == 0)
        dimnames(data)[[3]] <- periods - 2
    data
}

period.ratio <- function(data) {
    return(.do.period.ratio(data))
}

period.ratio1 <- function(data) {
    return(.do.period.ratio(data, annual = TRUE))
}

pop.apply <- function(data, fun, ..., split.along=c('None', 'age', 'traj', 'country')) {
	if(is.character(fun)) fun <- .remove.trailing.spaces(fun) 
	split.along <- match.arg(split.along)
	if(split.along == 'None') {
		res <- aaply(data, c(1,3,4), fun, ..., .drop=FALSE)
		# add age dimension
		dim(res) <- c(dim(data)[1], 1, dim(data)[3:4])
	} else {
		margin <- c(2,4,1)[which(split.along == c('age', 'traj', 'country'))]
		res <- aaply(data, margin, fun, ..., .drop=FALSE)
		res <- .add.dropped.dims.and.perm.after.aaply(data, res, margin)
	}
	for(i in c(1,3,4)) if(!is.null(dimnames(data)[[i]])) dimnames(res)[[i]] <- dimnames(data)[[i]]
	return(res)
}

.add.dropped.dims.and.perm.after.aaply <- function(data1, data2, margin) {
	if(length(dim(data2)) < length(dim(data1))) { # this is because sometimes aaply does drop a dimension
			# aaply puts the sliced dimension in front, so we have to find which dimensions were dropped
		i <- j <- 1
		ddim <- dim(data1)[-margin]
		rdim <- c()
		while(i <= length(ddim) && length(rdim) < length(ddim)) { 
			if(ddim[i]==1 && (j > length(dim(data2))-1 || dim(data2)[-1][j] > 1)) {
				rdim <- c(rdim, 1)
			} else {
				rdim <- c(rdim, dim(data2)[-1][j])
				j <- j+1
			}
			i <- i+1
		}
		dim(data2) <- c(dim(data2)[1], rdim)
	}
	if(margin != 1) {# must be rearanged because aaply sometimes puts the sliced dimension in front
		which.not.equal <- which(dim(data1)!=dim(data2))
		if(length(which.not.equal)>1 || (length(which.not.equal)==1 && which.not.equal==1)) {
			pos <- rep(NA, length(dim(data2)))
			pos[margin] <- 1
			pos[is.na(pos)] <- seq(2,length(dim(data2)))
			data2 <- aperm(data2, pos)
		}
	}
	return(data2)
}

.do.gmedian <- function(f, cats=NULL, single = FALSE) {
	if(all(is.na(f))) return(NA)
	# group median
    by <- if(single) 1 else 5
	if(is.null(cats)) cats <- seq(0, by=by, length=length(f)+1)
	nhalf <- sum(f)/2.
	cumsumf <- cumsum(f)
	medcat <- findInterval(nhalf, cumsumf) + 1
	med <- cats[medcat] + ((nhalf-cumsumf[medcat-1])/f[medcat])*(cats[medcat+1]-cats[medcat])
	return(med)
}

gmedian <- function(f, cats=NULL) {
    return(.do.gmedian(f, cats = cats, single = FALSE))
}

gmedian1 <- function(f, cats=NULL) {
    return(.do.gmedian(f, cats = cats, single = TRUE))
}

.do.gmean <- function(f, cats, single = FALSE) {
	if(all(is.na(f))) return(NA)
	# group mean
    by <- if(single) 1 else 5
	if(is.null(cats)) cats <- seq(0, by=by, length=length(f)+1)
	l <- min(length(cats), length(f)+1)
	mid.points <- cats[1:(l-1)] + (cats[2:l] - cats[1:(l-1)])/2.
	counts <- f*mid.points
	return(sum(counts)/sum(f))
}

gmean <- function(f, cats=NULL) {
    return(.do.gmean(f, cats = cats, single = FALSE))
}

gmean1 <- function(f, cats=NULL) {
    return(.do.gmean(f, cats = cats, single = TRUE))
}

.remove.trailing.spaces <- function(x) return(gsub("^[[:blank:]]|[[:blank:]]$", '', x))
.remove.all.spaces <- function(x) return(gsub("[[:blank:]]", '', x))

.do.age.func <- function(data, fun, mid.ages) {
    if(is.character(fun)) fun <- .remove.trailing.spaces(fun) 
    age <- as.integer(dimnames(data)[[2]])
    all.ages <- aperm(array(mid.ages, c(length(mid.ages), dim(data)[[1]],dim(data)[[3]],dim(data)[[4]])), 
                      c(2,1,3,4)) # to assure elementwise operations
    return(do.call(fun, list(data, all.ages[,age,,,drop=FALSE])))
}

age.func <- function(data, fun="*") {
	# data is expected to be 4-d array where the second dimension is age
	# It applies the given function to data and the corresponding age (middle of the age category)
    return(.do.age.func(data, fun, mid.ages = ages.all(annual = FALSE) + 2.5))
}

age.func1 <- function(data, fun="*") {
    # like age.func but for single-year age groups
    return(.do.age.func(data, fun, mid.ages = ages.all(annual = TRUE) + 0.5))
}

drop.age <- function(data) {
	dim(data) <- dim(data)[-2]
	return(data)
}

age.index01 <- function(end) return (c(-1,0,2:end))
age.index05 <- function(end) return (1:end)

mac.expression1 <- function(country) {
    factors <- ages.fert(annual = TRUE) + 0.5
    return(paste0("(", paste0(factors, "*R", country, "[", age.index.fert(annual = TRUE)-1, "]", collapse=" + "), ")/100"))
}

mac.expression5 <- function(country) mac.expression(country)

mac.expression <- function(country) {
	factors <- ages.fert(annual = FALSE) + 2.5
	return(paste0("(", paste0(factors, "*R", country, "[", age.index.fert(annual = FALSE), "]", collapse=" + "), ")/100"))
}

	
.solve.expression.for.country <- function(icountry, pop.pred, expression, adjust=FALSE, ...) {
    country <- pop.pred$countries$code[icountry]
    expr <- gsub('XXX', as.character(country), expression, fixed=TRUE)
    trajectories <- get.pop.trajectories.from.expression(expr, pop.pred, adjust=adjust, ...)
    return(get.pop.traj.quantiles(NULL, pop.pred, icountry, country, 
                                  trajectories=trajectories$trajectories,	q=get.quantiles.to.keep()))
}

.solve.observed.expression.for.country <- function(icountry, pop.pred, expression) {
    country <- pop.pred$countries$code[icountry]
    expr <- gsub('XXX', as.character(country), expression, fixed=TRUE)
    return(get.pop.observed.from.expression(expr, pop.pred))
}

get.pop.from.expression.all.countries <- function(expression, pop.pred, quantiles = NULL, 
                                                  time.index, observed = FALSE, 
                                                  adjust=FALSE, adj.to.file=NULL, 
                                                  allow.negative.adj = TRUE) {
    # get quantiles from expression for all countries
	compressed.expr <- gsub("[[:blank:]]*", "", expression) # remove spaces
	if(observed) compressed.expr <- paste0(compressed.expr, "_observed")
	if(!is.null(adj.to.file)) {
		adjust <- FALSE
		adjdata <- read.table(adj.to.file, header=TRUE, check.names=FALSE)
	}
	.adjust.to.dataset.if.needed <- function(dat, cidx) {
		if(is.null(adj.to.file)) return(dat)
		res <- dat
		for(i in 1:nrow(dat)) {
			cntry <- pop.pred$countries$code[cidx[i]]
			tmp <- adjust.to.dataset(country=cntry, q=dat[rownames(dat)==cntry,], adj.dataset=adjdata, 
			                         years=get.pop.prediction.periods(pop.pred, end.time.only=TRUE)[time.index], use="write",
			                         allow.negatives = allow.negative.adj)
			if(length(tmp)==0) next # no adjustment
			res[i,] <- tmp
		}
		return(res)
	}
	#if(adjust) compressed.expr <- paste0(compressed.expr, '_adjusted')
	if(!is.null(pop.pred$cache) && !is.null(pop.pred$cache[[compressed.expr]])) {
		data <- pop.pred$cache[[compressed.expr]][,,time.index, drop = FALSE]
		if(!observed) 
		    data <- data[,as.character(quantiles),, drop=FALSE]
		if(dim(data)[3] == 1) data <- adrop(data, 3)
		.all.is.na <- function(x) return(all(is.na(x)))
		countries.idx <- which(apply(data, 1, .all.is.na))
		if(length(countries.idx) <= 0) return(.adjust.to.dataset.if.needed(data, 1:nrow(pop.pred$countries)))
	} else {
		countries.idx <- 1:nrow(pop.pred$countries)
		data <- matrix(NA, nrow=dim(pop.pred$quantiles)[1], ncol=if(observed) 1 else length(quantiles))
		rownames(data) <- dimnames(pop.pred$quantiles)[[1]]
		if(!observed) {
		    colnames(data) <- quantiles
		    pop.pred$cache[[compressed.expr]] <- array(NA, dim(pop.pred$quantilesM), dimnames=dimnames(pop.pred$quantilesM))
		} else pop.pred$cache[[compressed.expr]] <- array(NA, 
		                                                  c(dim(pop.pred$quantilesM)[1], 1, dim(pop.pred$inputs$pop.matrix$male)[2]),
		                                                  dimnames=list(dimnames(pop.pred$quantilesM)[[1]], NULL, dimnames(pop.pred$inputs$pop.matrix$male)[[2]]))
	}
	ncores <- getOption("cl.cores", detectCores(logical = FALSE))
	if(ncores > 1 && length(countries.idx)>10) {
		# This can take lots of time. Run it in parallel
		cat('Evaluating expression for all countries in parallel on', ncores, 'cores.\n')
		cl <- create.pop.cluster(ncores)
		clusterExport(cl, c("pop.pred", "expression"), envir=environment())
		if(!observed)
		    quant.list <- parLapplyLB(cl, countries.idx, function(i) .solve.expression.for.country(i, pop.pred, expression, adjust=adjust, allow.negative.adj = allow.negative.adj))
		else 
		    quant.list <- parLapplyLB(cl, countries.idx, function(i) .solve.observed.expression.for.country(i, pop.pred, expression))
		stopCluster(cl)
		for(icountry in countries.idx) {
			pop.pred$cache[[compressed.expr]][icountry,,] <- quant.list[[icountry]]
			data[icountry,] <- if(observed) quant.list[[icountry]][time.index] else quant.list[[icountry]][paste0(quantiles*100, '%'), time.index]
		}
	} else { # run sequentially
		if(length(countries.idx)>10) cat('Evaluating expression for all countries sequentially. Please be patient.\n')
		for(icountry in countries.idx) {
		    if(!observed) {
			    quant <- .solve.expression.for.country(icountry, pop.pred, expression, adjust=adjust, allow.negative.adj = allow.negative.adj)
			    data[icountry,] <- quant[paste0(quantiles*100, '%'), time.index]
			    #stop('')
		    } else {
		        quant <- .solve.observed.expression.for.country(icountry, pop.pred, expression)
		        quant <- abind(quant, along = 2)
		        data[icountry,] <- quant[time.index]
		    }
			pop.pred$cache[[compressed.expr]][icountry,,] <- quant
		}
	}
	.save.cache(pop.pred)
	res <- if(!observed) .adjust.to.dataset.if.needed(data, countries.idx) else data
	return(res)	
}

get.pop.observed.all.countries <- function(pop.pred, time.index, sex='both', age='all') {
	data <- get.pop.observed.multiple.countries(pop.pred, countries=pop.pred$countries$code, sex=sex, 
														age=age, sum.over.ages=TRUE)$data
	return(data[,time.index])
}

get.pop.all.countries <- function(pop.pred, quantiles, projection.index, sex='both', age='all') {
	data <- NULL
	if(sex == 'both') {
		if(age[1]=='all') data <- pop.pred$quantiles[,as.character(quantiles), projection.index]
	} else {
		if(sex=='male') {
			if (age[1]=='all') data <- pop.pred$quantilesM[,as.character(quantiles), projection.index]
			else {if (length(age) == 1) data <- pop.pred$quantilesMage[,age,as.character(quantiles), projection.index]}
		} else {#female
			if (age[1]=='all') data <- pop.pred$quantilesF[,as.character(quantiles), projection.index]
			else {if (length(age) == 1) data <- pop.pred$quantilesFage[,age,as.character(quantiles), projection.index]}
		}
	}
	if(is.null(data)) { # create expression
		expr <- 'PXXX'
		if(sex=='male') expr <- paste(expr, 'M', sep='_')
		if(sex=='female') expr <- paste(expr, 'F', sep='_')
		if (age[1]!='all') expr <- paste0(expr, '[c(', gsub(" ", "", toString(age)), ')]')
		data <- get.pop.from.expression.all.countries(expr, pop.pred, quantiles, projection.index)
	}
	return(data)
}

litem <- function(i, x, default=NULL) { 
	# return element i of the list x if it exists otherwise default
	j <- match(i, names(x)) # this is suppose to be faster than i %in% names(x)
	if (is.na(j)) return(default) 
	x[[i]]
}

as.environment.bayesPop.prediction <- function(x){
	epred <- list2env(x)
	class(epred) <- c(class(x), class(epred))
	return(epred)
}

UNcountries <- function()
	return(UNlocations$country_code[UNlocations$location_type==4])
	
cohorts <- function(pop.pred, country=NULL, expression=NULL, pi=c(80, 95)) {
	.get.quantiles.from.cohort.data <- function(trajs) {
		res <- apply(trajs, 1, median)
		res <- rbind(res, apply(trajs, 1, quantile, quants))
		rownames(res) <- c("median", quants.char)
		res
	}
	if(is.null(country) && is.null(expression))
		stop("Either country or expression must be given.")
	
	if(!is.null(country)) {
		country.object <- get.country.object(country, country.table=pop.pred$countries)
		expression <- paste0("P", country.object$code, "{}")
	}
	alldata <- get.pop.trajectories.from.expression.multiple.age(expression, pop.pred)$trajectories # should be 3d array
	quants <- sort(c((1-pi/100)/2, (1+pi/100)/2))
	quants.char <- as.character(quants)
	result <- list()
	age.index <- as.integer(dimnames(alldata)[[1]])
	nage <- dim(alldata)[1]
	years <- dimnames(alldata)[[2]]
	last.observed.cohort <- pop.pred$proj.years.pop[1]-age.index[1]*5
	from.cohorts <- seq(last.observed.cohort, length=nage-1, by=-5)
	observed.cohorts <- paste(from.cohorts, '-', from.cohorts+5, sep="")
	for(cohort in length(observed.cohorts):1) {
		cohort.traj <- apply(alldata[cohort:nage,,,drop=FALSE], 3, 'diag')
		result[[observed.cohorts[cohort]]] <- .get.quantiles.from.cohort.data(cohort.traj)
		colnames(result[[observed.cohorts[cohort]]]) <- years[1:ncol(result[[observed.cohorts[cohort]]])]
	}
	nyears <- length(pop.pred$proj.years.pop)
	from.cohorts <- seq(last.observed.cohort + 5, length=nyears-2, by=5)
	projected.cohorts <- paste(from.cohorts, '-', from.cohorts+5, sep="")	
	for(cohort in 1:length(projected.cohorts)) {
		cohort.traj <- apply(alldata[,(cohort+1):nyears,,drop=FALSE], 3, 'diag')
		result[[projected.cohorts[cohort]]] <- .get.quantiles.from.cohort.data(cohort.traj)
		colnames(result[[projected.cohorts[cohort]]]) <- years[(cohort+1):(cohort+ncol(result[[projected.cohorts[cohort]]]))]
	}
	#stop('')
	result[["last.observed"]] <- last.observed.cohort
	return(result)
}

get.trajectory.indices <- function(pop.pred, country, what=c("TFR", "e0M", "e0F", "migM", "migF")) {
	# provides indices to the trajectories of TFR, e0 and migration
	country.object <- get.country.object(country, country.table=pop.pred$countries)
	e <- new.env()
	if (!.load.traj.file(pop.output.directory(pop.pred), country.object$code, e))
		return (NULL)
	par <- list(TFR="TFRpred", e0M="e0Mpred", e0F="e0Fpred", migM="migMpred", migF="migFpred")
	return(e$trajectory.indices[[par[[what[1]]]]])
}

extract.trajectories.eq <- function(pop.pred, country=NULL, expression=NULL, quant=0.5, values=NULL, nr.traj=1, ...) {
	# Return trajectories close to the given quantile, or close to the given values, including their index
	if(!is.null(country)) {
		country.object <- get.country.object(country, country.table=pop.pred$countries)
		trajectories <- get.pop.trajectories(pop.pred, country.object$code, ...)$trajectories
	} else {
		trajectories <- get.pop.trajectories.from.expression(expression, pop.pred, ...)$trajectories
	}
	if(is.null(trajectories)) return(NULL)
	if(is.null(values))
		values <- apply(trajectories[2:nrow(trajectories),], 1, quantile, quant, na.rm=TRUE)
	sumerrors <- apply(abs(trajectories[2:nrow(trajectories),] - values), 2, sum)
	sorterrors.idx <- order(sumerrors)
    return(list(trajectories=trajectories[,sorterrors.idx[1:nr.traj]], index=sorterrors.idx[1:nr.traj]))
}

.extract.trajectories <- function(fun, pop.pred, country=NULL, expression=NULL, quant=0.5, 
									values=NULL, all=TRUE, ...) {
	if(!is.null(country)) {
		country.object <- get.country.object(country, country.table=pop.pred$countries)
		trajectories <- get.pop.trajectories(pop.pred, country.object$code, ...)$trajectories
	} else {
		trajectories <- get.pop.trajectories.from.expression(expression, pop.pred, ...)$trajectories
	}
	if(is.null(trajectories)) return(NULL)
	if(is.null(values))
		values <- apply(trajectories[2:nrow(trajectories),], 1, quantile, quant, na.rm=TRUE)
	all.any <- if(all) 'all' else 'any'
	residx <- which(apply(do.call(fun, list(trajectories[2:nrow(trajectories),], values)), 2, all.any))
    return(list(trajectories=trajectories[,residx], index=residx))										
}

extract.trajectories.ge <- function(pop.pred, country=NULL, expression=NULL, quant=0.5, 
									values=NULL, all=TRUE, ...) {
	# Return trajectories greater than the given quantile, or greater than the given values
	return(.extract.trajectories('>=', pop.pred, country, expression, quant, values, all, ...))
}

extract.trajectories.le <- function(pop.pred, country=NULL, expression=NULL, quant=0.5, 
									values=NULL, all=TRUE, ...) {
	# Return trajectories greater than the given quantile, or greater than the given values
	return(.extract.trajectories('<=', pop.pred, country, expression, quant, values, all, ...))
}

