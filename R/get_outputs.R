has.pop.prediction <- function(sim.dir) {
	if(file.exists(file.path(sim.dir, 'predictions', 'prediction.rda'))) return(TRUE)
	return(FALSE)
}

get.pop.prediction <- function(sim.dir, aggregation=NULL) {
	############
	# Returns an object of class bayesPop.prediction
	############
	if(!is.null(aggregation)) return(get.pop.aggregation(sim.dir, name=aggregation))
	output.dir <- file.path(sim.dir, 'predictions')
	return(.get.prediction.object(output.dir))
}

has.pop.aggregation <- function(sim.dir=NULL, pop.pred=NULL, return.dirs=FALSE) {
	if(is.null(sim.dir)) { #remove the last subdirectory
		simdir.split <- strsplit(pop.pred$output.directory, .Platform$file.sep)[[1]]
		sim.dir <- file.path(simdir.split[-length(simdir.split)])
	}
	dirs <- list.files(sim.dir, pattern='^aggregations_', full.names=FALSE)
	if(return.dirs) return(dirs)
	return (length(dirs) > 0)
}

available.pop.aggregations <- function(pop.pred){
	dirs <- has.pop.aggregation(pop.pred=pop.pred, return.dirs=TRUE)
	if(length(dirs)<=0) return(c())
	return(substr(dirs, 14, nchar(dirs)))
}

get.pop.aggregation <- function(sim.dir=NULL, pop.pred=NULL, name=NULL) {
	############
	# Returns an object of class bayesPop.prediction created by aggregation
	############
	if(is.null(sim.dir)) { #remove the last subdirectory
		simdir.split <- strsplit(pop.pred$output.directory, .Platform$file.sep)[[1]]
		l <- length(simdir.split)-1
		lst <- vector('list', l)
		lst[1:l] <- simdir.split[-(l+1)]
		sim.dir <- do.call(file.path, lst)
	}
	dirs <- has.pop.aggregation(sim.dir=sim.dir, return.dirs=TRUE)
	if(length(dirs) == 0) {
		warning('No aggregation available in', sim.dir)
		return(NULL)
	}
	output.dir <- file.path(sim.dir, dirs)
	names <- substr(dirs, 14, nchar(dirs))
	if(length(names) == 1){
		if(!is.null(name) && name != names) 
			warning('Mismatch in aggregation names. Available aggregation is called', names)
		return(.get.prediction.object(output.dir))
	}
	idx <- which(names == name)
	if (length(idx) > 0) return(.get.prediction.object(output.dir[idx]))
	idx <- menu(names, title='Available aggregations:')
	return(.get.prediction.object(output.dir[idx]))
}

.get.prediction.object <- function(directory) {
	pred.file <- file.path(directory, 'prediction.rda')
	if(!file.exists(pred.file)) {
		warning('File ', pred.file, ' does not exist.')
		return(NULL)
	}
	load(file=pred.file)
	bayesPop.prediction$output.directory <- directory
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
		if (x$sex != 'all') cat(' for', x$sex)
		cat(':\n')
		print(x$projections, digits=digits, ...)
	}
}

get.pop.observed.with.age <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all') {
	sex <- match.arg(sex)
	data <- pop.pred$inputs$pop.matrix
	if(sex == 'both') {
		data <- data[['male']][,colnames(data[['male']])] + data[['female']][,colnames(data[['male']])]
	} else data <- data[[sex]]
	country.idx <- grep(paste('^', country, '_', sep=''), rownames(data), value=FALSE)
	data <- data[country.idx,]
	age.idx <- if(age[1]=='all' || age[1]=='psr') 1:nrow(data) else age
	age.idx <- age.idx[age.idx <= nrow(data)]
	return(list(data=data, age.idx=age.idx))
}


get.pop.observed <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all', sum.over.ages=TRUE) {
	data.age <- get.pop.observed.with.age(pop.pred, country, sex, age)
	data <- data.age$data
	age.idx <- data.age$age.idx
	if(age[1]=='psr')  # potential support ratio
		return(colSums(data[get.psr.nominator.index(),])/colSums(data[get.psr.denominator.startindex():nrow(data),]))
	if(sum.over.ages) return(colSums(data[age.idx,]))
	return(data[age.idx,])
}

get.psr.nominator.index <- function() return(5:13)
get.psr.denominator.startindex <- function() return(14)

get.pop.trajectories <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all',
 									nr.traj=NULL, typical.trajectory=FALSE) {
	traj.file <- file.path(pop.pred$output.directory, paste('totpop_country', country, '.rda', sep=''))
	quant <- hch <- age.idx <- traj <- traj.idx <-  NULL
	load.traj <- is.null(nr.traj) || nr.traj > 0 || typical.trajectory
	if (!file.exists(traj.file)) 
		return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch))
	load(traj.file)
	sex <- match.arg(sex)
	max.age <- dim(totpf)[1] # should be 27
	age.idx <- if(age[1]=='all' || age[1]=='psr') 1:max.age else age
	if(max(age.idx) > 27 || min(age.idx) < 1) stop('Age index must be between 1 (age 0-4) and 27 (age 130+).')
	if(sex == 'both' && age[1]=='all') {
		if(load.traj) traj <- totp
		quant <- pop.pred$quantiles
		hch <- totp.hch
	} else {
		if (age[1] == 'psr') { # potential support ratio
			if(sex == 'both' && load.traj)
				traj <- (colSums(totpm[get.psr.nominator.index(),,,drop=FALSE]) + 
								colSums(totpf[get.psr.nominator.index(),,,drop=FALSE]))/(
								colSums(totpm[get.psr.denominator.startindex():max.age,,,drop=FALSE]) + 
									colSums(totpf[get.psr.denominator.startindex():max.age,,,drop=FALSE]))
			if(sex == 'male' && load.traj) 
				traj <- colSums(totpm[get.psr.nominator.index(),,,drop=FALSE])/colSums(
											totpm[get.psr.denominator.startindex():max.age,,,drop=FALSE])
			if(sex == 'female' && load.traj) 
				traj <- colSums(totpf[get.psr.nominator.index(),,,drop=FALSE])/colSums(
											totpf[get.psr.denominator.startindex():max.age,,,drop=FALSE])
		} else {
			if(sex == 'both') {
				if(load.traj) traj <- colSums(totpm[age.idx,,,drop=FALSE]) + colSums(totpf[age.idx,,,drop=FALSE])
				hch <- colSums(totpm.hch[age.idx,,,drop=FALSE]) + colSums(totpf.hch[age.idx,,,drop=FALSE])
			} else {
				if(sex=='male') {
					if(load.traj) traj <- colSums(totpm[age.idx,,,drop=FALSE])
					hch <- colSums(totpm.hch[age.idx,,,drop=FALSE])
					if (length(age.idx) == max.age) quant <- pop.pred$quantilesM
					else {if (length(age.idx) == 1) quant <- pop.pred$quantilesMage[,age.idx,,]}
				} else { # female
					if(load.traj) traj <- colSums(totpf[age.idx,,,drop=FALSE])
					hch <- colSums(totpf.hch[age.idx,,,drop=FALSE])
					if (length(age.idx) == max.age) quant <- pop.pred$quantilesF
					else {if (length(age.idx) == 1) quant <- pop.pred$quantilesFage[,age.idx,,]}
				}
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
	 	rownames(traj) <- pop.pred$proj.years
	return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch))
}

get.pop.trajectories.multiple.age <- function(pop.pred, country, sex=c('both', 'male', 'female'), 
												age='all', nr.traj=NULL, proportion=FALSE) {
	# Like get.pop.trajectories() but it doesn't sum up over ages and it doesn't return quantiles
	# Called when creating pop pyramid. Doesn't handle potential support ratio.
	traj.file <- file.path(pop.pred$output.dir, paste('totpop_country', country, '.rda', sep=''))
	quant <- NULL
	age.idx <- NULL
	if (file.exists(traj.file)) {
		load(traj.file)
		sex <- match.arg(sex)
		max.age <- dim(totpm)[1] # should be 27
		age.idx <- if(age[1]=='all') 1:max.age else age
		if(sex == 'both') 
			traj <- totpm[age.idx,,,drop=FALSE] + totpf[age.idx,,,drop=FALSE]
		else {
			traj <- if(sex=='male') totpm[age.idx,,,drop=FALSE] else totpf[age.idx,,,drop=FALSE]

			if(proportion) {
				totpop <- (apply(totpm[,,,drop=FALSE], c(2,3), sum) + apply(totpf[,,,drop=FALSE], c(2,3), sum))
				for(iage in 1:dim(traj)[1])
					traj[iage,,] <- traj[iage,,]/totpop
			}
		}
		thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(traj)[3])
		if (thintraj$nr.points == 0) return(list(trajectories=NULL))
		traj.idx <- thintraj$index
	} else {
		traj <- NULL
		traj.idx <- NULL
	}	
	if(!is.null(traj)) 
	 	dimnames(traj)[[2]] <- pop.pred$proj.years
	return(list(trajectories=traj, index=traj.idx, age.idx=age.idx))
}

is.saved.pi <- function(pop.pred, pi, warning=TRUE) {
	if(length(pi) == 0) return(NULL)
	is.valid.pi <- rep(NA, length(pi))
	quantile.values <- as.numeric(dimnames(pop.pred$quantiles)[[2]])
	for (i in 1:length(pi)) {
		al <- 1-(1-pi[i]/100)/2		
		is.valid.pi[i] <- any(round(quantile.values,6)==round(al,6))
		if(!is.valid.pi[i] && warning)
			warning(pi[i], '% interval not available.')
	}
	return(is.valid.pi)
}

get.pop.traj.quantiles <- function(quantile.array, pop.pred, country.index, country.code, 
									trajectories=NULL, pi=80, q=NULL, reload=TRUE, ...) {
	al <- if(!is.null(q)) q else c((1-pi/100)/2, 1-(1-pi/100)/2)
	found <- FALSE
	if(!is.null(quantile.array)) {
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
		reload <- FALSE
		if (is.null(trajectories)) {
			if(pop.pred$nr.traj > 0) reload <- TRUE
		} else { 
			if (dim(trajectories)[2] < 2000 && pop.pred$nr.traj > dim(trajectories)[2]) reload <- TRUE
		}
		if(reload) {
			#load 2000 trajectories maximum for computing quantiles
			traj.reload <- get.pop.trajectories(pop.pred, country.code, nr.traj=2000, ...)
			trajectories <- traj.reload$trajectories
		}
		if (!is.null(trajectories)) {
			cqp <- apply(trajectories, 1, 
						quantile, al, na.rm = TRUE)
		} else {
			warning('Quantiles not found')
			return(NULL)
		}
	}
	return(cqp)
}

get.popVE.trajectories.and.quantiles <- function(pop.pred, country, event=c('births', 'deaths'), 
									sex=c('both', 'male', 'female'), age='all', sum.over.ages=TRUE,
 									nr.traj=NULL, typical.trajectory=FALSE) {
 	# get trajectories and quantiles for vital events
	traj.file <- file.path(pop.pred$output.directory, paste('vital_events_country', country, '.rda', sep=''))
	quant <- hch <- age.idx <- traj <- traj.idx <-  NULL
	if (!file.exists(traj.file)) 
		return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch))
	load(traj.file)
	sex <- match.arg(sex)
	event <- match.arg(event)
	alltraj <- if(event == 'births') list(male=btm, female=btf, male.hch=btm.hch, female.hch=btf.hch)
				else list(male=deathsm, female=deathsf, male.hch=deathsm.hch, female.hch=deathsf.hch)
	max.age <- dim(alltraj$male)[1] # should be 7 or 27
	age.idx <- age.idx.raw  <- if(age[1]=='all') 1:max.age else age	
	quantiles <- get.quantiles.to.keep()
	if(event == 'births') {
		if(age[1] != 'all') {
			age.idx <- age.idx - 3 # translate age index into mother's child-bearing age index
			if(max(age.idx) > max.age || min(age.idx) < 1) 
				stop('Age index for births must be between 4 (age 15-19) and 10 (age 45-49).')
		} else age.idx.raw <- age.idx + 3
	} 
	if(max(age.idx) > 27 || min(age.idx) < 1) stop('Age index must be between 1 (age 0-4) and 27 (age 130+).')
	
	if(sex == 'both') {
		if(sum.over.ages) {
			traj <- colSums(alltraj$male[age.idx,,,drop=FALSE]) + colSums(alltraj$female[age.idx,,,drop=FALSE])
			hch <- colSums(alltraj$male.hch[age.idx,,,drop=FALSE]) + colSums(alltraj$female.hch[age.idx,,,drop=FALSE])
		} else {
			traj <- alltraj$male[age.idx,,,drop=FALSE] + alltraj$female[age.idx,,,drop=FALSE]
			hch <- alltraj$male.hch[age.idx,,,drop=FALSE] + alltraj$female.hch[age.idx,,,drop=FALSE]
		}
	} else { # one sex
		if(sum.over.ages) {
			traj <- colSums(alltraj[[sex]][age.idx,,,drop=FALSE])
			hch <- colSums(alltraj[[paste(sex,'hch', sep='.')]][age.idx,,,drop=FALSE])
		} else {
			traj <- alltraj[[sex]][age.idx,,,drop=FALSE]
			hch <- alltraj[[paste(sex,'hch', sep='.')]][age.idx,,,drop=FALSE]
		}
	}
	if(sum.over.ages) { # quantiles are 2-d arrays
		quant <- apply(traj, 1, quantile, quantiles, na.rm = TRUE)
		dimnames(quant) <- list(quantiles, dimnames(alltraj$male)[[2]])
		traj.for.thinning <- traj
		year.dim <- 1
	} else { # quantiles are 3-d arrays: age x quantiles x period
		quant <- aperm(apply(traj, c(1,2), quantile, quantiles, na.rm = TRUE), c(2,1,3))
		dimnames(quant) <- list(pop.pred$ages[age.idx.raw], quantiles, dimnames(alltraj$male)[[2]])
		traj.for.thinning <- traj[1,,]
		year.dim <- 2
	}
	if(is.null(nr.traj) || nr.traj > 0) {
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
		dimnames(traj)[[year.dim]] <- dimnames(alltraj$male)[[2]]
	} else traj <- NULL
	 	
	return(list(trajectories=traj, index=traj.idx, quantiles=quant, 
				age.idx=age.idx, age.idx.raw=age.idx.raw, half.child=hch, event=event))
}


get.age.labels <- function(ages, collapsed=FALSE, age.is.index=FALSE) {
	all.ages <- c(seq(0, by=5, length=27), NA)
	ages.idx <- if(age.is.index) ages else which(is.element(all.ages, ages))
	ages.idx.shift <- ages.idx+1
	if(collapsed) {
		ages.idx.dif <- which(!is.element(ages.idx, ages.idx.shift))
		ages.idx.shift <- ages.idx.shift[!is.element(ages.idx.shift, ages.idx)]
		ages.idx <- ages.idx[ages.idx.dif]
	}
	lages <- all.ages[ages.idx]
	uages <- all.ages[ages.idx.shift]
	l <- length(lages)
	result <- paste(all.ages[ages.idx[1:(l-1)]], '-', all.ages[ages.idx.shift[1:(l-1)]]-1, sep='')
	if (l > 1) result <- c(result, if(is.na(all.ages[ages.idx.shift[l]])) paste(all.ages[ages.idx[l]], '+', sep='')
			else paste(all.ages[ages.idx[l]], '-', all.ages[ages.idx.shift[l]]-1, sep=''))
	return(result)
}	

.get.year.index <- function(year, years) {
	lyears <- length(years)
	res <- as.integer(cut(year, labels=1:lyears, breaks=c(years-3, years[lyears]+2)))
	return(res)
	
	#breaks <- c(years-3, years[lyears]+2)
	#h <- try(hist(year, breaks=breaks, plot=FALSE)$count, silent=TRUE)
	#return(if(inherits(h, "try-error")) NULL else which(h > 0)[1])
}
get.pop.prediction.periods <- function(pop.pred) {
	return(sapply(lapply(pop.pred$proj.years, '+', c(-3, 2)), paste, collapse='-'))
}
get.prediction.year.index <- function(pop.pred, year) {
	years <- pop.pred$proj.years
	return(.get.year.index(year, years))
}

get.observed.year.index <- function(pop.pred, year) {
	years <- as.integer(colnames(pop.pred$inputs$pop.matrix[['male']]))
	return(.get.year.index(year, years))
}
get.pop.observed.periods <- function(pop.pred) {
	return(sapply(lapply(as.integer(colnames(pop.pred$inputs$pop.matrix[['male']])), '+', c(-3, 2)), paste, collapse='-'))
}

get.predORobs.year.index <- function (pred, year) 
{
    projection.index <- get.prediction.year.index(pred, year)
    projection <- TRUE
    if (is.null(projection.index)) {
        projection <- FALSE
        projection.index <- get.observed.year.index(pred, year)
    }
    return(c(index = projection.index, is.projection = projection))
}

get.quantiles.to.keep <- function() 
	return(c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1))

get.countries.table.bayesPop.prediction <- function(object, ...) 
	return(object$countries)
	
get.pop <- function(object, pop.pred, aggregation=NULL, observed=FALSE, ...) {
	split.object <- strsplit(object, '_', fixed=TRUE)[[1]]
	country.part.idx <- grep('^C', split.object)
	if(length(country.part.idx) <= 0) stop('No country specified.')
	country.string <- regmatches(split.object[country.part.idx], 
						regexpr('^C[[:digit:]]*\\[|^C[[:digit:]]*\\{?', split.object[country.part.idx])) 
	country.code <- as.integer(gsub('C|\\[|\\{', '', country.string))
	split.object[country.part.idx] <- gsub(paste('C', country.code, sep=''), '', split.object[country.part.idx], fixed=TRUE)
	if(nchar(split.object[country.part.idx])<=0) split.object <- split.object[-country.part.idx]
	sex <- 'both'
	if(length(split.object) > 0) {
		for(sx in c('F', 'M')) {
			sx.idx <- grep(sx, split.object)
			if(length(sx.idx) > 0) {
				sex <- list(F='female', M='male')[[sx]]
				split.object[sx.idx] <- gsub(sx, '', split.object[sx.idx])
				break
			}
		}
	}
	sum.over.ages <- TRUE
	age <- 'all'
	age.part.idx <- grep("\\[|\\{", split.object)
	if(length(age.part.idx) > 0) {
		if(length(age.part.idx) > 1) stop('Only one age vector is allowed.')
		age <- eval(parse(text=gsub('\\[|\\]|\\{|\\}', '', split.object[age.part.idx])))
		if(is.null(age)) age <- 'all'
		if(grepl("{", split.object, fixed=TRUE)) sum.over.ages <- FALSE
	}
	country.object <- get.country.object(country.code, country.table=pop.pred$countries)
	if(is.null(country.object$code)) {
		av.aggrs <- available.pop.aggregations(pop.pred)
		indep.idx <- which(is.element('independence', av.aggrs))
		if(length(indep.idx) > 0)  # put independence aggregation first
			av.aggrs <- c('independence', av.aggrs[-indep.idx])
		if(is.null(aggregation)) aggregation <- av.aggrs
		for(aggr in aggregation) {
			if(!is.element(aggr, av.aggrs)) {warning('Aggregation', aggr, 'not available.'); next}
			aggr.obj <- get.pop.aggregation(pop.pred=pop.pred, name=aggr)
			country.object <- get.country.object(country.code, country.table=aggr.obj$countries)
			if(!is.null(country.object$code)) {pop.pred <- aggr.obj; break}
		}
	}
	if(observed) {
		traj <- get.pop.observed.with.age(pop.pred, country=country.object$code, sex=sex, age=age)
		d <- traj$data[traj$age.idx,]
		if(sum.over.ages) d <- colSums(d)
		data <- as.matrix(d)
		if(!is.vector(d)) {# only if it was not summed up, because then the as.matrix command adds an dimension
			if(age[1]=='all') { # extend to 27 age categories
				data <- rbind(data, matrix(0, nrow=27-nrow(data), ncol=ncol(data)))
				traj$age.idx <- c(traj$age.idx, (nrow(d)+1):27)
			}
			dim(data) <- c(dim(data), 1)
			dimnames(data)[2:length(dim(d))] <- dimnames(d)[2:length(dim(d))]
		}
	} else {
		traj <- .get.trajectories(sum.over.ages=sum.over.ages, pop.pred, country=country.object$code, sex=sex, age=age, ...)
		data <- traj$trajectories
	}
	if(length(dim(data)) > 2) dimnames(data)[[1]] <- traj$age.idx
	return(data)
}

.get.trajectories <- function(sum.over.ages=TRUE, ...){
	traj <- if(sum.over.ages) get.pop.trajectories(...) else get.pop.trajectories.multiple.age(...)
	return(traj)
}

.parse.pop.expression <- function(expression, args='...') 
	return (gsub('(C[[:graph:]]*[[:alnum:]]|C[[:graph:]]*\\]|C[[:graph:]]*\\})', 
			paste("get.pop('\\1', pop.pred,", args, ")"), expression))
	
get.pop.trajectories.from.expression <- function(expression, pop.pred, nr.traj=NULL, typical.trajectory=FALSE, ...) {
	# find parts that start with C and end either with number, letter or closed bracket 
	new.expression <- .parse.pop.expression(expression)
	result <- eval(parse(text=new.expression))
	traj.idx <- NULL
	if(typical.trajectory) {
		traj.idx <- bayesTFR:::get.typical.trajectory.index(result)
	} else {
		thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(result)[2])
		if (thintraj$nr.points > 0) 
		 	traj.idx <- thintraj$index
	}
	return(list(trajectories=result, index=traj.idx))
}

get.pop.observed.from.expression <- function(expression, pop.pred, ...) {
	result <- as.vector(eval(parse(text=.parse.pop.expression(expression, args='observed=TRUE, ...'))))
	names(result) <- colnames(pop.pred$inputs$pop.matrix[['male']])
	return(result)
}

gmedian <- function(f, cats) {
	nhalf <- sum(f)/2.
	cumsumf <- cumsum(f)
	medcat <- findInterval(nhalf, cumsumf) + 1
	med <- cats[medcat] + ((nhalf-cumsumf[medcat-1])/f[medcat])*(cats[medcat+1]-cats[medcat])
	return(med)
}

age.func <- function(data, func="*") {
	# data is expected to be 2- or 3-d array where the first dimension is age
	# It applies the given function to data and the corresponding age (middle of the age category)
	age <- as.integer(dimnames(data)[[1]])
	all.ages <- seq(2, by=5, length=27)
	return(do.call(func, list(data, all.ages[age])))
}