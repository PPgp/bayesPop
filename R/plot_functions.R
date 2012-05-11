get.pop.observed <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all') {
	sex <- match.arg(sex)
	data <- pop.pred$inputs$pop.matrix
	if(sex == 'both') {
		data <- data[['male']][,colnames(data[['male']])] + data[['female']][,colnames(data[['male']])]
	} else data <- data[[sex]]
	country.idx <- grep(paste('^', country, '_', sep=''), rownames(data), value=FALSE)
	data <- data[country.idx,]
	if(age[1]=='psr')  # potential support ratio
		return(colSums(data[get.psr.nominator.index(),])/colSums(data[get.psr.denominator.startindex():nrow(data),]))
	age.idx <- if(age[1]=='all') 1:nrow(data) else age
	return(colSums(data[age.idx,]))
}

get.psr.nominator.index <- function() return(5:13)
get.psr.denominator.startindex <- function() return(14)

get.pop.trajectories <- function(pop.pred, country, sex=c('both', 'male', 'female'), age='all',
 									nr.traj=NULL, typical.trajectory=FALSE) {
	traj.file <- file.path(pop.pred$output.directory, paste('totpop_country', country, '.rda', sep=''))
	quant <- hch <- age.idx <- traj <- traj.idx <-  NULL
	load.traj <- is.null(nr.traj) || nr.traj > 0
	if (!file.exists(traj.file)) 
		return(list(trajectories=traj, index=traj.idx, quantiles=quant, age.idx=age.idx, half.child=hch))
	load(traj.file)
	sex <- match.arg(sex)
	max.age <- dim(totp)[1] # should be 27
	age.idx <- if(age[1]=='all') 1:max.age else age
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
		max.age <- dim(totp)[1] # should be 27
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

get.age.labels <- function(ages, collapsed=FALSE) {
	all.ages <- c(seq(0, by=5, length=27), NA)
	ages.idx <- which(is.element(all.ages, ages))
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

pop.trajectories.plot.all <- function(pop.pred, 
									output.dir=file.path(getwd(), 'pop.trajectories'),
									output.type="png", verbose=FALSE, ...) {
	# plots pop trajectories for all countries
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.countries <- pop.pred$countries[,'code']
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	for (country in all.countries) {
		country.obj <- get.country.object(country, country.table=pop.pred$countries)
		if(verbose)
			cat('Creating population graph for', country.obj$name, '(', country.obj$code, ')\n')

		do.call(output.type, list(file.path(output.dir, 
										paste('pop.plot_c', country.obj$code, '.', postfix, sep=''))))
		pop.trajectories.plot(pop.pred, country=country.obj$code, ...)
		dev.off()
	}
	if(verbose)
		cat('\nTrajectory plots stored into', output.dir, '\n')
}


pop.trajectories.plot <- function(pop.pred, country, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), age='all',
								  sum.over.ages=FALSE,
								  half.child.variant=FALSE,
								  nr.traj=NULL, typical.trajectory=FALSE,
								  xlim=NULL, ylim=NULL, 
								  xlab='Year', ylab='Population projection', main=NULL,
								  dev.ncol=5, lwd=c(2,2,2,2,1), col=c('black', 'red', 'red', 'blue', 'gray'),
								  show.legend=TRUE, ...
								  ) {
	# lwd is a vector of 5 line widths for: 
	#	1. observed data, 2. median, 3. quantiles, 4. half child variant, 5. trajectories
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	if(length(lwd) < 5) {
		lwd <- rep(lwd, 5)
		lwd[5] <- 1
	}

	country <- get.country.object(country, country.table=pop.pred$countries)
	if(sum.over.ages || age[1]=='psr')
		do.pop.trajectories.plot(pop.pred, country, pi=pi, sex=sex, age=age,
									half.child.variant=half.child.variant, nr.traj=nr.traj,
									typical.trajectory=typical.trajectory,
									xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, lwd=lwd, col=col,
									show.legend=show.legend, ...)
	else {
		all.ages <- pop.pred$ages
		if(age=='all') age <- 1:20
		age.labels <- get.age.labels(pop.pred$ages)
		if(is.null(main)) {
			main <- country$name
			sex <- match.arg(sex) 
			if(sex != 'both') main <- paste(main, ': ', sex, sep='')
		}
		age.labels <- get.age.labels(pop.pred$ages)
		cur.mgp <- par('mgp')
		nplots <- length(age)
		if (nplots < dev.ncol) {
        	ncols <- nplots
			nrows <- 1
        } else {
			ncols <- dev.ncol
			nrows <- ceiling(nplots/dev.ncol)
        }		
		par(mfrow=c(nrows,ncols),  oma = c(0, 0, 2, 0))
		par(mar=c(2,2,1,0.4)+0.1, mgp=c(1,0.3,0))
		for(iage in age) {
			do.pop.trajectories.plot(pop.pred, country, pi=pi, sex=sex, age=iage,
									half.child.variant=half.child.variant, nr.traj=nr.traj,
									typical.trajectory=typical.trajectory,
									xlim=xlim, ylim=ylim, xlab='', ylab='', main=age.labels[iage], cex.main=0.9, 
									lwd=lwd, col=col, show.legend=show.legend, ...)
		}
		mtext(main, line = 0.5, outer = TRUE)
		par(mgp=cur.mgp)
	}
}

do.pop.trajectories.plot <- function(pop.pred, country, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), age='all',
								  half.child.variant=FALSE,
								  nr.traj=NULL, typical.trajectory=FALSE,
								  xlim=NULL, ylim=NULL, type='b', 
								  xlab='Year', ylab='Population projection', main=NULL, 
								  lwd=c(2,2,2,2,1), col=c('black', 'red', 'red', 'blue', 'gray'),
								  show.legend=TRUE, ...
								  ) {

	sex <- match.arg(sex)
	trajectories <- get.pop.trajectories(pop.pred, country$code, sex, age, nr.traj, 
										typical.trajectory=typical.trajectory)
	cqp <- list()
	for (i in 1:length(pi))
		cqp[[i]] <- get.pop.traj.quantiles(trajectories$quantiles, pop.pred, country$index, country$code, 
										trajectories=trajectories$trajectories, pi=pi[i], sex=sex, age=age)
	pop.observed <- get.pop.observed(pop.pred, country$code, sex=sex, age=age)
	obs.not.na <- !is.na(pop.observed)
	pop.observed <- if(sum(obs.not.na)==0) pop.observed[length(pop.observed)] else pop.observed[obs.not.na]
	x1 <- as.integer(names(pop.observed))
	x2 <- as.numeric(dimnames(pop.pred$quantiles)[[3]])
	y1 <- pop.observed
	if(is.null(xlim)) xlim <- c(min(x1, x2), max(x1, x2))
	if(is.null(ylim)) 
		ylim <- c(min(y1, if (!is.null(trajectories$trajectories))
							trajectories$trajectories[,trajectories$index]
						  else NULL, 
						  sapply(cqp, min, na.rm=TRUE), na.rm=TRUE), 
				  max(y1, if (!is.null(trajectories$trajectories))
				  			trajectories$trajectories[,trajectories$index] else NULL, 
				  		  sapply(cqp, max, na.rm=TRUE), na.rm=TRUE))
	if(is.null(main)) {
		main <- country$name 
		if(sex != 'both') main <- paste(main, ': ', sex, sep='')
		if(age[1] == 'psr') main <- paste(main, ' (Potential Support Ratio)', sep='')
		else {
			if(age[1] != 'all') {
				age.labels <- get.age.labels(pop.pred$ages[age], collapse=TRUE)
				main <- paste(main, ' (Age ', paste(age.labels, collapse=','), ')', sep='')
			}
		}
	}
	# plot historical data: observed
	plot(x1, y1, type=type, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
			panel.first = grid(), lwd=lwd[1], col=col[1], ...)
	# plot trajectories
	if(!is.null(trajectories$trajectories)) {
		for (i in 1:length(trajectories$index)) {
			lines(x2, trajectories$trajectories[,trajectories$index[i]], type='l', col=col[5], lwd=lwd[5])
		}
	}
	# plot median
	pop.median <- get.pop.traj.quantiles(trajectories$quantiles, pop.pred, country$index, country$code, 
										trajectories=trajectories$trajectories, q=0.5, sex=sex, age=age)
	lines(x2, pop.median, type='l', col=col[2], lwd=lwd[2]) 
	
	# plot given CIs
	lty <- 2:(length(pi)+1)
	for (i in 1:length(pi)) {		
		if (!is.null(cqp[[i]])) {
			lines(x2, cqp[[i]][1,], type='l', col=col[3], lty=lty[i], lwd=lwd[3])
			lines(x2, cqp[[i]][2,], type='l', col=col[3], lty=lty[i], lwd=lwd[3])
		}
	}
	legend <- c('median', paste(pi, '% PI', sep=''))
	lty <- c(1, lty)
	lwds <- c(lwd[2], rep(lwd[3], length(pi)))
	cols <- c(col[2], rep(col[3], length(pi)))
	if(sum(obs.not.na)>0) {
		legend <- c(legend, 'observed')
		lty <- c(lty, 1)
		lwds <- c(lwds, lwd[1])
		cols <- c(cols, col[1])
	}
	if (half.child.variant && !is.null(trajectories$half.child)) {
		lty <- c(lty, max(lty)+1)
		llty <- length(lty)
		lines(x2, trajectories$half.child[,1], type='l', col=col[4], lty=lty[llty], lwd=lwd[4])
		lines(x2, trajectories$half.child[,2], type='l', col=col[4], lty=lty[llty], lwd=lwd[4])
		legend <- c(legend, '+/- 0.5 child')
		cols <- c(cols, col[4])
		lwds <- c(lwds, lwd[4])
	}
	if(show.legend)
		legend('topleft', legend=legend, lty=lty, bty='n', col=cols, lwd=lwds)
}


pop.trajectories.table <- function(pop.pred, country, pi=c(80, 95),
								  sex=c('both', 'male', 'female'), age='all',
								  half.child.variant=FALSE) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	country <- get.country.object(country, country.table=pop.pred$countries)
	max.age.idx <- 27
	x <- pop.pred$proj.years
	sex <- match.arg(sex)
	l <- length(pop.pred$proj.years)
	pred.table <- matrix(NA, ncol=2*length(pi)+1, nrow=l)
	quant <- NULL
	if (age[1]=='all') age.idx <- 1:max.age.idx
	else {
		if(all(is.element(1:max.age.idx, age))) age.idx <- 1:max.age.idx
		else age.idx <- unique(age)
	}
	lage <- length(age.idx)
	if(lage==max.age.idx) {
		if(sex == 'both') quant <- pop.pred$quantiles
		else quant <- if(sex=='male') pop.pred$quantilesM else pop.pred$quantilesF
	}
	pred.table[,1] <- get.pop.traj.quantiles(quant, pop.pred, country$index, country$code, 
												q=0.5, sex=sex, age=age.idx)
	colnames(pred.table) <- c('median', rep(NA,ncol(pred.table)-1))
	idx <- 2
	for (i in 1:length(pi)) {
		cqp <- get.pop.traj.quantiles(quant, pop.pred, country$index, country$code, 
										pi=pi[i], sex=sex, age=age.idx)
		if (!is.null(cqp)) {
			pred.table[,idx:(idx+1)] <- t(cqp)
		} else{
			pred.table[,idx:(idx+1)] <- matrix(NA, nrow=l, ncol=2)
		}
		al <- (1-pi[i]/100)/2
		colnames(pred.table)[idx:(idx+1)] <- c(al, 1-al)
		idx <- idx+2
	}
	rownames(pred.table) <- x
	cn <- colnames(pred.table)[2:ncol(pred.table)]
	pred.table[,2:ncol(pred.table)] <- pred.table[,cn[order(cn)]]
	colnames(pred.table)[2:ncol(pred.table)] <- cn[order(cn)]
	if(half.child.variant) {
		# load the half child variants from trajectory file
		traj <- get.pop.trajectories(pop.pred, country$code, sex, age, nr.traj=0)
		if(!is.null(traj$half.child)) {
			pred.table <- cbind(pred.table, traj$half.child)
			colnames(pred.table)[(ncol(pred.table)-1):ncol(pred.table)] <- c('-0.5child', '+0.5child')
		}
	}
	return(pred.table)
}

pop.pyramid.all <- function(pop.pred, year=NULL,
									output.dir=file.path(getwd(), 'pop.pyramid'),
									output.type="png", verbose=FALSE, ...) {
	# plots pyramid for all countries and all years given by 'year'
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.countries <- pop.pred$countries[,'name']
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	if(is.null(year)) year <- pop.pred$present.year
	for (country in all.countries) {
		country.obj <- get.country.object(country, country.table=pop.pred$countries)
		if(verbose)
			cat('Creating pyramid(s) for', country, '(', country.obj$code, ')\n')

		for(y in year) {
			do.call(output.type, list(file.path(output.dir, 
										paste('pyr', y, '_c', country.obj$code, '.', postfix, sep=''))))
			pop.pyramid(pop.pred, country=country.obj$code, year=y, ...)
			dev.off()
		}
	}
	if(verbose)
		cat('\nPyramids stored into', output.dir, '\n')
}

.get.data.for.pyramid <- function(pop.pred, country, year, pi, proportion, age, draw.past.year) {
	country <- get.country.object(country, country.table=pop.pred$countries)
	year.idx <- c(if(is.null(year)) 1 else get.prediction.year.index(pop.pred, year), 1)
	draw.projection <- c(TRUE, TRUE)
	draw.observed <- FALSE
	pop.observed <- NULL
	if(draw.past.year) draw.observed <- TRUE
	if(is.na(year.idx)[1] || draw.past.year>1) {
		pop.observed <- get.pop.observed(pop.pred, country$code, sex='both')
		if(is.na(year.idx[1])) {
			year.idx[1] <- get.observed.year.index(pop.pred, year)
			if(is.na(year.idx[1])  || is.na(pop.observed[year.idx[1]])) stop('Unable to find data for year ', year)
			draw.projection[1] <- FALSE
		}
		if(draw.past.year>1) {
			year.idx[2] <- get.observed.year.index(pop.pred, draw.past.year)
			if(is.na(year.idx[2]) || is.na(pop.observed[year.idx[2]])) {
				warning('Unable to find data for year ', draw.past.year)
				draw.observed <- FALSE
			} else draw.projection[2] <- FALSE
		}
	}
	ages.idx <- age[age <=  length(pop.pred$ages)]
	lages <- length(ages.idx)
	pop.median <- pop.present <- list(male=rep(NA, lages), female=rep(NA, lages))
	pop.quant <- list(male=list(), female=list())
	nquant <- length(pi)
	if(!draw.projection[1] || (draw.projection[1] && year.idx[1]==1)) nquant <- 0
	if(nquant > 1) pi<-sort(pi, decr=TRUE)
	quantiles.table <- if(proportion) list(male=pop.pred$quantilesPropMage, female=pop.pred$quantilesPropFage)
                       else list(male=pop.pred$quantilesMage, female=pop.pred$quantilesFage)
    is.valid.pi <- if(proportion && nquant>0) is.saved.pi(pop.pred, pi)
                   else rep(TRUE, nquant)
	maxx<-0
	for(sex in c('male', 'female')) {
		for(iage in 1:lages) {
			if(sum(draw.projection)>0) {
				dimt <- dim(quantiles.table[[sex]])
				dimn <- dimnames(quantiles.table[[sex]])
				table <- drop(quantiles.table[[sex]][,ages.idx[iage],,])
				table <- array(table, dimt[c(1,3:4)], dimnames=c(list(NULL), dimn[3], dimn[4]))
				med <- get.pop.traj.quantiles(table, pop.pred, country$index, country$code, 
												q=0.5, sex=sex, age=ages.idx[iage])
			}
			if(sum(draw.projection) < 2) observed.data <- get.pop.observed(pop.pred, country$code, sex=sex, age=iage)
			pop.median[[sex]][iage] <- if(draw.projection[1]) med[year.idx[1]] 
										else observed.data[year.idx[1]]/(if(proportion) pop.observed[year.idx[1]] else 1)
			if(is.na(pop.median[[sex]][iage])) pop.median[[sex]][iage] <- 0
			maxx <- max(maxx, pop.median[[sex]][iage])
			if(draw.observed) {
				pop.present[[sex]][iage] <- if(draw.projection[2]) med[year.idx[2]] 
											else observed.data[year.idx[2]]/(if(proportion) pop.observed[year.idx[2]] else 1)
				if(is.na(pop.present[[sex]][iage])) pop.present[[sex]][iage] <- 0
				maxx <- max(maxx, pop.present[[sex]][iage])
			}
			if (nquant == 0) next
			for (i in 1:nquant) {
				if (!is.valid.pi[i]) {
					pop.quant[[sex]][[i]] <- NA
					next
				}
				quant <- get.pop.traj.quantiles(table, 
												pop.pred, country$index, country$code, 
												pi=pi[i], sex=sex, age=ages.idx[iage])
				if(length(pop.quant[[sex]]) < i) 
					pop.quant[[sex]][[i]] <- array(NA, c(2, lages))
				pop.quant[[sex]][[i]][,iage] <- quant[,year.idx[1]]
				maxx <- max(maxx, pop.quant[[sex]][[i]][1:2,iage])
			}
		}
	}

	age.labels <- get.age.labels(pop.pred$ages[ages.idx])
	return(list(country=country, draw.projection=draw.projection, draw.observed=draw.observed, nquant=nquant, 
		pop.median=pop.median, pop.quant=pop.quant, pop.present=pop.present, maxx=maxx, is.valid.pi=is.valid.pi, 
		age.labels=age.labels, lages=lages, ages.idx=ages.idx, year.idx=year.idx, pop.observed=pop.observed))
}

pop.pyramid <- function(pop.pred, country, year=NULL, pi=c(80, 95), proportion=FALSE,
							main=NULL, age=1:21, draw.past.year=FALSE) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	data <- .get.data.for.pyramid(pop.pred, country, year=year, pi=pi, proportion=proportion, age=age, draw.past.year=draw.past.year)
	mgp <- par('mgp')
	oma <- par('oma')
	mar <- par('mar')
	par(mfrow=c(1,2),  oma = c(0, 0, 2, 0))
	par(mar=c(5,6,2,-0.1)+0.1, mgp=c(3,0.5,0))
	present.year.col <- "violet"
	with(data, {
	plot(c(-1,0), c(0, lages), type='n', axes=FALSE, xlab = "", ylab = "", main='Male', first.panel=grid(),
			cex.main=0.9)
	cols <- rainbow(max(nquant, 5), start=0.15)[1:nquant]
	if(nquant > 0) {
		for(i in 1:nquant) {
			if(!is.valid.pi[i]) next
			rect(-pop.quant[['male']][[i]][2,]/maxx, (1:lages)-0.45, 
				-pop.quant[['male']][[i]][1,]/maxx, (1:lages)+0.45, col=cols[i],
				border= NA)
		}
	}
	if(draw.observed) {
		rect(-pop.present[['male']]/maxx, (1:lages)-0.15, rep(0, lages), (1:lages)+0.15, #lwd=2,
			col=present.year.col, border=present.year.col, density=20)
		for (i in 1:lages)
			lines(rep(-pop.present[['male']][i]/maxx, 2), c(i-0.15, i+0.15), lwd=3, col=present.year.col)
	}
	rect(-pop.median[['male']]/maxx, (1:lages)-0.45, rep(0, lages), (1:lages)+0.45, #lwd=2,
			border='black')
	for (i in 1:lages)
		lines(rep(-pop.median[['male']][i]/maxx, 2), c(i-0.45, i+0.45), lwd=3)

	labels <- if(proportion) round(seq(maxx, 0, length=11),2) else round(signif(seq(maxx, 0, length=11),2))
	axis(1, at=-labels/maxx, labels=labels)
	axis(2, at=1:lages, labels=age.labels, las=2)
	par(mar=c(5,-0.1,2,6)+0.1)
	plot(c(0,1), c(0, lages), type='n', axes=FALSE, xlab = "", ylab = "", main='Female', first.panel=grid(),
			cex.main=0.9)
	if(nquant > 0) {
		for(i in 1:nquant) {
			if(!is.valid.pi[i]) next
			rect(pop.quant[['female']][[i]][1,]/maxx, (1:lages)-0.45, 
				pop.quant[['female']][[i]][2,]/maxx, 
				 (1:lages)+0.45, col=cols[i], border=NA)
		}
	}
	if(draw.observed) {
		rect(rep(0, lages), (1:lages)-0.15, pop.present[['female']]/maxx, (1:lages)+0.15, #lwd=2,
			col=present.year.col, border=present.year.col, density=20)
		for (i in 1:lages)
			lines(rep(pop.present[['female']][i]/maxx, 2), c(i-0.15, i+0.15), lwd=3, col=present.year.col)
	}
	rect(rep(0, lages), (1:lages)-0.45, pop.median[['female']]/maxx, (1:lages)+0.45, #lwd=2,
			border='black')
	for (i in 1:lages)
		lines(rep(pop.median[['female']][i]/maxx, 2), c(i-0.45, i+0.45), lwd=3)
	
	labels <- if(proportion) round(seq(0, maxx, length=11),2) else round(signif(seq(0, maxx, length=11),2))
	axis(1, at=labels/maxx, labels=labels)
	axis(4, at=1:lages, labels=age.labels, las=2)
	legend <- if(draw.projection[1] && year.idx[1] > 1) 'median' else 'observed'
    if (any(is.valid.pi)) legend <- c(legend, paste(pi[is.valid.pi], '% PI', sep=''))
    cols <- cols[is.valid.pi]
    lwd <- c(3, rep(5, sum(is.valid.pi)))
    if(draw.observed) {
    	legend <- c(legend, paste((if(draw.projection[2]) pop.pred$proj.years[year.idx[2]] else as.integer(names(pop.observed)[year.idx[2]])) + c(-3, 2), collapse='-'))
    	cols <- c(cols, present.year.col)
    	lwd <- c(lwd, 3)
    }
	legend('topright', legend=legend, bty='n', col=c('black', cols), lwd=lwd)
	if(is.null(main)) 
		main <- paste(country$name, ': ', 
			paste((if(draw.projection[1]) pop.pred$proj.years[year.idx[1]] else as.integer(names(pop.observed)[year.idx[1]])) + c(-3, 2), collapse='-'), sep='')
	mtext(main, line = 0.5, outer = TRUE)
	})	
	par(mgp=mgp, oma=oma, mar=mar)
}

pop.trajectories.pyramid.all <- function(pop.pred, year=NULL,
									output.dir=file.path(getwd(), 'pop.traj.pyramid'),
									output.type="png", verbose=FALSE, ...) {
	# plots pyramid for all countries and all years given by 'year'
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.countries <- pop.pred$countries[,'name']
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	if(is.null(year)) year <- pop.pred$present.year
	for (country in all.countries) {
		country.obj <- get.country.object(country, country.table=pop.pred$countries)
		if(verbose)
			cat('Creating trajectory pyramid(s) for', country, '(', country.obj$code, ')\n')

		for(y in year) {
			do.call(output.type, list(file.path(output.dir, 
										paste('pyr', y, '_c', country.obj$code, '.', postfix, sep=''))))
			pop.trajectories.pyramid(pop.pred, country=country.obj$code, year=y, ...)
			dev.off()
		}
	}
	if(verbose)
		cat('\nTrajectory pyramids stored into', output.dir, '\n')
}


pop.trajectories.pyramid <- function(pop.pred, country, year=NULL, pi=c(80, 95), 
					nr.traj=NULL, proportion=FALSE, main=NULL, age=1:21, draw.past.year=FALSE) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	data <- .get.data.for.pyramid(pop.pred, country, year=year, pi=pi, proportion=proportion, age=age, draw.past.year=draw.past.year)
	mgp <- par('mgp')
	oma <- par('oma')
	mar <- par('mar')
	par(mfrow=c(1,2),  oma = c(0, 0, 2, 0))
	par(mar=c(5,6,2,-0.1)+0.1, mgp=c(3,0.5,0))
	male.trajectories <- female.trajectories <- NULL
	with(data, {
		if(draw.projection[1] && year.idx[1] > 1) {
			mtraj <- get.pop.trajectories.multiple.age(pop.pred, country$code, sex='male', 
										age=age, nr.traj, proportion=proportion)
			if(!is.null(dim(mtraj$trajectories))) {
				male.trajectories <- drop(mtraj$trajectories[,year.idx[1],mtraj$index])
				male.trajectories <- array(male.trajectories, c(dim(mtraj$trajectories)[1],length(mtraj$index)))
				ftraj <- get.pop.trajectories.multiple.age(pop.pred, country$code, sex='female', 
										age=ages.idx, nr.traj, proportion=proportion)
				female.trajectories <- drop(ftraj$trajectories[,year.idx[1],ftraj$index])
				female.trajectories <- array(female.trajectories, c(dim(ftraj$trajectories)[1],length(ftraj$index)))
			}
			if(!is.null(male.trajectories)) maxx <- max(maxx, male.trajectories, female.trajectories)
		}
		plot(c(-1,0), c(0, lages), type='n', axes=FALSE, xlab = "", ylab = "", main='Male', first.panel=grid(),
			cex.main=0.9)
		if(!is.null(male.trajectories)) {
			for(i in 1:dim(male.trajectories)[2]) lines(-male.trajectories[,i]/maxx, 1:lages, col='grey')
		}
		lines(-pop.median[['male']]/maxx, 1:lages, col='red')
		lty <- 2:(nquant+1)
		if(nquant > 0) {
			for(i in 1:nquant) {
				if(!is.valid.pi[i]) next
				lines(-pop.quant[['male']][[i]][1,]/maxx, 1:lages, col='red', lty=lty[i])
				lines(-pop.quant[['male']][[i]][2,]/maxx, 1:lages, col='red', lty=lty[i])
			}
		}
		if(draw.observed) lines(-pop.present[['male']]/maxx, 1:lages, col='black')

		labels <- if(proportion) round(seq(maxx, 0, length=11),2) else round(signif(seq(maxx, 0, length=11),2))
		axis(1, at=-labels/maxx, labels=labels)
		axis(2, at=1:lages, labels=age.labels, las=2)
		abline(v=0)
		
		par(mar=c(5,-0.1,2,6)+0.1)
		plot(c(0,1), c(0, lages), type='n', axes=FALSE, xlab = "", ylab = "", main='Female', first.panel=grid(),
			cex.main=0.9)
		if(!is.null(female.trajectories)) {
			for(i in 1:dim(female.trajectories)[2]) lines(female.trajectories[,i]/maxx, 1:lages, col='grey')
		}
		lines(pop.median[['female']]/maxx, 1:lages, col='red')
		if(nquant > 0) {
			for(i in 1:nquant) {
				if(!is.valid.pi[i]) next
				lines(pop.quant[['female']][[i]][1,]/maxx, 1:lages, col='red', lty=lty[i])
				lines(pop.quant[['female']][[i]][2,]/maxx, 1:lages, col='red', lty=lty[i])
			}
		}
		if(draw.observed) lines(pop.present[['female']]/maxx, 1:lages, col='black')
		labels <- if(proportion) round(seq(0, maxx, length=11),2) else round(signif(seq(0, maxx, length=11),2))
		axis(1, at=labels/maxx, labels=labels)
		axis(4, at=1:lages, labels=age.labels, las=2)
		abline(v=0)
		if(is.null(main)) 
			main <- paste(country$name, ': ', 
			paste((if(draw.projection[1]) pop.pred$proj.years[year.idx[1]] else as.integer(names(pop.observed)[year.idx[1]])) + c(-3, 2), collapse='-'), sep='')
		mtext(main, line = 0.5, outer = TRUE)
		legend <- if(draw.projection[1] && year.idx[1] > 1) 'median' else 'observed'
    	if (any(is.valid.pi)) legend <- c(legend, paste(pi[is.valid.pi], '% PI', sep=''))
    	cols <- rep('red', sum(is.valid.pi)+1)
    	lty <- c(1,lty[is.valid.pi])
    	if(draw.observed) {
    		legend <- c(legend, paste((if(draw.projection[2]) pop.pred$proj.years[year.idx[2]] else as.integer(names(pop.observed)[year.idx[2]])) + c(-3, 2), collapse='-'))
    		cols <- c(cols, 'black')
    		lty <- c(lty, 1)
    	}
		legend('topright', legend=legend, lty=lty, bty='n', col=cols)
	})
	par(mgp=mgp, oma=oma, mar=mar)
}

.get.year.index <- function(year, years) {
	lyears <- length(years)
	res <- as.integer(cut(year, labels=1:lyears, breaks=c(years-3, years[lyears]+2)))
	return(res)
	
	#breaks <- c(years-3, years[lyears]+2)
	#h <- try(hist(year, breaks=breaks, plot=FALSE)$count, silent=TRUE)
	#return(if(inherits(h, "try-error")) NULL else which(h > 0)[1])
}
get.prediction.year.index <- function(pred, year) {
	years <- pred$proj.years
	return(.get.year.index(year, years))
}

get.observed.year.index <- function(pred, year) {
	years <- as.integer(colnames(pred$inputs$pop.matrix[['male']]))
	return(.get.year.index(year, years))
}