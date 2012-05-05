pop.aggregate <- function(pop.pred, regions, method=c('independence', 'regional'),
						name = method,
						inputs=list(e0F.sim.dir=NULL, e0M.sim.dir='joint_', tfr.sim.dir=NULL), 
						verbose=FALSE) {
	data(LOCATIONS)
	regions <- unique(regions)
	method <- match.arg(method)
	if(missing(name)) name <- method
	if(method == 'independence') 
		aggr.pred <- pop.aggregate.independence(pop.pred, regions, name, verbose=verbose)
	if(method == 'regional')
		aggr.pred <- pop.aggregate.regional(pop.pred, regions, name, inputs=inputs, verbose=verbose)
	invisible(aggr.pred)
}

get.countries.for.region <- function(region, pop.pred) {
	reg.idx <- which(LOCATIONS[,'country_code'] == region)
	all.countries <- LOCATIONS[LOCATIONS[,'location_type'] == 4,]
	if(length(reg.idx) <= 0) {
		warning('No region code ', region, ' available.')
		return(c())
	}
	location.type <- LOCATIONS[reg.idx,'location_type']
	if(!is.element(location.type, c(0,2,3))) {
		warning('Invalid location type for region ', region,'. Allowed types: 0,2,3. But is ', location.type)
		return(c())
	}
	if(location.type == 0)  # corresponds to world, i.e. all countries
		countries <- all.countries[,'country_code']
	else {
		what.code <- if(location.type == 2) 'area_code' else 'reg_code'
		countries <- all.countries[is.element(all.countries[,what.code], region),'country_code']
	}
	return(countries[is.element(countries, pop.pred$countries[,'code'])])
}

pop.aggregate.regional <- function(pop.pred, regions, name,
						inputs=list(e0F.sim.dir=NULL, e0M.sim.dir='joint_', tfr.sim.dir=NULL), 
						verbose=FALSE) {
	inp <- pop.pred$inputs
	for (item in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MIGm', 'MIGf', 'SRB', 'PASFR', 'MIGtype'))
		inp[[item]] <- NULL
	aggregated.countries <- list()
	if(verbose) cat('\nAggregating inputs using regional method.')
	for(id in regions) {
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
		inp$PASFR <- rbind(inp$PASFR, .aggregate.pasfr(pop.pred, countries, countries.index, id))
		inp$MIGtype <- rbind(inp$MIGtype, .aggregate.migtype(pop.pred, countries, countries.index, id))
		aggregated.countries[[as.character(id)]] <- countries
	}
	dir <- if(is.null(inputs$tfr.sim.dir)) pop.pred$inputs$TFRpred$mcmc.set$meta$parent.meta$output.dir else inputs$tfr.sim.dir
	if(!is.null(dir)) inp$TFRpred <- get.tfr.prediction(dir)
	dir <- if(is.null(inputs$e0F.sim.dir)) pop.pred$inputs$e0Fpred$mcmc.set$meta$parent.meta$output.dir else inputs$e0F.sim.dir
	if(!is.null(dir)) inp$e0Fpred <- get.e0.prediction(dir)
	if(is.null(inputs$e0M.sim.dir) || inputs$e0M.sim.dir == 'joint_') {
		if(has.e0.jmale.prediction(inp$e0Fpred)) inp$e0Mpred <- get.e0.jmale.prediction(inp$e0Fpred)
	} else inp$e0Mpred <- get.e0.prediction(inputs$e0M.sim.dir)
	outdir <- gsub('predictions', paste('aggregations', name, sep='_'), pop.pred$output.directory)
	aggr.pred <- do.pop.predict(regions, inp=inp, outdir=outdir, nr.traj=pop.pred$nr.traj, ages=pop.pred$ages, verbose=verbose)
	aggr.pred$aggregation.method <- 'regional'
	aggr.pred$aggregated.countries <- aggregated.countries
	return(aggr.pred)
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
		pop.colidx <- which(is.element(as.integer(colnames(popmatrix))-3, unlist(strsplit(colnames(pop.pred$inputs[[item]])[3:ncol(pop.pred$inputs[[item]])], '-'))))
		res[[item]] <- as.data.frame(matrix(NA, nrow=22, ncol=ncol(pop.pred$inputs[[item]]), 
						dimnames=list(c(0, 1, seq(5, 95, by=5), '100+'), colnames(pop.pred$inputs[[item]]))))
		for(age in rownames(res[[item]])) {
			mort.age.idx <- mort.idx & gsub(' ', '', pop.pred$inputs[[item]][,'age']) == age
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

pop.aggregate.independence <- function(pop.pred, regions, name, verbose=verbose) {
	if(verbose) cat('\nAggregating using independence method.')
	nreg <- length(regions)
	quantiles.to.keep <- as.numeric(dimnames(pop.pred$quantiles)[[2]])
	quant <- quantM <- quantF <- array(NA, c(nreg, dim(pop.pred$quantiles)[2:3]), dimnames=dimnames(pop.pred$quantiles))
	quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(nreg, dim(pop.pred$quantilesMage)[2:4]),
						dimnames=dimnames(pop.pred$quantilesMage))
	mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(nreg,dim(pop.pred$traj.mean.sd)[2:3]))
	outdir <- gsub('predictions', paste('aggregations', name, sep='_'), pop.pred$output.directory)
	if(file.exists(outdir)) unlink(outdir, recursive=TRUE)
	dir.create(outdir, recursive=TRUE)
	aggregated.countries <- list()
	id.idx <- 0
	valid.regions <- rep(FALSE, length(regions))
	for(reg.idx in 1:length(regions)) {
		id <- regions[reg.idx]
		if(verbose) cat('\nAggregating region ', id)
		countries <- get.countries.for.region(id, pop.pred)
		if(length(countries)==0) next
		id.idx <- id.idx + 1
		valid.regions[reg.idx] <- TRUE
		countries.index <- which(is.element(pop.pred$countries[,'code'], countries))
		for(cidx in 1:length(countries.index)) {
			traj.file <- file.path(pop.pred$output.directory, paste('totpop_country', countries[cidx], '.rda', sep=''))
			load(traj.file)
			if(cidx == 1) {
				stotp <- totp
				stotpm <- totpm
				stotpf <- totpf
				stotp.hch <- totp.hch
				stotpm.hch <- totpm.hch
				stotpf.hch <- totpf.hch
			} else {
				stotp <- stotp + totp
				stotpm <- stotpm + totpm
				stotpf <- stotpf + totpf
				stotp.hch <- stotp.hch + totp.hch
				stotpm.hch <- stotpm.hch + totpm.hch
				stotpf.hch <- stotpf.hch + totpf.hch
			}
		}
		totp <- stotp
		totpm <- stotpm
		totpf <- stotpf
		totp.hch <- stotp.hch
		totpm.hch <- stotpm.hch
		totpf.hch <- stotpf.hch
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch,
			 file = file.path(outdir, paste('totpop_country', id, '.rda', sep='')))
		quant[id.idx,,] = apply(stotp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sd[id.idx,1,] <- apply(stotp, 1, mean, na.rm = TRUE)
		mean_sd[id.idx,2,] = apply(stotp, 1, sd, na.rm = TRUE)
		for (i in 1:dim(stotpm)[1]) {
			quantMage[id.idx,i,,] <- apply(stotpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
			quantFage[id.idx,i,,] = apply(stotpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
			quantPropMage[id.idx,i,,] <- apply(stotpm[i,,]/stotp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			quantPropFage[id.idx,i,,] <- apply(stotpf[i,,]/stotp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		}
		sstotpm <- colSums(stotpm)
		quantM[id.idx,,] = apply(sstotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdM[id.idx,1,] <- apply(sstotpm, 1, mean, na.rm = TRUE)
		mean_sdM[id.idx,2,] = apply(sstotpm, 1, sd, na.rm = TRUE)
		sstotpf <- colSums(stotpf)
		quantF[id.idx,,] = apply(sstotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdF[id.idx,1,] <- apply(sstotpf, 1, mean, na.rm = TRUE)
		mean_sdF[id.idx,2,] = apply(sstotpf, 1, sd, na.rm = TRUE)
		aggregated.countries[[as.character(id)]] <- countries
	}
	aggr.pred <- pop.pred
	which.reg.index <- function(x, set) return(which(set == x))
	reg.idx <- sapply(regions[valid.regions], which.reg.index, set=LOCATIONS[,'country_code']) 
	aggr.pred$countries=data.frame(code=LOCATIONS[reg.idx, 'country_code'], name=LOCATIONS[reg.idx, 'name'])
	aggr.pred$output.directory <- outdir
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
	aggr.pred$aggregation.method <- 'independence'
	aggr.pred$aggregated.countries <- aggregated.countries
	bayesPop.prediction <- aggr.pred
	save(bayesPop.prediction, file=file.path(outdir, 'prediction.rda'))
	cat('\nAggregations stored into', outdir, '\n')
	return(bayesPop.prediction)
}