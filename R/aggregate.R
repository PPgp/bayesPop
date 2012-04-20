pop.aggregate <- function(regions, pop.pred, 
						inputs=list(e0F.sim.dir=NULL, e0M.sim.dir='joint_', tfr.sim.dir=NULL), 
						verbose=FALSE) {
	data(LOCATIONS)
	regions <- unique(regions)
	inp <- pop.pred$inputs
	for (item in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MIGm', 'MIGf', 'SRB', 'PASFR', 'MIGtype'))
		inp[[item]] <- NULL
	if(verbose) cat('\nAggregating inputs.')
	for(id in regions) {
		countries <- LOCATIONS[is.element(LOCATIONS[,'area_code'], id),'country_code']
		countries <- countries[is.element(countries, pop.pred$countries[,'code'])]
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
	}
	if(!is.null(inputs$tfr.sim.dir))
		inp$TFRpred <- get.tfr.prediction(inputs$tfr.sim.dir)
	if(!is.null(inputs$e0F.sim.dir))
		inp$e0Fpred <- get.e0.prediction(inputs$e0F.sim.dir)
	if(is.null(inputs$e0M.sim.dir) || inputs$e0M.sim.dir == 'joint_') 
		if(has.e0.jmale.prediction(inp$e0Fpred)) inp$e0Mpred <- get.e0.jmale.prediction(inp$e0Fpred)

	outdir <- gsub('predictions', 'aggregations', pred$output.dir)
	invisible(do.pop.predict(regions, inp=inp, outdir=outdir, nr.traj=pop.pred$nr.traj, ages=pop.pred$ages, verbose=verbose))
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