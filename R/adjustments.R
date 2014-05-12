adjust.trajectories <- function(country, env, quant.env, adj.env=NULL) {
	if(is.null(adj.env)) adj.env <- new.env()
	datasets <- list(totp='', totp.hch='', totpf='Fage', totpm='Mage', totpm.hch='Mage', totpf.hch='Fage')
	country.char <- as.character(country)
	for(traj.name in names(datasets)) {
		adj.name <- datasets[[traj.name]]
		dif.name <- paste0('AdjDpop', adj.name)
		if(is.null(adj.env[[dif.name]])) {
			#print(c(dif.name, adj.name))
			q <- quant.env[[paste0('quantiles', adj.name)]]
			adjust.quantiles(q, adj.name, adj.env)
		}
		dif <- if(length(dim(adj.env[[dif.name]]))>2) adj.env[[dif.name]][country.char,,] else adj.env[[dif.name]][country.char,,drop=FALSE]

		res <- env[[traj.name]]		
		if(length(dim(res))>2) { # includes age
			res21 <- aaply(res[1:21,,], 3, '-', dif)
			res21 <- aperm(res21, c(2,3,1))
			res[1:21,,] <- res21
		} else {
			res <- aaply(res, 2, '-', dif)
			res <- aperm(res, c(2,1))
		}
		if(is.list(res)) stop('')
		env[[traj.name]] <- res
	}
}

adjust.quantiles <- function(q, what, env=NULL) {
	dif <- NULL
	if(!is.null(env)) {
		if(!is.null(env[[paste0('AdjQpop', what)]])) return(env[[paste0('AdjQpop', what)]])
		if(!is.null(env[[paste0('AdjDpop', what)]])) dif <- env[[paste0('AdjDpop', what)]]
	}
	if(is.null(dif)) {
		if(is.null(env)) env <- new.env()
		countries <- dimnames(q)[[1]]
		ages <- NULL
		if(length(dim(q))>3) { # includes age dimension
			ages <- dimnames(q)[[2]]
			ages <- ages[as.numeric(ages)<=100]
		}
		wpp <- .get.wpp(env, what, countries, ages)
		if(length(dim(q))>3) { # includes age dimension
			years <- as.numeric(dimnames(q)[[4]])
			if((years[1] %% 5) != 0) years <- years+2 
			med.raw <- q[,,'0.5',as.character(years)%in%dimnames(wpp)[[3]]]
			med <- med.raw[,1:21,] # collapse to 21 age categories
			med[,21,] <- med.raw[,21,] + apply(med.raw[,22:27,], c(1,3), sum) 
			dif <- abind(matrix(0, nrow=dim(med)[1], ncol=21), med-wpp, along=3)

		} else {
			years <- as.numeric(dimnames(q)[[3]])
			if((years[1] %% 5) != 0) years <- years+2 
			med <- q[,'0.5',as.character(years)%in%colnames(wpp)]
			dif <- as.matrix(cbind(0, med-wpp))
		}
	}
	if(length(dim(q))>3) {
		res21 <- aaply(q[,1:21,,], 3, '-', dif)
		res21 <- aperm(res21, c(2,3,1,4))
		res <- q
		res[,1:21,,] <- res21
	} else {
		res <- aaply(q, 2, '-', dif)
		res <- aperm(res, c(2,1,3))
	}
	env[[paste0('AdjDpop', what)]] <- dif
	env[[paste0('AdjQpop', what)]] <- res
	return(res)
}

.get.wpp <- function(env, what, countries=NULL, ages=NULL) {
	switch(which(c('', 'M', 'F', 'Mage', 'Fage') == what), 
				tpop(countries, prediction.only=TRUE, e=env),
				tpopM(countries, prediction.only=TRUE, e=env),
				tpopF(countries, prediction.only=TRUE, e=env),
				tpopM(countries, prediction.only=TRUE, sum.over.ages=FALSE, ages=ages, e=env),
				tpopF(countries, prediction.only=TRUE, sum.over.ages=FALSE, ages=ages, e=env)
			)
}

if.not.exists.load <- function(name, env) {
	if(!exists(name, where=env, inherits=FALSE))
		do.call('data', list(name, package='wpp2012', envir=env))
}

tpop <- function(countries, prediction.only=FALSE, e=NULL) {
	# Create a dataset of total population
	if(is.null(e)) e <- new.env()
	if(!prediction.only) {
		if.not.exists.load('popM', e)
		if.not.exists.load('popF', e)
		tpop.obs <- sumMFbycountry('popM', 'popF', e)
	}
	#projection stored separately from observations
	if.not.exists.load('popMprojMed', e)
	if.not.exists.load('popFprojMed', e)
	tpopp <- sumMFbycountry('popMprojMed', 'popFprojMed', e)
	if(!prediction.only) tpopp <- merge(tpop.obs, tpopp, by='country_code')
	return(.reduce.to.countries(tpopp, countries))
}

tpopF <- function(...) return(tpop.sex('F', ...))
tpopM <- function(...) return(tpop.sex('M', ...))

tpop.sex <- function(sex, countries, sum.over.ages=TRUE, ages=NULL, prediction.only=FALSE, e=NULL) {
	# Create a dataset of total population by sex
	if(is.null(e)) e <- new.env()
	if(!prediction.only) {
		dataset <- paste0('pop', sex)
		do.call('data', list(dataset, package='wpp2012', envir=e))
		pop.obs <- if(sum.over.ages) sum.by.country(dataset) else sum.by.country.and.age(dataset)
	}
	dataset <- paste0('pop', sex, 'projMed')
	if.not.exists.load(dataset, e)
	popp <- if(sum.over.ages) sum.by.country(e[[dataset]]) else sum.by.country.and.age(e[[dataset]])
	if(!prediction.only)  popp <- merge(pop.obs, popp, by='country_code')
	if(sum.over.ages) return(.reduce.to.countries(popp, countries))
	.reduce.to.countries.and.ages(popp, countries, ages)
}

.reduce.to.countries <- function(dataset, countries){
	tpop <- dataset[,-which(colnames(dataset)=='country_code')]
	rownames(tpop) <- dataset$country_code
	tpop <- tpop[countries,]
	tpop
}

.reduce.to.countries.and.ages <- function(dataset, countries, ages){
	dataset <- dataset[dataset$country_code %in% as.integer(countries),]
	if(is.null(ages)) ages <- as.character(seq(0,100, by=5))
	age.vector <- as.character(dataset$age[1:21])
	age.vector <- unlist(strsplit(gsub('\\+', '-130', age.vector), '-'))
	age.vector <- age.vector[seq(1,length(age.vector), by=2)]
	colidx <- (1:ncol(dataset))[-which(colnames(dataset) %in% c('country_code', 'age'))]
	res <- array(NA, c(length(countries), length(ages), ncol(dataset)-2))
	for(i in 1:length(countries)) {
		idx <- which(dataset$country_code==countries[i])
		tmp <- dataset[idx,colidx]
		rownames(tmp) <- age.vector
		res[i,,] <- as.matrix(tmp[ages,])
	}
	dimnames(res)<- list(countries, ages, colnames(dataset)[colidx])
	res
}


sum.by.country <- function(dataset) {
	year.cols.idx <- grep('^[0-9]{4}', colnames(dataset))
	ddply(dataset[,c(which(colnames(dataset)=='country_code'), year.cols.idx)], "country_code", .fun=colwise(sum))
}

sum.by.country.and.age <- function(dataset) {
	year.cols.idx <- grep('^[0-9]{4}', colnames(dataset))
	ddply(dataset[,c(which(colnames(dataset) %in% c('country_code', 'age')), year.cols.idx)], c("country_code", "age"), .fun=colwise(sum))
}

sumMFbycountry <- function(datasetM, datasetF, e) {
	tpopM <- sum.by.country(e[[datasetM]])
	tpopF <- sum.by.country(e[[datasetF]])
	cbind(country_code=tpopM[,'country_code'], tpopM[,2:ncol(tpopM)] + tpopF[,2:ncol(tpopF)])
}

