adjust.trajectories <- function(country, env, quant.env, adj.env=NULL, allow.negatives = TRUE,
                                adj.to.file = NULL) {
	if(is.null(adj.env)) adj.env <- new.env()
	wpp.year <- quant.env$wpp.year
	annual <- quant.env$annual
	if(is.null(annual)) annual <- FALSE
	datasets <- list(totp='', totp.hch='', totpf='Fage', totpm='Mage', totpm.hch='Mage', totpf.hch='Fage')
	country.char <- as.character(country)
	for(traj.name in names(datasets)) {
		adj.name <- datasets[[traj.name]]
		dif.name <- paste0('AdjDpop', adj.name)
		
		if(is.null(adj.env[[dif.name]])) {
			#print(c(dif.name, adj.name))
			q <- quant.env[[paste0('quantiles', adj.name)]]
			adjust.quantiles(q, adj.name, wpp.year=wpp.year, env=adj.env, annual = annual, adj.to.file = adj.to.file)
		}
		dif <- if(length(dim(adj.env[[dif.name]]))>2) adj.env[[dif.name]][country.char,,] else adj.env[[dif.name]][country.char,,drop=FALSE]
		res <- env[[traj.name]]		
		if(length(dim(res))>2) { # includes age
		    age.idx <- age.index.all(annual = annual, observed = TRUE)
			res21 <- aaply(res[age.idx,,], 3, '-', dif)
            if(!allow.negatives) res21 <- pmax(0, res21)
			res21 <- aperm(res21, c(2,3,1))
			res[age.idx,,] <- res21
		} else {
			res <- aaply(res, 2, '-', dif)
			if(!allow.negatives) res <- pmax(0, res)
			res <- aperm(res, c(2,1))
		}
		env[[traj.name]] <- res
	}
}

adjust.quantiles <- function(q, what, wpp.year, annual = FALSE, env=NULL, allow.negatives = TRUE,
                             adj.to.file = NULL) {
	dif <- NULL
	if(!is.null(env)) {
		if(!is.null(env[[paste0('AdjQpop', what)]])) return(env[[paste0('AdjQpop', what)]])
		if(!is.null(env[[paste0('AdjDpop', what)]])) dif <- env[[paste0('AdjDpop', what)]]
	}
	age.idx <- age.index.all(annual = annual, observed = TRUE)
	age.idx.all <- age.index.all(annual = annual, observed = FALSE)
	age.idx.old <- age.idx.all[age.idx.all > max(age.idx)]
	if(is.null(dif)) {
		if(is.null(env)) env <- new.env()
		countries <- dimnames(q)[[1]]
		ages <- NULL
		if(length(dim(q))>3) { # includes age dimension
			ages <- dimnames(q)[[2]]
			ages <- ages[as.numeric(ages)<=100]
		}
		wpp <- if(!is.null(adj.to.file)) .get.adjustments.from.file(
		    adj.to.file, env, what, countries, ages, wpp.year=wpp.year, annual = annual) else .get.wpp(
		        env, what, countries, ages, wpp.year=wpp.year, annual = annual)
		if(length(dim(q))>3) { # includes age dimension
			years <- as.numeric(dimnames(q)[[4]])
			if(!annual && (years[1] %% 5) != 0) years <- years+2 
			med.raw <- q[,,'0.5',as.character(years)%in%dimnames(wpp)[[3]]]
			if(length(dim(med.raw))>2) { # multiple countries
				med <- med.raw[,age.idx,] # collapse to 21 age categories
				med[,age.idx.old[1]-1,] <- med.raw[,age.idx.old[1]-1,] + apply(med.raw[,age.idx.old,], c(1,3), sum) 
			} else { #1 country
				med <- med.raw[age.idx,]
				med[age.idx.old[1]-1,] <- med.raw[age.idx.old[1]-1,] + apply(med.raw[age.idx.old,], 2, sum) 
				med <- abind(med, along=0) # add dimension
			}
			dif <- med-wpp
			if(! years[1] %in% dimnames(dif)[[3]])
			    dif <- abind(matrix(0, nrow=dim(med)[1], ncol=length(age.idx)), dif, along=3)		
		} else {
			years <- as.numeric(dimnames(q)[[3]])
			if(!annual && (years[1] %% 5) != 0) years <- years+2 
			med <- q[,'0.5',as.character(years)%in%colnames(wpp)]
			dif <- med-wpp
			if(! years[1] %in% colnames(dif))
			    dif <- cbind(0, dif) # add a column for present year
			dif <- as.matrix(dif)
		}
	} else countries <- dimnames(dif)
	if(length(dim(q))>3) {
		if(dim(q)[1]==1){ # one country - the generic aaply fails because of some dimension dropping 
			res21 <- aaply(q[,age.idx,,], 2, '-', dif[1,,], .drop=FALSE)
			res21 <- aperm(res21, c(2,1,3))
		} else {
			res21 <- aaply(q[,age.idx,,], 3, '-', dif, .drop=FALSE)
			res21 <- aperm(res21, c(2,3,1,4))
		}
		res <- q
		res[,age.idx,,] <- res21
	} else { # no age dimension
		res <- aaply(q, 2, '-', dif, .drop=FALSE)
		res <- aperm(res, c(2,1,3))
	}
	if(is.null(dimnames(dif)[[1]])) dimnames(dif)[[1]] <- countries
	if(is.null(dimnames(res)[[1]])) dimnames(res)[[1]] <- countries
	if(!allow.negatives) res <- pmax(0, res)
	env[[paste0('AdjDpop', what)]] <- dif
	env[[paste0('AdjQpop', what)]] <- res
	return(res)
}

.get.wpp <- function(env, what, countries=NULL, ages=NULL, annual = FALSE, ...) {
	switch(which(c('', 'M', 'F', 'Mage', 'Fage') == what), 
				tpop(countries, prediction.only=TRUE, e=env, annual = annual, ...),
				tpopM(countries, prediction.only=TRUE, e=env, annual = annual, ...),
				tpopF(countries, prediction.only=TRUE, e=env, annual = annual, ...),
				tpopM(countries, prediction.only=TRUE, sum.over.ages=FALSE, ages=ages, e=env, annual = annual, ...),
				tpopF(countries, prediction.only=TRUE, sum.over.ages=FALSE, ages=ages, e=env, annual = annual, ...)
			)
}

.get.adjustments.from.file <- function(file, env, what, countries=NULL, ages=NULL, annual = FALSE, ...) {
    adj.dataset <- fread(file)
    year.cols <- grep('^[0-9]{4}', colnames(adj.dataset), value = TRUE)
    if("reg_code" %in% colnames(adj.dataset)) adj.dataset[, country_code := reg_code]
    idx <- if(!is.null(countries)) which(adj.dataset$country_code %in% countries) else 1:nrow(adj.dataset)
    adj.dataset.long <- data.table::melt(adj.dataset[idx, c("country_code", "sex", "age", year.cols), with = FALSE], 
                                         id.vars = c("country_code", "sex", "age"), variable.name = "year",
                                         variable.factor = FALSE)
    if(what == "") sx <- c("male", "female")
    else sx <- list(M = "male", F = "female", Mage = "male", Fage = "female")[[what]]
    adj.dataset.long <- adj.dataset.long[tolower(sex) %in% sx][, sex := NULL]
    
    if(!what %in% c("Mage", "Fage")) {
        adj.data <- adj.dataset.long[, list(value = sum(value)), by = c("country_code", "year")]
        adj.data.wide <- data.table::dcast(adj.data, country_code ~ year, value.var = "value")
        adj.data.mat <- as.matrix(adj.data.wide[, -1, with = FALSE])
        rownames(adj.data.mat) <- adj.data.wide$country_code
        return(adj.data.mat)
    }
    # age-specific 
    adj.data.wide <- data.table::dcast(adj.dataset.long, country_code + age ~ year, value.var = "value")
    return(.reduce.to.countries.and.ages(adj.data.wide, countries, ages, annual = annual))
}

if.not.exists.load <- function(name, env, wpp.year=2012) {
	if(!exists(name, where=env, inherits=FALSE)) {
		do.call('data', list(name, package=paste0('wpp', wpp.year), envir=env))
	    env[[name]] <- as.data.table(env[[name]])
	}
}

tpop <- function(countries, prediction.only=FALSE, e=NULL, annual = FALSE, ...) {
	# Create a dataset of total population
	if(is.null(e)) e <- new.env()
	suffix <- if(annual) "1" else ""
	if(!prediction.only) {
		if.not.exists.load(paste0('popM', suffix), e, ...)
		if.not.exists.load(paste0('popF', suffix), e, ...)
		tpop.obs <- sumMFbycountry(paste0('popM', suffix), paste0('popF', suffix), e)
	}
	#projection stored separately from observations
	if.not.exists.load(paste0('popMprojMed', suffix), e, ...)
	if.not.exists.load(paste0('popFprojMed', suffix), e, ...)
	tpopp <- sumMFbycountry(paste0('popMprojMed', suffix), paste0('popFprojMed', suffix), e)
	if(!prediction.only) tpopp <- merge(tpop.obs, tpopp, by='country_code')
	return(.reduce.to.countries(tpopp, countries))
}

tpopF <- function(...) return(tpop.sex('F', ...))
tpopM <- function(...) return(tpop.sex('M', ...))

tpop.sex <- function(sex, countries, sum.over.ages=TRUE, ages=NULL, prediction.only=FALSE, e=NULL, annual = FALSE, ...) {
	# Create a dataset of total population by sex
	if(is.null(e)) e <- new.env()
	suffix <- if(annual) "1" else ""
	if(!prediction.only) {
		dataset <- paste0('pop', sex, suffix)
		if.not.exists.load(dataset, e, ...)
		#do.call('data', list(dataset, package='wpp2012', envir=e))
		pop.obs <- if(sum.over.ages) .sum.by.country(dataset) else .sum.by.country.and.age(dataset)
	}
	dataset <- paste0('pop', sex, 'projMed', suffix)
	if.not.exists.load(dataset, e, ...)
	popp <- if(sum.over.ages) .sum.by.country(e[[dataset]]) else .sum.by.country.and.age(e[[dataset]])
	if(!prediction.only)  popp <- merge(pop.obs, popp, by='country_code')
	if(sum.over.ages) return(.reduce.to.countries(popp, countries))
	.reduce.to.countries.and.ages(popp, countries, ages, annual = annual)
}

.reduce.to.countries <- function(dataset, countries){
    tpop <- as.data.frame(dataset)
	tpop <- tpop[,-which(colnames(dataset)=='country_code')]
	rownames(tpop) <- dataset$country_code
	tpop[countries,]
}

.reduce.to.countries.and.ages <- function(dataset, countries, ages, annual = FALSE){
	dataset <- as.data.frame(dataset[dataset$country_code %in% as.integer(countries),])
	if(annual){
	    by <- 1
	} else {
	    by <- 5
	}
	nage <- length(unique(dataset$age))
	if(is.null(ages)) ages <- as.character(seq(0,100, by=by))
	age.vector <- as.character(dataset$age[1:nage])
	if(!annual) {
	    age.vector <- unlist(strsplit(gsub('\\+', '-130', age.vector), '-'))
	    age.vector <- age.vector[seq(1,length(age.vector), by=2)]
	} else age.vector <- gsub('\\+', '', age.vector)
	colidx <- (1:ncol(dataset))[-which(colnames(dataset) %in% c('country_code', 'age'))]
	res <- array(NA, c(length(countries), length(ages), ncol(dataset)-2))
	for(i in 1:length(countries)) {
		idx <- which(dataset$country_code==countries[i])
		if(length(idx) == 0) next
		tmp <- dataset[idx,colidx]
		rownames(tmp) <- age.vector
		res[i,,] <- as.matrix(tmp[ages,])
	}
	dimnames(res)<- list(countries, ages, colnames(dataset)[colidx])
	res
}

.sum.by.country <- function(dataset) {
	year.cols <- grep('^[0-9]{4}', colnames(dataset), value = TRUE)
	dataset[, c("country_code", year.cols), with = FALSE][, lapply(.SD, sum, na.rm = TRUE), by = "country_code"]
}

.sum.by.country.and.age <- function(dataset) {
	year.cols <- grep('^[0-9]{4}', colnames(dataset), value = TRUE)
	dataset[, c("country_code", "age", year.cols), with = FALSE][, lapply(.SD, sum, na.rm = TRUE), by = c("country_code", "age")]
}

sumMFbycountry <- function(datasetM, datasetF, e) {
	tpopM <- .sum.by.country(e[[datasetM]])
	tpopF <- .sum.by.country(e[[datasetF]])
	tpopM[, 2:ncol(tpopM)] <- tpopM[,2:ncol(tpopM)] + tpopF[,2:ncol(tpopF)]
    tpopM
}

adjust.to.dataset <- function(country, q, adj.dataset=NULL, adj.file=NULL, years=NULL, 
                              use=c('write', 'trajectories'), allow.negatives = TRUE) {
	if(is.null(adj.dataset)) {
		adj.dataset <- read.table(adj.file, header=TRUE, check.names=FALSE)
	}
	colidx <- if(is.null(years)) (1:ncol(adj.dataset))[-which(colnames(adj.dataset)%in%c('country_code', 'country', 'name'))] else as.character(years)
	idx1 <- which(adj.dataset$country_code == country)
	if(use=='write') {
		med <- q['0.5']
		dif <- med - adj.dataset[idx1,colidx]
		if(!allow.negatives) dif <- pmin(q, dif)
		return(q-dif)
	}
	if(use=='trajectories') {
		med <- apply(q, 1, 'median', na.rm = TRUE)[colnames(adj.dataset[,colidx])]
		dif <- as.matrix(med - adj.dataset[idx1,colidx])
		res <- aaply(q[colnames(adj.dataset[,colidx]),], 2, '-', dif)
		if(!allow.negatives) res <- pmax(res, 0)
		res <- aperm(res, c(2,1))
		if(! rownames(q)[1] %in% rownames(res))
			res <- rbind(q[1,], res) # add current year
		rownames(res) <- rownames(q)
		return(res)
	}
	return(NULL)
}

adjust.to.aggregation <- function(pop.pred, target.file, output.file = "adjusted_population.txt", 
                                        target.code = NULL, variant.name = "median", 
                                        target.id.column = "country_code", output.id.column = "reg_code",
                                        stat = "mean", exclude.codes = NULL){
    
    if(!has.pop.aggregation(pop.pred = pop.pred))
        stop("The pop.pred object does not contain any aggregation. Consider running pop.aggregate().")
    
    if(!stat %in% c("median", "mean"))
        stop("Argument 'stat' must be either 'mean' or 'median'.")
    
    agdat <- data.table::fread(target.file)
    
    # select the desired variant
    if("variant" %in% colnames(agdat)){ 
        if(! tolower(variant.name) %in% tolower(agdat$variant))
            stop("Value '", variant.name, "' not found in the column 'variant' of the target.file. Consider changing the 'variant.name' argument.")
        agdat <- agdat[variant == variant.name][, variant := NULL]
    }
    
    # select the desired aggregated location
    if(target.id.column %in% colnames(agdat)){
        if(!is.null(target.code)){
            agdat <- agdat[get(target.id.column) == target.code][, (target.id.column) := NULL]
            if(nrow(agdat) == 0) stop("Value ", target.code, " not found in target.file.")
        } else {
            if(length((target.code <- unique(agdat[[target.id.column]]))) > 1){
                stop("target.file contains more than one location. Set 'target.code' to specifify which location to use.")
            }
        }
    } else {
        if(!is.null(target.code) && any(duplicated(agdat[, c("sex", "age"), with = FALSE]))) 
            stop("Column ", target.id.column, " not found in target.file. Consider changing the 'target.id.column' argument.")
    }
    
    if(any(duplicated(agdat[, c("sex", "age"), with = FALSE])))
        stop("The target.file contains duplicates of sex and age combinations.")
    
    # remove + from age and convert to numeric
    if(pop.pred$annual) agdat[, age := as.integer(gsub("+", "", age, fixed = TRUE))]
    
    # make it a long format
    year.cols <- grep('^[0-9]{4}', colnames(agdat), value = TRUE)
    agdat.long <- data.table::melt(agdat, id.vars = c("sex", "age"), measure.vars = year.cols,
                                   variable.name = "year", variable.factor = FALSE, 
                                   value.name = "target")
    agdat.long[, year := as.integer(year)]
    
    # get simulated aggregation
    pop.aggr <- get.pop.aggregation(pop.pred = pop.pred)
    
    # determine the id of the (simulated) aggregated location
    agsimlocs <- pop.aggr$countries
    if(nrow(agsimlocs) > 1){ # there is more than 1 aggregation
        if(is.null(target.code))
            stop("The aggregation of pop.pred contain more than one aggregated locations. Consider setting 'target.code'.")
        agsim.id <- target.code
    } else agsim.id <- agsimlocs$code
    
    # get simulated aggregated pop trajectories
    popaggr.traj <- rbind(get.pop.exba(paste0("P", agsim.id, "_M{}"), pop.aggr, as.dt = TRUE)[, sex := "male"],
                          get.pop.exba(paste0("P", agsim.id, "_F{}"), pop.aggr, as.dt = TRUE)[, sex := "female"]
                        )
    
    # derive the desired statistics from simulated aggregation
    popaggr.stat <- popaggr.traj[, list(sim = do.call(stat, list(indicator))), by = c("year", "sex", "age")]
    
    # merge together and determine the differences
    targets <- merge(agdat.long, popaggr.stat, by = c("year", "sex", "age"))[, shift := sim - target]
    
    # get simulated trajectories and statistics for all locations
    pop.traj <- rbind(get.pop.exba(paste0("PXXX_M{}"), pop.pred, as.dt = TRUE)[, sex := "male"],
                      get.pop.exba(paste0("PXXX_F{}"), pop.pred, as.dt = TRUE)[, sex := "female"]
                      )
    pop.stat <- pop.traj[, list(sim = do.call(stat, list(indicator))), by = c("country_code", "year", "sex", "age")]
    
    # compute pop shares of each location within the aggregation,
    # to be used for distributing the adjustments
    pop.stat[, share := sim / sum(sim), by = c("year", "sex", "age")]
    
    # for zero shares, use shares derived from total population
    totpop <- pop.stat[, list(totpop = sum(sim)), by = c("country_code", "year")][, totshare := totpop / sum(totpop), by = c("year")]
    pop.stat[totpop, totshare := i.totshare, on = c("country_code", "year")]
    pop.stat[is.na(share), share := totshare]
    
    # exclude locations
    if(!is.null(exclude.codes)){
        pop.stat[country_id %in% exclude.codes, share := 0]
    }
    
    # rescale
    pop.stat[, share := share / sum(share), by = c("year", "sex", "age")]
    
    # merge in the total shift
    pop.stat[targets, totshift := i.shift, on = c("year", "sex", "age")]
    # if the target is smaller than simulated set the share for zero age groups to 0
    pop.stat[totshift > 0 & sim == 0, share := 0] 

    # rescale again
    pop.stat[, share := share / sum(share), by = c("year", "sex", "age")][is.na(share), share := 0]
    
    
    # compute adjusted (target) pop for each location
    pop.stat[, shift := totshift * share][, simadj := pmax(0, sim - shift)]
    
    respop <- data.table::dcast(pop.stat, country_code + sex + age ~ year, value.var = "simadj")
    data.table::setnames(respop, "country_code", output.id.column)
    
    data.table::fwrite(respop, file = output.file, sep = "\t")
    cat("\nAdjustment file written into", output.file, "\n")
    return(output.file)
}