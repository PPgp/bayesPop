
project.pasfr <- function(inputs = NULL, present.year = 2020, end.year = 2100, 
                          wpp.year = 2019, annual = FALSE, nr.est.points = if(annual) 15 else 3,
                          digits = 2, out.file.name = "percentASFR.txt", verbose = FALSE) {
  # Generate PASFR for given TFR using the kantorova.pasfr method.
  # This function allows to generate PASFR outside of a bayesPop simulation.
  # Results are written into a file given by out.file.name.
  # ########
    bayesTFR:::load.bdem.dataset('UNlocations', wpp.year, envir=globalenv(), verbose = verbose)
    observed <- list()
    
    # determine projection years and periods
    if(annual) {
        proj.periods <- proj.years <- (present.year + 1):end.year
        obs.periods <- obs.years <- (present.year - 300):present.year  # for now take all observed data; it will be filtered later
        period.step <- 1
    } else {
        # find a 5-year bracket for present year
        possible.years <- seq(1000, 3000, by = 5)
        period <- cut(present.year, possible.years, labels = FALSE, right = FALSE)
        # generate vectors of observed and projection years
        proj.years <- seq(possible.years[period]+3, end.year-2, by = 5)
        proj.periods <- paste(proj.years-3, proj.years+2, sep = "-")
        obs.years <- seq(possible.years[period-50]+3, proj.years[1]-5, by = 5)
        obs.periods <- paste(obs.years-3, obs.years+2, sep = "-")
        period.step <- 5
    }
    # get PASFR data
    pasfrlist <- .get.pasfr.data(inputs$pasfr, wpp.year, obs.periods, include.projection = FALSE)
    PASFR <- pasfrlist$pasfr
    observed$PASFR <- pasfrlist$obs.pasfr
    
    # get TFR
    TFRpred <- .get.tfr.data(inputs, wpp.year, verbose = verbose)
    if(!is.data.frame(TFRpred)) { # object of class bayesTFR.prediction
        country.table <- get.countries.table(TFRpred)
    } else {
        TFRpred.sub <- TFRpred[TFRpred$year %in% c(obs.years, proj.years),]
        if(nrow(TFRpred.sub) == 0 && nchar(TFRpred$year[1]==9)) {# year column might be in the form of period instead of year
            TFRpred$year <- as.integer(substr(TFRpred$year, 1,4))+3
            TFRpred.sub <- TFRpred[TFRpred$year %in% c(obs.years, proj.years),]
        }
        TFRpred <- aggregate(value ~ year + country_code, TFRpred.sub, median) # in case there are more trajectories
        #TFRpred <- aggregate(formula("value ~ year + country_code"), TFRpred, median) # using formula to trick the CRAN check
        TFRpred <- TFRpred[order(TFRpred$year),]
        country.table <- unique(merge(TFRpred[,'country_code', drop=FALSE], 
                           UNlocations[, c("country_code", "name")], by="country_code"))
        colnames(country.table)[1] <- "code"
    }
    # get pasfr patterns
    patterns <- .get.mig.mx.pasfr.patterns(inputs, wpp.year)
    PASFRpattern <- patterns$pasfr.pattern
    
    # prepare inputs
    inp <- new.env()
    for(par in c('PASFR', 'PASFRpattern', 'TFRpred', 'present.year', 'end.year', 'observed'))
        assign(par, get(par), envir=inp)
    
    # compute global norm
    inp$PASFRnorms <- compute.pasfr.global.norms(inp)
    
    # get last.observed
    lastobs <- NULL
    if(!is.null(inputs$last.observed))
        lastobs <- read.table(inputs$last.observed, sep = "\t", check.names = FALSE, header = TRUE)
    
    reslist <- NULL
    country.inputs <- new.env()
    # iterate over countries
    #for(icountry in which(country.table$name == "Fiji")) {
    for(icountry in 1:nrow(country.table)) {
        if(verbose && interactive()) cat('\r', "Computing PASFR: ", round(icountry/nrow(country.table)*100), '%')
        this.proj.years <- proj.years
        this.proj.periods <- proj.periods
        # extract country-specific inputs
        country <- country.table[icountry,]
        if (!is.data.frame(TFRpred)) {
            country.obj <- get.country.object(country$code, TFRpred$mcmc.set$meta)
            country.medians.TFRpred <- bayesTFR::get.median.from.prediction(TFRpred, country.obj$index, country.obj$code)[-1]
            country.obs.tfr <- bayesTFR:::get.tfr.reconstructed(TFRpred$tfr_matrix_reconstructed, TFRpred$mcmc.set$meta)
            country.obs.TFRpred <- country.obs.tfr[1:if(!is.null(TFRpred$present.year.index)) TFRpred$present.year.index else nrow(country.obs.tfr),country.obj$index]
            country.tfr <- c(country.obs.TFRpred, country.medians.TFRpred)
            match.years <- which(proj.years %in% names(country.tfr))
        } else {
            country.tfr.df <- TFRpred[TFRpred$country_code==country$code,]
            country.tfr <- country.tfr.df$value
            names(country.tfr) <- country.tfr.df$year
            match.years <- which(proj.years %in% country.tfr.df$year)
        }
        this.proj.years <- this.proj.years[match.years]
        this.proj.periods <- this.proj.periods[match.years]
        country.inputs$PASFRpattern <- .get.par.from.inputs('PASFRpattern', inp, country$code)
        country.inputs$observed <- list(PASFR=.get.par.from.inputs('PASFR', inp$observed, country$code))
        this.age <- inp$observed$PASFR[inp$observed$PASFR$country_code == country$code, "age"]
        if(is.null(country.inputs$observed$PASFR) || is.null(country.inputs$PASFRpattern)) next

        if(!is.null(lastobs)) { # add projection years if there is a gap between last.observed and first projection
            lastobs.cntry <- lastobs[lastobs$country_code == country$code, "last.observed"]
            if(lastobs.cntry + period.step < this.proj.years[1]) {
                add.years <- seq(lastobs.cntry + period.step, this.proj.years[1] - period.step, by = period.step)
                this.proj.years <- c(add.years, this.proj.years)
                if(annual) this.proj.periods <- this.proj.years
                else {
                    add.periods <- sort(seq(proj.years[1]-8, length = length(add.years), by = -period.step))
                    this.proj.periods <- c(paste(add.periods, add.periods+5, sep="-"), proj.periods)
                }
                country.inputs$observed$PASFR <- country.inputs$observed$PASFR[, !colnames(country.inputs$observed$PASFR) %in% this.proj.periods]
            }
        }
        # compute pasfr for the median TFR
        pasfr <- kantorova.pasfr(country.tfr, country.inputs, norms=inp$PASFRnorms, proj.years = this.proj.years, 
                                tfr.med=country.tfr[length(country.tfr)], annual = annual, nr.est.points = nr.est.points)
        # as percentage
        pasfr <- round(pasfr*100, digits)
        
        colnames(pasfr) <- this.proj.periods
        
        # collect results
        reslist[[icountry]] <- cbind(data.table(country_code=country$code, name=country$name, age=this.age),
                                     pasfr)
    }
    if(verbose && interactive()) cat('\n')
    
    # create one data frame
    resdf <- NULL
    for(i in 1:length(reslist)) 
        resdf <- rbind(resdf, reslist[[i]], fill = TRUE)
    non.int.cols <- c("country_code", "name", "age")
    resdf <- resdf[, c(non.int.cols, sort(setdiff(colnames(resdf), non.int.cols))), 
                   with = FALSE]
    # export
    if(!is.null(out.file.name)) {
      write.table(resdf, out.file.name, sep="\t", row.names=FALSE)
      cat("\nResults written into ", out.file.name, "\n")
    }
    invisible(data.frame(resdf, check.names = FALSE))
}


project.pasfr.traj <- function(inputs=NULL, countries=NULL, nr.traj=NULL, present.year=2020, end.year=2100, wpp.year=2019, 
                               digits=2, out.file.name="percentASFRtraj.txt") {
    # Generate trajectories of PASFR for given TFR using the kantorova.pasfr method.
    # ########
    observed <- list()
    ages <- seq(15, 45, by=5)
    age.labels <- paste(ages, ages+4, sep="-")
    
    # get SRB data
    srblist <- .get.srb.data.and.time.periods(inputs$srb, present.year, end.year, wpp.year)
    SRB <- srblist$srb
    observed$SRB <- srblist$obs.srb
    proj.periods <- srblist$proj.periods
    proj.years <- srblist$proj.years
    
    # get PASFR data
    pasfrlist <- .get.pasfr.data(inputs$pasfr, wpp.year, srblist$obs.periods, 
                                 srblist$proj.periods)
    PASFR <- pasfrlist$pasfr
    observed$PASFR <- pasfrlist$obs.pasfr
    
    # get TFR
    TFRpred <- .get.tfr.data(inputs, wpp.year)
    if(!is.data.frame(TFRpred)) { # object of class bayesTFR.prediction
        country.table <- get.countries.table(TFRpred)
    } else {
        obs.years <- as.integer(substr(srblist$obs.periods, 1, 4))+3
        TFRpred <- TFRpred[TFRpred$year %in% c(obs.years, proj.years),]
        #TFRpred <- aggregate(value ~ year + country_code, TFRpred, median)
        #TFRpred <- aggregate(formula("value ~ year + country_code"), TFRpred, median) # using formula to trick the CRAN check
        TFRpred <- TFRpred[order(TFRpred$year),]
        locs <- bayesTFR:::load.bdem.dataset('UNlocations', wpp.year)
        country.table <- unique(merge(TFRpred[,'country_code', drop=FALSE], 
                                      locs[, c("country_code", "name")], by="country_code"))
        colnames(country.table)[1] <- "code"
    }
    # get pasfr patterns
    patterns <- .get.mig.mx.pasfr.patterns(inputs, wpp.year)
    PASFRpattern <- patterns$pasfr.pattern
    
    inp <- new.env()
    for(par in c('SRB', 'PASFR', 'PASFRpattern', 'TFRpred', 'present.year', 'end.year', 'observed'))
        assign(par, get(par), envir=inp)
    inp$PASFRnorms <- compute.pasfr.global.norms(inp)
    
    resdf <- NULL
    country.inputs <- new.env()
    
    countries.idx <- if(is.null(countries)) 1:nrow(country.table) else which(country.table$code %in% countries)
    lcountries <- length(countries.idx)
    counter <- 0
    # iterate over countries
    for(icountry in countries.idx) {
        counter <- counter + 1
        # extract country-specific inputs
        country <- country.table[icountry,]
        if(interactive()) cat('\r', "Computing PASFR: ", round(counter/lcountries*100), '%')
        if (!is.data.frame(TFRpred)) {
            country.obj <- get.country.object(country$code, TFRpred$mcmc.set$meta)
            country.medians.TFRpred <- bayesTFR::get.median.from.prediction(TFRpred, country.obj$index, country.obj$code)
            last.median <- country.medians.TFRpred[length(country.medians.TFRpred)]
            country.obs.tfr <- bayesTFR:::get.tfr.reconstructed(TFRpred$tfr_matrix_reconstructed, TFRpred$mcmc.set$meta)[,country.obj$index]
            #country.obs.TFRpred <- country.obs.tfr[1:if(!is.null(TFRpred$present.year.index)) TFRpred$present.year.index else nrow(country.obs.tfr),country.obj$index]
            country.tfr.traj <- bayesTFR::get.tfr.trajectories(TFRpred, country=country.obj$code)
            #country.tfr <- c(country.obs.TFRpred, country.medians.TFRpred)
        } else {
            country.tfr.df <- TFRpred[TFRpred$country_code==country$code, c('year', 'trajectory', 'value')]
            country.tfr.traj <- country.tfr.df[country.tfr.df$year %in% c(max(obs.years), proj.years),]
            country.tfr.traj <- reshape(country.tfr.traj, direction="wide", idvar="year", timevar="trajectory")
            rownames(country.tfr.traj) <- country.tfr.traj$year
            country.tfr.traj <- country.tfr.traj[,-which(colnames(country.tfr.traj)=="year")]
            last.median <- apply(country.tfr.traj, 1, median)[nrow(country.tfr.traj)]
            country.obs <- country.tfr.df[country.tfr.df$year %in% obs.years,]
            # should be the same values over trajectories, so take mean to collapse to one record per year
            country.obs <- aggregate(value ~ year, country.obs, mean) 
            country.obs.tfr <- country.obs$value
            names(country.obs.tfr) <- country.obs$year
        }
        country.inputs$PASFRpattern <- .get.par.from.inputs('PASFRpattern', inp, country$code)
        country.inputs$observed <- list(PASFR=.get.par.from.inputs('PASFR', inp$observed, country$code))
        if(is.null(country.inputs$observed$PASFR) || is.null(country.inputs$PASFRpattern)) next
        # compute pasfr for each TFR trajectory
        country.pasfr <- NULL
        traj.idx <- if(is.null(nr.traj)) 1:ncol(country.tfr.traj) else unique(as.integer(seq(1, ncol(country.tfr.traj), length=nr.traj)))
        for (itraj in traj.idx) {
            proj.tfr.traj <- country.tfr.traj[-1,itraj]
            names(proj.tfr.traj) <- rownames(country.tfr.traj)[-1]
            pasfr <- kantorova.pasfr(c(country.obs.tfr, proj.tfr.traj), country.inputs, norms=inp$PASFRnorms, proj.years=proj.years, 
                                 tfr.med=last.median)
            # as percentage
            pasfr <- round(pasfr*100, digits)
            colnames(pasfr) <- proj.years
            rownames(pasfr) <- ages
            pasfr.long <- reshape2::melt(pasfr, value.name=as.character(itraj))
            colnames(pasfr.long)[1:2] <- c("age", "year")
            country.pasfr <- if(is.null(country.pasfr)) pasfr.long else merge(country.pasfr, 
                                                pasfr.long, by=c('age', 'year'), sort=FALSE)
        }
        country.pasfr <- cbind(data.frame(country_code=country$code, name=country$name), country.pasfr)
        # add to other results
        resdf <- rbind(resdf, country.pasfr)
    }
    if(interactive()) cat('\n')
    # export
    if(!is.null(out.file.name)) {
        write.table(resdf, out.file.name, sep="\t", row.names=FALSE)
        cat("\nResults written into ", out.file.name, "\n")
    }
    invisible(resdf)
}
