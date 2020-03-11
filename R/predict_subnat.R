if(getRversion() >= "2.15.1") utils::globalVariables(c("UNlocations"))

pop.predict.subnat <- function(end.year = 2060, start.year = 1950, present.year = 2020, wpp.year = 2019,
                                output.dir = file.path(getwd(), "bayesPop.output"),
                               locations = NULL, default.country = NULL,
                        inputs = list(
                            popM = NULL,
                            popF = NULL,
                            mxM = NULL,
                            mxF = NULL,
                            srb = NULL,
                            pasfr = NULL,
                            patterns = NULL,
                            migM = NULL,
                            migF = NULL,	
                            e0F.file = NULL, e0M.file = NULL, 
                            tfr.file = NULL, 
                            e0F.sim.dir = NULL, e0M.sim.dir = NULL, 
                            tfr.sim.dir = NULL,
                            migMtraj = NULL, migFtraj = NULL	
                        ), nr.traj = 1000, keep.vital.events = FALSE,
                        fixed.mx = FALSE, fixed.pasfr = FALSE, 
                        replace.output = FALSE, verbose = TRUE) {
    #prediction.exist <- FALSE
    ages <- seq(0, by = 5, length = 27)

    if(is.null(locations)) 
        stop("Argument locations must be given.")
        
    #UNlocations <- NULL # needed for R check not to complain
    locs <- process.reg.locations(locations, verbose = verbose)
    if (is.null(default.country) && any(locs$location_type == 0))
        default.country <- locs[locs$location_type == 0, "country_code"][1]
    env <- globalenv()
    assign("UNlocations", locs, envir=env)
    
    outdir <- file.path(output.dir, 'predictions')
    if(has.pop.prediction(sim.dir = output.dir)) {
        if(!replace.output)
            stop('Prediction in ', output.dir,
                ' already exists.\nSet replace.output = TRUE if you want to overwrite existing projections.')
        unlink(outdir, recursive=TRUE)
        .remove.cache.file(outdir)
    }
    inp <- load.subnat.inputs(inputs, start.year, present.year, end.year, wpp.year, 
                              default.country = default.country,
                              fixed.mx = fixed.mx, fixed.pasfr = fixed.pasfr, verbose = verbose)

    reg.codes <- intersect(unique(inp$POPm0[,'country_code']), UNcountries())
    
    do.pop.predict(reg.codes, inp, outdir, nr.traj, ages, 
                   keep.vital.events = keep.vital.events, fixed.mx = inp$fixed.mx, fixed.pasfr = fixed.pasfr, 
                   function.inputs = inputs, verbose = verbose)
    invisible(get.pop.prediction(output.dir))
}

process.reg.locations <- function(locations, verbose = FALSE) {
    locs <- locations
    if(is.character(locations)) {
        locs <- read.delim(file = locations, comment.char = '#', check.names = FALSE)
        if(verbose) cat('Loading ', locations, '.\n')
    }
    # swap reg_code for country_code
    locs <- swap.reg.code(locs, table.name = "locations")

    # add location_type column if not present
    if(! "location_type" %in% colnames(locs))
        locs$location_type <- 4
    locs
}

swap.reg.code <- function(ds, table.name = "input") {
    # swap reg_code for country_code
      if(!"reg_code" %in% colnames(ds))
          stop(table.name, " table must have a column 'reg_code'.")
    if("country_code" %in% colnames(ds)) ds$country_code <- NULL # delete if exists
    colnames(ds)[colnames(ds) == "reg_code"] <- "country_code" # swap
    ds
}

repeat.dataset.for.regions <- function(ds, region.codes) {
    for(col in c("country_code", "country", "name")) # delete country columns
        ds[[col]] <- NULL
    res <- NULL
    for(code in region.codes) {
        res <- rbind(res, cbind(country_code = code, ds))
    }
    res
}

create.dataset.from.wpp <- function(country_code, region.codes, ...) {
    ds <- load.wpp.dataset.for.country(country_code, ...)
    repeat.dataset.for.regions(ds, region.codes)
}

load.wpp.dataset.for.country <- function(cntry_id, ...) {
    ds <- load.wpp.dataset(...)
    return(ds[ds$country_code == cntry_id, ])
}

create.trajs.from.wpp <- function(country_code, region.codes, ...) {
    ds <- load.wpp.traj.for.country(country_code, ...)
    repeat.dataset.for.regions(ds, region.codes)
}

load.wpp.traj.for.country <- function(cntry_id, ...) {
    ds <- .load.wpp.traj(...)
    return(ds[ds$country_code == cntry_id, ])
}

.get.tfr.data.subnat <- function(inputs,  wpp.year, default.country, region.codes, 
                                 verbose = FALSE) {
    if(!is.null(inputs$tfr.file)) {
        if(inputs$tfr.file == 'median_')
            TFRpred <- create.trajs.from.wpp(default.country, 
                                             region.codes, 'tfr', wpp.year, 
                                             median.only = TRUE)
        else {
            file.name <- inputs$tfr.file
            if(!file.exists(file.name))
                stop('File ', file.name, 
                     ' does not exist.\nSet tfr.sim.dir, tfr.file or change WPP year.')
            if(verbose) cat('\nLoading ', file.name, '\n')
            TFRpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
            TFRpred <- TFRpred[,c('LocID', 'Year', 'Trajectory', 'TF')]
            colnames(TFRpred) <- c('country_code', 'year', 'trajectory', 'value')
        } 
    } else {
        if(!is.null(inputs$tfr.sim.dir)) 
            TFRpred <- get.tfr.prediction(inputs$tfr.sim.dir, mcmc.dir=NA)
        else TFRpred <- create.trajs.from.wpp(default.country, 
                                              region.codes, 'tfr', wpp.year)
    }
    return(TFRpred)
}

.get.mig.shares <- function(ds, region.codes) {
    if(length(grep("migration[[:alpha:]]*_share", colnames(ds))) == 0) return(NULL)
    # put ds into the same order as region.codes
    ds <- merge(data.frame(country_code = region.codes), ds, by = "country_code", sort = FALSE)
    MIGshare <- list(inmig = list(M = NULL, F = NULL))
    MIGshare$outmig  <- MIGshare$inmig
    colname.def1 <- paste0("migration_share")
    default <- NULL
    if(colname.def1 %in% colnames(ds))
        default <- ds[, colname.def1]
    for(sex in c("M", "F")) {
        this.def <- default
        colname.def2 <- paste0("migration", sex, "_share")
        if(colname.def2 %in% colnames(ds))
            this.def <- ds[, colname]
        for(what in c("in", "out")) {
            whatname <- paste0(what, "mig")
            assigned <- FALSE
            for(colname in c(paste0(what, "migration", sex, "_share"), paste0(what, "migration_share"))) {
                if(colname %in% colnames(ds)) {
                    MIGshare[[whatname]][[sex]] <- ds[, colname]
                    assigned <- TRUE
                    break
                }
            }
            if(!assigned)
                MIGshare[[whatname]][[sex]] <- this.def
        }
    }
    # scale
    for(what in c("inmig", "outmig")) {
        for(sex in c("M", "F")) {
            if(is.null(MIGshare[[what]][[sex]]))
                stop(what, "ration shares for sex ", sex, " missing while other shares available!")
        }
        s <- sum(MIGshare[[what]]$M + MIGshare[[what]]$F)
        MIGshare[[what]]$M <- MIGshare[[what]]$M/s
        MIGshare[[what]]$F <- MIGshare[[what]]$F/s
    }
    return(MIGshare)
}

load.subnat.inputs <- function(inputs, start.year, present.year, end.year, wpp.year, 
                               default.country = NULL, fixed.mx = FALSE, 
                        fixed.pasfr = FALSE, existing.mig = NULL, verbose = FALSE) {
    observed <- list()
    # check inputs
    for(item in c('popM', 'popF')) {
        if(is.null(inputs[[item]]))
            stop("Item ", item, " is missing in the 'inputs' argument.")
    }
    pop.ini.matrix <- pop.ini <- list(M = NULL, F = NULL)
    # Get initial population counts
    for(sex in c('M', 'F')) {
        dataset.name <- paste0('pop', sex)
        POP0 <- read.pop.file(inputs[[dataset.name]])
        num.columns <- grep('^[0-9]{4}$', colnames(POP0), value = TRUE) # values of year-columns
        if(!is.element(as.character(present.year), colnames(POP0))) {
            stop('Wrong present.year. ', present.year, ' not available in the ', dataset.name, ' file.\nAvailable years: ',
                 paste(num.columns, collapse = ', '))
        }
        POP0[is.na(POP0)] <- 0
        num.columns <- num.columns[which(as.integer(num.columns) <= present.year)]
        pop.ini.matrix[[sex]] <- POP0[, num.columns, drop = FALSE]
        dimnames(pop.ini.matrix[[sex]]) <- list(paste(POP0[,'reg_code'], POP0[,'age'], sep = '_'), 
                                                as.character(as.integer(num.columns)))
        pop.ini[[sex]] <- POP0[ ,c('reg_code', 'age', present.year)]
        colnames(pop.ini[[sex]])[1] <- 'country_code'
    }
    POPm0 <- pop.ini[['M']]
    POPf0 <- pop.ini[['F']]
    region.codes <- unique(POPm0$country_code)
    # Get death rates
    MXm.pred <- MXf.pred <- NULL
    if(is.null(inputs$mxM)) {
        if(!is.null(default.country))
            MXm <- create.dataset.from.wpp(default.country, region.codes, 'mxM', wpp.year)
        else stop("mxM must be given if there is no default.country.")
    } else MXm <- swap.reg.code(read.pop.file(inputs$mxM), table.name = "mxM")
    names.MXm.data <- names(MXm)
    num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.MXm.data) # index of year-columns
    cols.starty <- as.integer(substr(names.MXm.data[num.columns], 1,4))
    cols.endy <- as.integer(substr(names.MXm.data[num.columns], 6,9))
    start.index <- which((cols.starty <= start.year) & (cols.endy > start.year))
    present.index <- which((cols.endy >= present.year) & (cols.starty < present.year))
    #estim.periods <- names.MXm.data[num.columns[start.index:present.index]]
    estim.periods <- names.MXm.data[num.columns[1:present.index]] # ?why
    start.year <- cols.starty[start.index]
    
    if(fixed.mx) {
        end.index <- which((cols.endy >= end.year) & (cols.starty < end.year))
        proj.periods <- names.MXm.data[num.columns[(present.index+1):end.index]]
        MXm.pred <- MXm[,c('country_code', 'age', proj.periods)]
    }
    MXm <- MXm[,c('country_code', 'age', estim.periods)]
    if(is.null(inputs$mxF)) {
      if(!is.null(default.country))
          MXf <- create.dataset.from.wpp(default.country, region.codes, 'mxF', wpp.year)
      else stop("mxF must be given if there is no default.country.")
    } else MXf <- swap.reg.code(read.pop.file(inputs$mxF), table.name = "mxF")
    if(fixed.mx) MXf.pred <- MXf[,c('country_code', 'age', proj.periods)]
    MXf <- MXf[,c('country_code', 'age', estim.periods)]
    
    estim.years <- cols.starty[start.index:present.index] + 3
    
    # Get sex ratio at birth
    srb.data <- inputs$srb
    if(is.null(srb.data)) 
        srb.data <- create.dataset.from.wpp(default.country, region.codes, 'sexRatio', wpp.year)
    else srb.data <- swap.reg.code(read.pop.file(srb.data), table.name = "srb")
    srblist <- .get.srb.data.and.time.periods(srb.data, present.year, end.year, wpp.year)
    SRB <- srblist$srb
    observed$SRB <- srblist$obs.srb
    proj.periods <- srblist$proj.periods
    obs.periods <- srblist$obs.periods
    proj.years <- srblist$proj.years
    
    # Get percentage age-specific fertility rate
    pasfr.data <- inputs$pasfr
    if(is.null(pasfr.data)) 
        pasfr.data <- create.dataset.from.wpp(default.country, region.codes, 'percentASFR', wpp.year)
    else pasfr.data <- swap.reg.code(read.pop.file(inputs$pasfr), table.name = "pasfr")
    pasfrlist <- .get.pasfr.data(pasfr.data, wpp.year, obs.periods, proj.periods, 
                                 include.projection = fixed.pasfr)
    PASFR <- pasfrlist$pasfr
    observed$PASFR <- pasfrlist$obs.pasfr
    
    # Get migration type, migration base year, mx & pasfr patterns
    pattern.data.def <- get(paste0('vwBaseYear', wpp.year))
    pattern.data.def <- repeat.dataset.for.regions(
                            pattern.data.def[pattern.data.def$country_code == default.country,], 
                            region.codes)
    if(!is.null(inputs$patterns)) {
        pattern.data <- swap.reg.code(read.pop.file(inputs$patterns), table.name = "patterns")
        for(col in colnames(pattern.data.def)) {
            # fill-in missing columns
            if(! col %in% colnames(pattern.data))
                pattern.data <- merge(pattern.data, pattern.data.def[, c("country_code", col)], by = "country_code")
        }
    } else pattern.data <- pattern.data.def
    patterns <- .get.mig.mx.pasfr.patterns(inputs, wpp.year, pattern.data = pattern.data)
    MIGtype <- patterns$mig.type
    MXpattern <- patterns$mx.pattern
    PASFRpattern <- patterns$pasfr.pattern
    MIGshare <- .get.mig.shares(pattern.data, region.codes)
    
    # Get age-specific migration
    wppds <- data(package = paste0('wpp', wpp.year))
    recon.mig <- NULL
    if(is.null(inputs[['migM']]) || is.null(inputs[['migF']])) {
        #MIGm <- MIGf <- create.dataset.from.wpp(default.country, region.codes, 'migrationM', 2012) # to have the right dataset format
        #numcols <- colnames(MIGm)[-which(colnames(MIGm) %in% c("country", "name", "country_code", "age"))]
        mignat <- load.wpp.dataset.for.country(default.country, 'migration', wpp.year)
        mignat <- mignat[,-which(colnames(mignat) %in% c("country", "name", "country_code"))]
        rogcastro <- abs(load.wpp.dataset.for.country(156, 'migrationM', 2012)[, "2015-2020"])
        rc <- rogcastro/sum(rogcastro)
        migdistr <- mignat[rep(1, 21),] * rc
        nreg <- length(region.codes)
        nrc <- length(rc)
        nT <- length(mignat)
        popprop <- matrix(0, nrow = 2*nreg, ncol = nrc)
        migmtxM <- migmtxF <- NULL
        migages <- as.character(POPm0$age[1:21])
        which.mig.share <- rep("inmig", length(mignat))
        which.mig.share[mignat < 0] <- "outmig"
        if(!is.null(MIGshare)) { # migration shares given
            scaleM <- scaleF <- matrix(0, nrow = nreg, ncol = nT)
            for(iy in 1:nT) {
                scaleM[, iy] <- MIGshare[[which.mig.share[iy]]]$M
                scaleF[, iy] <- MIGshare[[which.mig.share[iy]]]$F
            }
        }
        #stop('')
        for(iage in 1:21) {
            popidx <- which(POPm0$age == migages[iage])
            if(is.null(MIGshare)) { # migration shares not given
                popprop[1:nreg, iage] <- POPm0[popidx, 3]
                popprop[(nreg+1):(2*nreg), iage] <- POPf0[POPf0$age == migages[iage], 3]
                popprop[, iage] <- popprop[, iage]/sum(popprop[, iage])
                scaleM <- popprop[1:nreg, rep(iage, nT)]
                scaleF <- popprop[(nreg+1):(2*nreg), rep(iage, nT)]
            }
            thismigM <- migdistr[rep(iage, nreg), ] * scaleM
            thismigF <- migdistr[rep(iage, nreg), ] * scaleF
                
            #thismigM <- thismigF <- NULL
            #for(y in colnames(mignat)) {
            #    thismigM <- cbind(thismigM, migdistr[iage,y] * popprop[1:nreg, iage])
            #    thismigF <- cbind(thismigF, migdistr[iage,y] * popprop[(nreg+1):(2*nreg), iage])
            #}
            #colnames(thismigM) <- colnames(thismigF) <- colnames(mignat)
            migmtxM <- rbind(migmtxM, data.frame(country_code = POPm0[popidx, "country_code"],
                                                 age = migages[iage], thismigM))
            migmtxF <- rbind(migmtxF, data.frame(country_code = POPm0[popidx, "country_code"],
                                                 age = migages[iage], thismigF))
        }
        MIGm <- migmtxM[order(migmtxM$country_code),]
        MIGf <- migmtxF[order(migmtxF$country_code),]
        colnames(MIGm)[3:ncol(MIGm)] <- colnames(MIGf)[3:ncol(MIGf)] <- colnames(mignat)
        #MIGm[,3:ncol(MIGm)] <- 0
        #MIGf[,3:ncol(MIGf)] <- 0
        #stop('')
    }
    if(!is.null(inputs[['migM']])) # cannot be inputs$migM because it would take the value of inputs$migMpred
        MIGm <- swap.reg.code(read.pop.file(inputs[['migM']]), table.name = "migM")
    if(!is.null(inputs[['migF']])) 
        MIGf <- swap.reg.code(read.pop.file(inputs[['migF']]), table.name = "migF")
    if(!is.null(obs.periods)) {
        avail.obs.periods <- is.element(obs.periods, colnames(MIGm))
        observed$MIGm <- MIGm[,c('country_code', 'age', obs.periods[avail.obs.periods])]
        observed$MIGf <- MIGf[,c('country_code', 'age', obs.periods[avail.obs.periods])]
    }
    MIGm <- MIGm[,c('country_code', 'age', proj.periods)]
    MIGf <- MIGf[,c('country_code', 'age', proj.periods)]
    # Get migration trajectories if available
    migMpred <- migFpred <- NULL
    for(sex in c('M', 'F')) {
        if(is.null(inputs[[paste0('mig',sex,'traj')]])) next
        file.name <- inputs[[paste0('mig',sex,'traj')]]
        if(!file.exists(file.name))
            stop('File ', file.name, ' does not exist.')
        # comma separated trajectories file
        var.name <- paste0('mig',sex, 'pred')
        if(verbose) cat('\nLoading ', file.name)
        migpred.raw <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
        migpred <- migpred.raw[,c('LocID', 'Year', 'Trajectory', 'Age', 'Migration')]
        colnames(migpred) <- c('country_code', 'year', 'trajectory', 'age', 'value')
        assign(var.name, migpred)
    }
    
    # Get life expectancy
    e0F.wpp.median.loaded <- FALSE
    e0Fpred <- e0Mpred <- NULL
    if(!fixed.mx){
        if(!is.null(inputs$e0F.file)) { # female
            if(inputs$e0F.file == 'median_') {
                e0Fpred <- create.trajs.from.wpp(default.country, 
                                                 region.codes, 'e0F', wpp.year, 
                                                 median.only = TRUE)
                e0F.wpp.median.loaded <- TRUE
            } else {
                file.name <-  inputs$e0F.file
                if(!file.exists(file.name))
                    stop('File ', file.name, ' does not exist.')
                # comma separated trajectories file
                if(verbose) cat('\nLoading ', file.name)
                e0Fpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
                e0Fpred <- e0Fpred[,c('LocID', 'Year', 'Trajectory', 'e0')]
                colnames(e0Fpred) <- c('country_code', 'year', 'trajectory', 'value')
            }
        } else {
            if(!is.null(inputs$e0F.sim.dir)) { 
                if(inputs$e0F.sim.dir == 'median_') {
                    e0Fpred <- create.trajs.from.wpp(default.country, 
                                                     region.codes, 'e0F', wpp.year, 
                                                     median.only = TRUE)
                    e0F.wpp.median.loaded <- TRUE
                } else 
                    e0Fpred <- get.e0.prediction(inputs$e0F.sim.dir, mcmc.dir=NA)
            } else e0Fpred <- create.trajs.from.wpp(default.country, 
                                                    region.codes, 'e0F', wpp.year)			
        }
        
        if(!is.null(inputs$e0M.file)) { # male
            if(inputs$e0M.file == 'median_')
                e0Mpred <- create.trajs.from.wpp(default.country, 
                                                 region.codes, 'e0M', wpp.year, 
                                                 median.only = TRUE)
            else {
                file.name <-  inputs$e0M.file
                if(!file.exists(file.name)) 
                    stop('File ', file.name, ' does not exist.')
                if(verbose) cat('\nLoading ', file.name)
                e0Mpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
                e0Mpred <- e0Mpred[,c('LocID', 'Year', 'Trajectory', 'e0')]
                colnames(e0Mpred) <- c('country_code', 'year', 'trajectory', 'value')
            }
        } else {
            if(!is.null(inputs$e0M.sim.dir)) { 
                if(inputs$e0M.sim.dir == 'joint_') {
                    if(e0F.wpp.median.loaded) 
                        e0Mpred <- create.trajs.from.wpp(default.country, region.codes, 
                                                         'e0M', wpp.year)
                    else {
                        if(!has.e0.jmale.prediction(e0Fpred))
                            stop('No joint prediction for female and male available. Correct the e0M.sim.dir argument.' )
                        e0Mpred <- get.e0.jmale.prediction(e0Fpred)
                    }
                } else e0Mpred <- get.e0.prediction(inputs$e0M.sim.dir, mcmc.dir=NA) # independent from female
            } else
                e0Mpred <- create.trajs.from.wpp(default.country, 
                                                 region.codes, 'e0M', wpp.year)
        }
    } # end if(!fixed.mx)
    
    # Get TFR
    TFRpred <- .get.tfr.data.subnat(inputs, wpp.year, default.country, region.codes, verbose=verbose)
    inp <- new.env()
    for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MXm.pred', 'MXf.pred', 'MXpattern', 'SRB',
                 'PASFR', 'PASFRpattern', 'MIGtype', 'MIGm', 'MIGf',
                 'e0Mpred', 'e0Fpred', 'TFRpred', 'migMpred', 'migFpred', 'estim.years', 'proj.years', 'wpp.year', 
                 'start.year', 'present.year', 'end.year', 'fixed.mx', 'fixed.pasfr', 'observed'))
        assign(par, get(par), envir=inp)
    inp$pop.matrix <- list(male=pop.ini.matrix[['M']], female=pop.ini.matrix[['F']])
    #stop('')
    env <- new.env()
    do.call("data", list("pasfr_global_norms", envir = env))
    inp$PASFRnorms <- env$pasfr.glob.norms
    inp$lc.for.hiv <- TRUE
    inp$lc.for.all <- TRUE
    return(inp)
}

pop.aggregate.subnat <- function(pop.pred, regions, locations, ..., verbose = FALSE) {
    locs <- process.reg.locations(locations, verbose = verbose)
    pop.aggregate(pop.pred, regions, my.location.file = locs, ...)
}