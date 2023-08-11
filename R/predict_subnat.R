if(getRversion() >= "2.15.1") utils::globalVariables(c("UNlocations"))

pop.predict.subnat <- function(end.year = 2060, start.year = 1950, present.year = 2020, wpp.year = 2019,
                                output.dir = file.path(getwd(), "bayesPop.output"),
                               locations = NULL, default.country = NULL, annual = FALSE,
                        inputs = list(
                            popM = NULL, popF = NULL,
                            mxM = NULL, mxF = NULL, srb = NULL,
                            pasfr = NULL, patterns = NULL,
                            migM = NULL, migF = NULL,	
                            migMt = NULL, migFt = NULL, mig = NULL,
                            e0F.file = NULL, e0M.file = NULL, tfr.file = NULL, 
                            e0F.sim.dir = NULL, e0M.sim.dir = NULL, 
                            tfr.sim.dir = NULL,
                            migMtraj = NULL, migFtraj = NULL, migtraj = NULL,
                            GQpopM = NULL, GQpopF = NULL, average.annual = NULL
                        ), nr.traj = 1000, keep.vital.events = FALSE,
                        fixed.mx = FALSE, fixed.pasfr = FALSE, lc.for.all = TRUE, mig.is.rate = FALSE,
                        replace.output = FALSE, verbose = TRUE) {
    #prediction.exist <- FALSE
    ages <- ages.all(annual, observed = FALSE)

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
                              fixed.mx = fixed.mx, fixed.pasfr = fixed.pasfr, 
                              lc.for.all = lc.for.all, mig.is.rate = mig.is.rate,
                              annual = annual, verbose = verbose)

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
                        fixed.pasfr = FALSE, lc.for.all = TRUE, mig.is.rate = FALSE,
                        annual = FALSE, verbose = FALSE) {
    observed <- list()
    # check inputs
    for(item in c('popM', 'popF')) {
        if(is.null(inputs[[item]]))
            stop("Item ", item, " is missing in the 'inputs' argument.")
    }
    pop.ini.matrix <- pop.ini <- GQ <- list(M = NULL, F = NULL)
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
        # Group quarters
        dataset.name <- paste0('GQpop', sex)
        if(!is.null(inputs[[dataset.name]])) {
            GQ[[sex]] <- read.pop.file(inputs[[dataset.name]])
            colnames(GQ[[sex]]) <- tolower(colnames(GQ[[sex]]))
            if(! 'age' %in% colnames(GQ[[sex]]) || ! 'reg_code' %in% colnames(GQ[[sex]]) || ! 'gq' %in% colnames(GQ[[sex]]))
                stop('Columns "age", "reg_code" and "gq" must be present in the GQpop datasets.')
            GQ[[sex]] <- GQ[[sex]][, c("reg_code", "age", "gq")]
            colnames(GQ[[sex]])[1] <- 'country_code'
        }
    }
    POPm0 <- pop.ini[['M']]
    POPf0 <- pop.ini[['F']]
    GQm <- GQ[['M']]
    GQf <- GQ[['F']]
    
    region.codes <- unique(POPm0$country_code)
    # Get death rates
    MXm.pred <- MXf.pred <- NULL
    if(is.null(inputs$mxM)) {
        if(!is.null(default.country))
            MXm <- create.dataset.from.wpp(default.country, region.codes, 'mxM', wpp.year)
        else stop("mxM must be given if there is no default.country.")
    } else MXm <- swap.reg.code(read.pop.file(inputs$mxM), table.name = "mxM")
    names.MXm.data <- names(MXm)
    numcol.expr <- if(annual) '^[0-9]{4}$' else '^[0-9]{4}.[0-9]{4}$'
    num.columns <- grep(numcol.expr, names.MXm.data) # index of year-columns
    if(length(num.columns) == 0) stop("Column names of numeric columns of mx are not in the right format. Use ",
                                      if(annual) "XXXX, e.g. 1963" else "XXXX-XXXX, e.g. 1960-1965" )
    cols.starty <- as.integer(substr(names.MXm.data[num.columns], 1,4))
    if(!annual) {
        cols.endy <- as.integer(substr(names.MXm.data[num.columns], 6,9))
        start.index <- which((cols.starty <= start.year) & (cols.endy > start.year))
        if(length(start.index) == 0) {
          if(start.year <= cols.starty[1]) start.index <- 1
          else stop("start.year not available in the mx dataset")
        }
        present.index <- which((cols.endy >= present.year) & (cols.starty < present.year))
    } else {
        cols.endy <- cols.starty
        start.index <- which(cols.starty >= start.year)[1]
        present.index <- which(cols.starty == present.year)
    }
    if(length(present.index) == 0) stop("present.year not available in the mx dataset")
    #estim.periods <- names.MXm.data[num.columns[start.index:present.index]]
    estim.periods <- names.MXm.data[num.columns[1:present.index]] # ?why
    start.year <- cols.starty[start.index]
    
    if(fixed.mx) {
        end.index <- if(annual) which(cols.endy == end.year) else which((cols.endy >= end.year) & (cols.starty < end.year))
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
    
    estim.years <- cols.starty[start.index:present.index]
    if(!annual) estim.years <- estim.years + 3
    
    # Get sex ratio at birth
    srb.data <- inputs$srb
    if(is.null(srb.data)) 
        srb.data <- create.dataset.from.wpp(default.country, region.codes, 'sexRatio', wpp.year)
    else srb.data <- swap.reg.code(read.pop.file(srb.data), table.name = "srb")
    srblist <- .get.srb.data.and.time.periods(srb.data, present.year, end.year, wpp.year, annual = annual)
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
    
    miginp <- .get.mig.data.subnat(inputs, wpp.year, annual, periods = c(estim.periods, proj.periods), 
                            default.country = default.country, region.codes = region.codes,
                            pop0 = list(M = POPm0, F = POPf0), MIGshare = MIGshare,
                            verbose = verbose)
    
    MIGm <- miginp[["migM"]]
    MIGf <- miginp[["migF"]]
    
    if(!is.null(obs.periods)) {
        avail.obs.periods <- is.element(obs.periods, colnames(MIGm))
        observed$MIGm <- MIGm[,c('country_code', 'age', obs.periods[avail.obs.periods])]
        observed$MIGf <- MIGf[,c('country_code', 'age', obs.periods[avail.obs.periods])]
    }
    MIGm <- MIGm[,c('country_code', 'age', proj.periods)]
    MIGf <- MIGf[,c('country_code', 'age', proj.periods)]
    # Get migration trajectories if available
    migpr <- .load.mig.traj(inputs, verbose = verbose)
    migMpred <- migpr$M
    migFpred <- migpr$F
    
    if(length(mig.is.rate) < 2) mig.is.rate <- rep(mig.is.rate, 2) # one for observed, one for projection
    mig.rate.code <- c(miginp[["migcode"]]*mig.is.rate[1], migpr[["migcode"]]*mig.is.rate[2])
    
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
                 'PASFR', 'PASFRpattern', 'MIGtype', 'MIGm', 'MIGf', 'GQm', 'GQf',
                 'e0Mpred', 'e0Fpred', 'TFRpred', 'migMpred', 'migFpred', 'estim.years', 'proj.years', 'wpp.year', 
                 'start.year', 'present.year', 'end.year', 'annual', 'fixed.mx', 'fixed.pasfr', 
                 'lc.for.all', 'mig.rate.code', 'observed'))
        assign(par, get(par), envir=inp)
    inp$pop.matrix <- list(male=pop.ini.matrix[['M']], female=pop.ini.matrix[['F']])
    env <- new.env()
    do.call("data", list("pasfr_global_norms", envir = env))
    inp$PASFRnorms <- env$pasfr.glob.norms
    inp$lc.for.hiv <- TRUE
    inp$average.annual <- inputs$average.annual
    return(inp)
}

.get.mig.data.subnat <- function(inputs, wpp.year, annual, periods, default.country, region.codes, 
                                 pop0 = NULL, MIGshare = NULL, verbose = FALSE) {
  # Get age-specific migration
  wppds <- data(package=paste0('wpp', wpp.year))
  recon.mig <- NULL
  miginp <- list()
  migtempl <- NULL
  migcode <- 0
  for(sex in c("M", "F")) {
    inpname <- paste0('mig', sex) 
    if(!is.null(inputs[[inpname]])) { # migration given by sex and age
      miginp[[inpname]] <- swap.reg.code(read.pop.file(inputs[[inpname]]), table.name = inpname)
      next
    }
    # create a template to be filled with derived migration
    if(is.null(migtempl)) {
        # Here create only a dataframe filled with NAs 
        migtempl <- .get.mig.template(unique(pop0[[sex]]$country_code), 
                                        ages = pop0[[sex]]$age[age.index.all(annual, observed = TRUE)],
                                        time.periods = periods)
        migtempl[,which(!colnames(migtempl) %in% c("country", "country_code", "age"))] <- NA
        migtempl <- data.table(migtempl)
    }
    fname <- paste0(inpname, 't')
    fnametot <- paste0("mig")
    if(!is.null(inputs[[fname]]) || !is.null(inputs[[fnametot]])) { # migration given as totals or totals by sex
      if(!is.null(inputs[[fname]])) {
        totmig <- data.table(swap.reg.code(read.pop.file(inputs[[fname]])), table.name = fname)
        migcode <- 2
      } else {
        totmig <- data.table(swap.reg.code(read.pop.file(inputs[[fnametot]])), table.name = fnametot)
        migcode <- 3
      }
      migcols <- intersect(colnames(totmig), periods)
      miginp[[inpname]] <- data.frame(migration.totals2age(totmig, ages = migtempl$age[age.index.all(annual, observed = TRUE)],
                                                           annual = annual, time.periods = migcols, 
                                                           scale = if(is.null(inputs[[fname]])) 0.5 else 1, # since the totals are sums over sexes
                                                           template = migtempl), check.names = FALSE)
      next
    }
    # If we get here, migration is not given. Thus, get it from the national values and given shares
    if(annual) stop("Migration must be given.")
    mignat <- load.wpp.dataset.for.country(default.country, 'migration', wpp.year) # TODO: Allow input of annual national values 
    migdistr <- migration.totals2age(mignat, annual = annual)
    migages <- migdistr$age
    mignat <- mignat[,-which(colnames(mignat) %in% c("country", "name", "country_code"))]
    migdistr <- as.matrix(migdistr[,! colnames(migdistr) %in% c("country", "name", "country_code", "age"), with = FALSE])
    nreg <- length(region.codes)
    nrc <- nrow(migdistr)
    nT <- ncol(migdistr)
    popprop <- matrix(0, nrow = 2*nreg, ncol = nrc)
    migmtxM <- migmtxF <- NULL
    which.mig.share <- rep("inmig", nT)
    which.mig.share[mignat < 0] <- "outmig"
    if(!is.null(MIGshare)) { # migration shares given
      scaleM <- scaleF <- matrix(0, nrow = nreg, ncol = nT)
      for(iy in 1:nT) {
        scaleM[, iy] <- MIGshare[[which.mig.share[iy]]]$M
        scaleF[, iy] <- MIGshare[[which.mig.share[iy]]]$F
      }
    }
    for(iage in seq_along(migages)) {
      popidx <- which(pop0[["M"]]$age == migages[iage])
      if(is.null(MIGshare)) { # migration shares not given, take population
        popprop[1:nreg, iage] <- pop0[["M"]][popidx, 3]
        popprop[(nreg+1):(2*nreg), iage] <- pop0[["F"]][pop0[["F"]]$age == migages[iage], 3]
        popprop[, iage] <- popprop[, iage]/sum(popprop[, iage])
        scaleM <- popprop[1:nreg, rep(iage, nT)]
        scaleF <- popprop[(nreg+1):(2*nreg), rep(iage, nT)]
      }
      thismigM <- migdistr[rep(iage, nreg), ] * scaleM
      thismigF <- migdistr[rep(iage, nreg), ] * scaleF
      
      migmtxM <- rbind(migmtxM, data.frame(country_code = pop0[["M"]][popidx, "country_code"],
                                           age = migages[iage], thismigM))
      migmtxF <- rbind(migmtxF, data.frame(country_code = pop0[["M"]][popidx, "country_code"],
                                           age = migages[iage], thismigF))
    }
    MIGm <- migmtxM[order(migmtxM$country_code),]
    MIGf <- migmtxF[order(migmtxF$country_code),]
    colnames(MIGm)[3:ncol(MIGm)] <- colnames(MIGf)[3:ncol(MIGf)] <- colnames(mignat)
    miginp[["migM"]] <- MIGm
    miginp[["migF"]] <- MIGf
    break # no need to go to the next iteration as we have everything we need
  }
  miginp[["migcode"]] <- migcode
  return(miginp)
}


pop.aggregate.subnat <- function(pop.pred, regions, locations, ..., verbose = FALSE) {
    locs <- process.reg.locations(locations, verbose = verbose)
    pop.aggregate(pop.pred, regions, my.location.file = locs, ...)
}