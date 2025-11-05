if(getRversion() >= "2.15.1") utils::globalVariables(c("UNlocations"))

pop.predict.subnat <- function(end.year = 2060, start.year = 1950, present.year = 2020, wpp.year = 2019,
                                output.dir = file.path(getwd(), "bayesPop.output"),
                               locations = NULL, default.country = NULL, annual = FALSE,
                        inputs = list(
                            popM = NULL, popF = NULL,
                            mxM = NULL, mxF = NULL, srb = NULL,
                            pasfr = NULL, patterns = NULL,
                            migM = NULL, migF = NULL,	
                            migMt = NULL, migFt = NULL, mig = NULL, mig.fdm = NULL,
                            e0F.file = NULL, e0M.file = NULL, tfr.file = NULL, 
                            e0F.sim.dir = NULL, e0M.sim.dir = NULL, 
                            tfr.sim.dir = NULL,
                            migMtraj = NULL, migFtraj = NULL, migtraj = NULL, migFDMtraj = NULL,
                            GQpopM = NULL, GQpopF = NULL, average.annual = NULL
                        ), nr.traj = 1000, keep.vital.events = FALSE,
                        fixed.mx = FALSE, fixed.pasfr = FALSE, lc.for.all = TRUE, 
                        mig.is.rate = FALSE, mig.age.method = c("rc", "fdmp", "fdmnop"),
                        mig.rc.fam = NULL, pasfr.ignore.phase2 = FALSE, 
                        replace.output = FALSE, verbose = TRUE) {
    
    ages <- ages.all(annual, observed = FALSE)
    mig.age.method <- match.arg(mig.age.method)
    
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
                              mig.age.method = mig.age.method, mig.rc.fam = mig.rc.fam,
                              annual = annual, verbose = verbose)

    reg.codes <- intersect(unique(inp$POPm0[,'country_code']), UNcountries())
    
    do.pop.predict(reg.codes, inp, outdir, nr.traj, ages, 
                   keep.vital.events = keep.vital.events, fixed.mx = inp$fixed.mx, fixed.pasfr = fixed.pasfr, 
                   function.inputs = inputs, pasfr.ignore.phase2 = pasfr.ignore.phase2, verbose = verbose)
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
    ds <- bayesTFR:::load.from.wpp(...)
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
                                 annual = FALSE, verbose = FALSE) {
    if(!is.null(inputs$tfr.file)) {
        if(inputs$tfr.file == 'median_')
            TFRpred <- create.trajs.from.wpp(default.country, 
                                             region.codes, 'tfr', wpp.year, 
                                             median.only = TRUE, annual = annual)
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
                                              region.codes, 'tfr', wpp.year, 
                                              annual = annual)
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
                        annual = FALSE, mig.age.method = "rc", mig.rc.fam = NULL, 
                        verbose = FALSE) {
    out <- `in` <- country_code <- NULL
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
    }
    POPm0 <- pop.ini[['M']]
    POPf0 <- pop.ini[['F']]
    
    region.codes <- unique(POPm0$country_code)
    # Get death rates
    MXm.pred <- MXf.pred <- NULL
    if(is.null(inputs$mxM)) {
        if(!is.null(default.country))
            MXm <- create.dataset.from.wpp(default.country, region.codes, 'mxM', wpp.year, annual = annual)
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
          MXf <- create.dataset.from.wpp(default.country, region.codes, 'mxF', wpp.year, annual = annual)
      else stop("mxF must be given if there is no default.country.")
    } else MXf <- swap.reg.code(read.pop.file(inputs$mxF), table.name = "mxF")
    if(fixed.mx) MXf.pred <- MXf[,c('country_code', 'age', proj.periods)]
    MXf <- MXf[,c('country_code', 'age', estim.periods)]
    
    estim.years <- cols.starty[start.index:present.index]
    if(!annual) estim.years <- estim.years + 3
    
    # Get sex ratio at birth
    srb.data <- inputs$srb
    if(is.null(srb.data)) 
        srb.data <- create.dataset.from.wpp(default.country, region.codes, 'sexRatio', wpp.year, annual = annual)
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
        pasfr.data <- create.dataset.from.wpp(default.country, region.codes, 'percentASFR', wpp.year, annual = annual)
    else pasfr.data <- swap.reg.code(read.pop.file(inputs$pasfr), table.name = "pasfr")

    pasfrlist <- .get.pasfr.data(pasfr.data, wpp.year, obs.periods, proj.periods, 
                                 include.projection = fixed.pasfr, annual = annual)
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
    patterns <- .get.mig.mx.pasfr.patterns(inputs, wpp.year, pattern.data = pattern.data,
                                           annual = annual)
    MIGtype <- patterns$mig.type
    MXpattern <- patterns$mx.pattern
    PASFRpattern <- patterns$pasfr.pattern
    MIGshare <- .get.mig.shares(pattern.data, region.codes)
    
    # Get age-specific migration
    mig.rc.inout <- NULL
    is.fdm <- startsWith(mig.age.method, "fdm")
    if(is.fdm){
      if(is.null(inputs[["mig.fdm"]])){
        # create a dataset of model Rogers-Castro
        migio <- data.table(age = ages.all(annual, observed = TRUE),
                            `in` = rcastro.schedule(annual))[, out := `in`]
        mig.rc.inout <- migio[rep(migio[, .I], length(region.codes))]
        mig.rc.inout[, country_code := rep(region.codes, each = nrow(migio))]
        warning("Item 'mig.fdm' is missing for the migration FDM method. Using model Rogers-Castro for both, in- and out-migration.")
      } else {
        mig.rc.inout <- data.table(swap.reg.code(read.pop.file(inputs[["mig.fdm"]])))
        if(! "in" %in% colnames(mig.rc.inout) || ! "out" %in% colnames(mig.rc.inout))
          stop("Column 'in' or 'out' is missing in the mig.fdm dataset.")
      }
      fdmMIGtype.names <- list(MigFDMb0 = "beta0", MigFDMb1 = "beta1", MigFDMmin = "min", 
                               MigFDMsrin = "in_sex_ratio", MigFDMsrout = "out_sex_ratio")
      fdmMIGtype.defaults <- list(beta0 = if(annual) 0.07 else 0.35, beta1 = 0.5, in_sex_ratio = 0.5, out_sex_ratio = 0.5,
                                  min = if(annual) 0.02 else 0.2)
      for(fdmvar in names(fdmMIGtype.names)){
        default.value <- fdmMIGtype.defaults[[fdmMIGtype.names[[fdmvar]]]]
        default.values <- if(fdmvar == "min") pmin(default.value, MIGtype[["MigFDMb0"]]/10) else rep(default.value, nrow(MIGtype))
        if(fdmvar %in% colnames(MIGtype)){
          if(any((fdm.na <- is.na(MIGtype[[fdmvar]]))))  # NAs can appear if running for a different wpp.year than 2024 (due to merging) and there is a mismatch in countries between the two revisions
            MIGtype[[fdmvar]][fdm.na] <- default.values[fdm.na]
        } else MIGtype[[fdmvar]] <- default.values
        mig.rc.inout <- merge(mig.rc.inout, MIGtype[, c("country_code", fdmvar)], by = "country_code")
        setnames(mig.rc.inout, fdmvar, fdmMIGtype.names[[fdmvar]])
      }
    } # end FDM settings
    if(!is.null(mig.rc.fam))  mig.rc.fam <- data.table(mig.rc.fam)
    miginp <- .get.mig.data.subnat(inputs, wpp.year, annual, periods = c(estim.periods, proj.periods), 
                            default.country = default.country, region.codes = region.codes,
                            pop0 = list(M = POPm0, F = POPf0), MIGshare = MIGshare,
                            mig.is.rate = mig.is.rate, mig.age.method = mig.age.method,
                            rc.data = if(is.fdm) mig.rc.inout else mig.rc.fam, 
                            pop = if(is.fdm) pop.ini.matrix else NULL,
                            verbose = verbose)
    
    MIGm <- miginp[["migM"]]
    MIGf <- miginp[["migF"]]
    
    if(!is.null(obs.periods)) {
        avail.obs.periods <- is.element(obs.periods, colnames(MIGm))
        observed$MIGm <- MIGm[,c('country_code', 'age', obs.periods[avail.obs.periods])]
        observed$MIGf <- MIGf[,c('country_code', 'age', obs.periods[avail.obs.periods])]
    }
    MIGm <- MIGm[,c('country_code', 'age', proj.periods[proj.periods %in% colnames(MIGm)])]
    MIGf <- MIGf[,c('country_code', 'age', proj.periods[proj.periods %in% colnames(MIGf)])]

    # assign some migrate-specific attributes, since they get lost by slicing above
    if(!is.null((rates <- attr(miginp[["migM"]], "rate")))){
      if(length(intersect(proj.periods, colnames(rates))) > 0) {
        attr(MIGm, "rate") <- rates[, c('country_code', intersect(proj.periods, colnames(rates))), with = FALSE]
        attr(MIGm, "code") <- attr(miginp[["migM"]], "code")[, c('country_code', intersect(proj.periods, colnames(attr(miginp[["migM"]], "code")))), with = FALSE]
      } else {
        attr(MIGm, "rate") <- NULL
        attr(MIGm, "code") <- NULL
      }
      if(!is.null(obs.periods)) {
        attr(observed$MIGm, "rate") <- rates[, c('country_code', obs.periods[avail.obs.periods]), with = FALSE]
        attr(observed$MIGm, "code") <- attr(miginp[["migM"]], "code")[, c('country_code', obs.periods[avail.obs.periods]), with = FALSE]
      }
    }
    if(!is.null((rcout <- attr(miginp[["migM"]], "rc.out")))){
      attr(MIGm, "rc.out") <- rcout #[, c('country_code', proj.periods), with = FALSE]
      if(!is.null(obs.periods))
        attr(observed$MIGm, "rc.out") <- rcout #[, c('country_code', obs.periods[avail.obs.periods]), with = FALSE]
    }
    if(!is.null((rates <- attr(miginp[["migF"]], "rate")))){
      if(length(intersect(proj.periods, colnames(rates))) > 0) {
        attr(MIGf, "rate") <- rates[, c('country_code', intersect(proj.periods, colnames(rates))), with = FALSE]
        attr(MIGf, "code") <- attr(miginp[["migF"]], "code")[, c('country_code', intersect(proj.periods, attr(miginp[["migF"]], "code"))), with = FALSE]
      } else {
        attr(MIGf, "rate") <- NULL
        attr(MIGf, "code") <- NULL
      }
      if(!is.null(obs.periods)) {
        attr(observed$MIGf, "rate") <- rates[, c('country_code', obs.periods[avail.obs.periods]), with = FALSE]
        attr(observed$MIGf, "code") <- attr(miginp[["migF"]], "code")[, c('country_code', obs.periods[avail.obs.periods]), with = FALSE]
      }
    }
    if(!is.null((rcout <- attr(miginp[["migF"]], "rc.out")))){
      attr(MIGf, "rc.out") <- rcout #[, c('country_code', proj.periods), with = FALSE]
      if(!is.null(obs.periods))
        attr(observed$MIGf, "rc.out") <- rcout #[, c('country_code', obs.periods[avail.obs.periods]), with = FALSE]
    }
    # Get migration trajectories if available
    migpr <- .load.mig.traj(inputs, mig.age.method = mig.age.method, verbose = verbose)
    migMpred <- migpr$M
    migFpred <- migpr$F
    migBpred <- migpr$B
    migFDMpred <- migpr$FDM
    has.mig.traj <- !is.null(migMpred) || !is.null(migFpred) || !is.null(migBpred)
    if(length(mig.is.rate) < 2) mig.is.rate <- rep(mig.is.rate, 2) # one for observed, one for projection
    mig.rate.code <- c(miginp[["migcode"]]*mig.is.rate[1],  
                       (if(has.mig.traj) migpr[["migcode"]] else miginp[["migcode"]])*mig.is.rate[2])
    
    # Get life expectancy
    e0F.wpp.median.loaded <- FALSE
    e0Fpred <- e0Mpred <- NULL
    if(!fixed.mx){
        if(!is.null(inputs$e0F.file)) { # female
            if(inputs$e0F.file == 'median_') {
                e0Fpred <- create.trajs.from.wpp(default.country, 
                                                 region.codes, 'e0F', wpp.year, 
                                                 median.only = TRUE, annual = annual)
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
                                                     median.only = TRUE, annual = annual)
                    e0F.wpp.median.loaded <- TRUE
                } else 
                    e0Fpred <- get.e0.prediction(inputs$e0F.sim.dir, mcmc.dir=NA)
            } else e0Fpred <- create.trajs.from.wpp(default.country, 
                                                    region.codes, 'e0F', wpp.year, annual = annual)			
        }
        
        if(!is.null(inputs$e0M.file)) { # male
            if(inputs$e0M.file == 'median_')
                e0Mpred <- create.trajs.from.wpp(default.country, 
                                                 region.codes, 'e0M', wpp.year, 
                                                 median.only = TRUE, annual = annual)
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
                                                         'e0M', wpp.year, annual = annual)
                    else {
                        if(!has.e0.jmale.prediction(e0Fpred))
                            stop('No joint prediction for female and male available. Correct the e0M.sim.dir argument.' )
                        e0Mpred <- get.e0.jmale.prediction(e0Fpred)
                    }
                } else e0Mpred <- get.e0.prediction(inputs$e0M.sim.dir, mcmc.dir=NA) # independent from female
            } else
                e0Mpred <- create.trajs.from.wpp(default.country, 
                                                 region.codes, 'e0M', wpp.year, annual = annual)
        }
    } # end if(!fixed.mx)
    
    # Get TFR
    TFRpred <- .get.tfr.data.subnat(inputs, wpp.year, default.country, region.codes, 
                                    annual = annual, verbose=verbose)
    
    # Group quarters
    for(sex in c('M', 'F')) {
      dataset.name <- paste0('GQpop', sex)
      if(is.null(inputs[[dataset.name]])) next
      GQ[[sex]] <- .format.gq(inputs[[dataset.name]], annual,  
                              c(estim.years[length(estim.years)], proj.years), # need to include the current year
                              c(obs.periods[length(obs.periods)], proj.periods), 
                              "reg_code", what = dataset.name, verbose)
      colnames(GQ[[sex]])[1] <- 'country_code'
    }
    GQm <- GQ[['M']]
    GQf <- GQ[['F']]

    inp <- new.env()
    for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MXm.pred', 'MXf.pred', 'MXpattern', 'SRB',
                 'PASFR', 'PASFRpattern', 'MIGtype', 'MIGm', 'MIGf', 'GQm', 'GQf',
                 'e0Mpred', 'e0Fpred', 'TFRpred', 'migMpred', 'migFpred', 'migBpred', 'migFDMpred', 'estim.years', 'proj.years', 'wpp.year', 
                 'start.year', 'present.year', 'end.year', 'annual', 'fixed.mx', 'fixed.pasfr', 
                 'lc.for.all', 'mig.rate.code', 'mig.age.method', 'mig.rc.fam', 'mig.rc.inout', 'observed'))
        assign(par, get(par), envir=inp)
    inp$pop.matrix <- list(male=pop.ini.matrix[['M']], female=pop.ini.matrix[['F']])
    env <- new.env()
    do.call("data", list("pasfr_global_norms", envir = env))
    inp$PASFRnorms <- if(annual) env$pasfr.glob.norms1 else env$pasfr.glob.norms5
    if(!is.null(obs.periods)) { 
      # if any of the observed years are missing in the global norm, use the latest norm for those time periods
      missing.years <- obs.periods[! obs.periods %in% colnames(inp$PASFRnorms$PasfrGlobalNorm)]
      if(length(missing.years) > 0) {
        last.norm <- inp$PASFRnorms$PasfrGlobalNorm[, rep(ncol(inp$PASFRnorms$PasfrGlobalNorm), length(missing.years)), drop = FALSE]
        colnames(last.norm) <- missing.years
        inp$PASFRnorms$PasfrGlobalNorm <- cbind(inp$PASFRnorms$PasfrGlobalNorm, last.norm)
      }
    }
    inp$lc.for.hiv <- TRUE
    inp$average.annual <- inputs$average.annual
    return(inp)
}

.get.mig.data.subnat <- function(inputs, wpp.year, annual, periods, default.country, region.codes, 
                                 pop0 = NULL, mig.age.method = "rc",
                                 MIGshare = NULL, mig.is.rate = c(FALSE, FALSE), 
                                 rc.data = NULL, pop = NULL, verbose = FALSE) {
  # Get age-specific migration
  wppds <- data(package=paste0('wpp', wpp.year))
  if(startsWith(mig.age.method, "fdm")){ # need also population 
    popdt <- NULL
    for(sex in c("M", "F")){
      popdt <- rbind(popdt, 
                     cbind(data.table(country_code = as.integer(sapply(strsplit(rownames(pop[[sex]]), "_"), function(x) x[1])),
                              age = sapply(strsplit(rownames(pop[[sex]]), "_"), function(x) x[2]))[, sx := sex],
                            data.table(pop[[sex]]))
                    )
    }
    if(annual) # if age is in annual form, convert to integers 
      popdt$age <- as.integer(popdt$age)
    popdt <- melt(popdt, id.vars = c("country_code", "age", "sx"), variable.name = "year", value.name = "pop", variable.factor = FALSE)
    popdtt <- popdt[, list(pop = sum(pop)), by = c("country_code", "year", "age")] # sum over sexes
    globpop <- popdt[, list(pop = sum(pop)), by = c("year", "age", "sx")]
    globpopt <- globpop[, list(pop = sum(pop)), by = c("year", "age")] # sum over sexes
  }
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
      if(!mig.age.method %in% c("rc")) migcode <- 4 # TODO: this would be wrong if migcode is 2.
      if(mig.age.method == "fdmp") migcode <- 5 # migcode is only used if migration is given as a rate
      migcols <- intersect(colnames(totmig), periods)
      # disaggregate into ages
      migmtx <- migration.totals2age(totmig, ages = migtempl$age[age.index.all(annual, observed = TRUE)],
                                                           annual = annual, time.periods = migcols, 
                                                           scale = if(is.null(inputs[[fname]])) 0.5 else 1, # since the totals are sums over sexes
                                                            method = mig.age.method, mig.is.rate = mig.is.rate[1], 
                                                           template = migtempl, rc.data = rc.data, 
                                     pop = if(mig.is.rate[1]) popdtt else popdt[sx == sex], 
                                     pop.glob = if(mig.is.rate[1]) globpopt else globpop[sx == sex])
      miginp[[inpname]] <- data.frame(migmtx, check.names = FALSE)
      for(attrib in c("rate", "code", "rc.out")) {
        if(!is.null((val <- attr(migmtx, attrib))))
          attr(miginp[[inpname]], attrib) <- val
      }
      next
    }
    # If we get here, migration is not given. Thus, get it from the national values and given shares
    if(annual && wpp.year < 2022) stop("Migration must be given for an annual simulation and wpp.year < 2022.")
    migdsname <- 'migration'
    if(wpp.year >= 2022) migdsname <- paste0(migdsname, if(annual) 1 else 5)
    
    mignat <- load.wpp.dataset.for.country(default.country, migdsname, wpp.year, annual = annual)
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
    if(annual && "100+" %in% migages && ! "100+" %in% pop0[["M"]]$age) # make the name of the open age group consistent between migration and population
      for (sx in c("M", "F")) pop0[[sx]]$age[pop0[[sx]]$age == 100] <- "100+"
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