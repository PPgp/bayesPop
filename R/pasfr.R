compute.pasfr <- function(inputs, present.year, end.year, wpp.year) {
    observed <- list()
    # get SRB data
    srblist <- .get.srb.data.and.time.periods(inputs$srb, present.year, end.year, wpp.year)
    SRB <- srblist$srb
    observed$SRB <- srblist$obs.srb
    proj.years <- srblist$proj.years
    
    # get PASFR data
    pasfrlist <- .get.pasfr.data(inputs$pasfr, wpp.year, srblist$obs.periods, 
                                 srblist$proj.periods)
    PASFR <- pasfrlist$pasfr
    observed$PASFR <- pasfrlist$obs.pasfr
    
    # get TFR
    if(!is.null(inputs$tfr.sim.dir)) {
        TFRpred <- get.tfr.prediction(inputs$tfr.sim.dir, mcmc.dir=NA)
        country.table <- get.countries.table(TFRpred)
    } else {
        TFRpred <- subset(.load.wpp.traj('tfr', wpp.year), trajectory==1)
        locs <- bayesTFR:::load.bdem.dataset('UNlocations', wpp.year)
        country.table <- merge(TFRpred[,'country_code', drop=FALSE], 
                           locs[, c("country_code", "name")], by="country_code")
        colnames(country.table)[1] <- "code"
    }
    # get pasfr patterns
    patterns <- .get.mig.mx.pasfr.patterns(inputs, wpp.year)
    PASFRpattern <- patterns$pasfr.pattern
    
    inp <- new.env()
    for(par in c('SRB', 'PASFR', 'PASFRpattern', 'TFRpred', 'present.year', 'end.year', 'observed'))
        assign(par, get(par), envir=inp)
    inp$PASFRnorms <- compute.pasfr.global.norms(inp)
    
    # iterate over countries
    resdf <- NULL
    for(icountry in 1:nrow(country.table)) {
        # extract country-specific inputs
        country <- country.table[icountry,]
        #cat(paste0("Computing PASFR for index ", icountry, " country ", country$name, "\n"))
        
        #country.input <- bayesPop:::get.country.inputs(country, inputs, nr.traj, country.table$name[icountry])
        if (!is.data.frame(TFRpred)) {
            country.TFRpred <- get.tfr.trajectories(tfr, country$code)
            country.obj <- get.country.object(country$code, tfr$mcmc.set$meta)
            country.medians.TFRpred <- bayesTFR::get.median.from.prediction(tfr, country.obj$index, country.obj$code)[-1]
            country.obs.tfr <- bayesTFR:::get.tfr.reconstructed(tfr$tfr_matrix_reconstructed, tfr$mcmc.set$meta)
            country.obs.TFRpred <- country.obs.tfr[1:if(!is.null(tfr$present.year.index)) tfr$present.year.index else nrow(country.obs.tfr),country.obj$index]
        } else {
            
        }
        #tfr.median <- apply(country.TFRpred, 1, median)
        
        country.PASFRpattern <- bayesPop:::.get.par.from.inputs('PASFRpattern', inputs, country$code)
        
        country.observed <- list()
        country.observed$PASFR <- bayesPop:::.get.par.from.inputs('PASFR', observed, country$code)
        
        country.inputs <- new.env()
        assign('PASFRpattern', get('country.PASFRpattern'), envir=country.inputs)
        assign('observed', get('country.observed'), envir=country.inputs)
        
        # get pasfr for the median TFR
        pasfr <- bayesPop:::kantorova.pasfr(c(country.obs.TFRpred, tfr.median), 
                                            country.inputs, norms=PASFRnorms, proj.years=proj.years, 
                                            tfr.med=tfr.median[nrow(country.TFRpred)]) 
        
        # alternatively iterate over trajectories
        # would need to somehow bind the resulting pasfr objects together
        #for(itraj in 1:nr.traj) {
        #    pasfr <- bayesPop:::kantorova.pasfr(c(country.input$observed$TFRpred, 
        #                               country.input$TFRpred[,itraj]), 
        #                            country.input, norms=inputs$PASFRnorms, proj.years=inputs$proj.years, 
        #                                 tfr.med=tfr.median[nrow(country.input$TFRpred)])                                             
        #}
        
        # as percentage
        pasfr <- round(pasfr*100, 2)
        
        colnames(pasfr) <- period.columns
        
        # add to other results
        resdf <- rbind(resdf, cbind(data.frame(country_code=country$code, name=country$name, age=age.labels), pasfr))
    }
    
    # export
    write.table(resdf, "percentASFR.txt", sep="\t", row.names=FALSE)
}

start <- function() {
	## inputs
	rootDir <- 'H:\\PPP2017'
	TFRDir <- paste(rootDir, '\\TFR',sep='')
	## directory with TFR predictions
	tfdir <- paste(TFRDir,'\\results20170417bigsmall',sep='')
	
	pasfrObservedFile <- paste0(TFRDir, "\\pasfr.txt")
	
	setwd(tfdir)
	
	wpp.year <- 2017
	start.year <- 1950
	present.year <- 2015
	end.year <- 2100

	ages <- seq(15, 45, by=5)
	age.labels <- paste(ages, ages+4, sep="-")

	#load input data for all countries
	#inputs <- bayesPop:::load.inputs(list(tfr.sim.dir=tfdir, pasfr=pasfrObservedFile), start.year, present.year, end.year, wpp.year)
	SRB <- bayesPop:::load.wpp.dataset('sexRatio', wpp.year)
	names.SRB.data <- names(SRB)
	num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.SRB.data) # index of year-columns
	cols.starty <- as.integer(substr(names.SRB.data[num.columns], 1,4))
    cols.endy <- as.integer(substr(names.SRB.data[num.columns], 6,9))
	start.index <- which((cols.starty <= present.year) & (cols.endy > present.year))
	end.index <- which((cols.endy >= end.year) & (cols.starty < end.year))
	proj.years <- cols.starty[start.index:end.index] + 3
	proj.periods <- names.SRB.data[num.columns[start.index:end.index]]
	obs.periods <- names.SRB.data[num.columns[1:(start.index-1)]]

	tfr <- get.tfr.prediction(tfdir, mcmc.dir=NA)
	
	observed <- list()

	PASFRdata <- bayesPop:::read.pop.file(pasfrObservedFile)
	PASFR <- PASFRdata[,c('country_code', 'age', proj.periods)]	
	avail.obs.periods <- is.element(obs.periods, colnames(PASFRdata))
	observed$PASFR <- PASFRdata[,c('country_code', 'age', obs.periods[avail.obs.periods])]	
	
	# ATTENTION, the vwBaseYearXXXX file will indicate, for each country, which Norm to use
	vwBase <- bayesPop:::read.bayesPop.file(paste('vwBaseYear', wpp.year, '.txt', sep=''))
	create.pattern <- function(dataset, columns) {
		pattern <- data.frame(dataset[,'country_code'])
		for(col in columns)
			if(col %in% colnames(dataset))
				pattern <- cbind(pattern, dataset[,col])
		if(ncol(pattern)==1) pattern <- NULL
		else colnames(pattern) <- c('country_code', columns)[1:ncol(pattern)]
		return(pattern)
	}
	PASFRpattern <- create.pattern(vwBase, c("PasfrNorm", paste0("Pasfr", bayesPop:::.remove.all.spaces(levels(vwBase$PasfrNorm)))))
	
	# Compute Global Norms
	inputs <- new.env()
	assign('PASFRpattern', get('PASFRpattern'), envir=inputs)
	assign('observed', get('observed'), envir=inputs)
	PASFRnorms <- bayesPop:::compute.pasfr.global.norms(inputs)
	
	# projection periods and age groups
	period.columns <- paste(proj.years-3, proj.years+2, sep="-")

	# determine available countries
	country.table <- get.countries.table(tfr)

	#nr.traj <- tfr$nr.traj

	# resulting data frame
	resdf <- NULL

	# iterate over countries
	for(icountry in 1:nrow(country.table)) {
		# extract country-specific inputs
		country <- country.table[icountry,]
		cat(paste0("Computing PASFR for index ", icountry, " country ", country$name, "\n"))
		
		#country.input <- bayesPop:::get.country.inputs(country, inputs, nr.traj, country.table$name[icountry])
		country.TFRpred <- get.tfr.trajectories(tfr, country$code)
		country.obj <- get.country.object(country$code, tfr$mcmc.set$meta)
		country.medians.TFRpred <- bayesTFR::get.median.from.prediction(tfr, country.obj$index, country.obj$code)[-1]
		country.obs.tfr <- bayesTFR:::get.tfr.reconstructed(tfr$tfr_matrix_reconstructed, tfr$mcmc.set$meta)
		country.obs.TFRpred <- country.obs.tfr[1:if(!is.null(tfr$present.year.index)) tfr$present.year.index else nrow(country.obs.tfr),country.obj$index]
		
		tfr.median <- apply(country.TFRpred, 1, median)
		
		country.PASFRpattern <- bayesPop:::.get.par.from.inputs('PASFRpattern', inputs, country$code)
		
		country.observed <- list()
		country.observed$PASFR <- bayesPop:::.get.par.from.inputs('PASFR', observed, country$code)
		
		country.inputs <- new.env()
		assign('PASFRpattern', get('country.PASFRpattern'), envir=country.inputs)
		assign('observed', get('country.observed'), envir=country.inputs)
		
		# get pasfr for the median TFR
		pasfr <- bayesPop:::kantorova.pasfr(c(country.obs.TFRpred, tfr.median), 
											country.inputs, norms=PASFRnorms, proj.years=proj.years, 
											tfr.med=tfr.median[nrow(country.TFRpred)]) 

		# alternatively iterate over trajectories
		# would need to somehow bind the resulting pasfr objects together
		#for(itraj in 1:nr.traj) {
		#    pasfr <- bayesPop:::kantorova.pasfr(c(country.input$observed$TFRpred, 
		#                               country.input$TFRpred[,itraj]), 
		#                            country.input, norms=inputs$PASFRnorms, proj.years=inputs$proj.years, 
		#                                 tfr.med=tfr.median[nrow(country.input$TFRpred)])                                             
		#}
		
		# as percentage
		pasfr <- round(pasfr*100, 2)
		
		colnames(pasfr) <- period.columns
		
		# add to other results
		resdf <- rbind(resdf, cbind(data.frame(country_code=country$code, name=country$name, age=age.labels), pasfr))
	}
	
	# export
	write.table(resdf, "percentASFR.txt", sep="\t", row.names=FALSE)
}

start()