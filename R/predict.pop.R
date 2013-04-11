if(getRversion() >= "2.15.1") utils::globalVariables("LOCATIONS")
data('LOCATIONS', package='bayesPop', envir=environment())

pop.predict <- function(end.year=2100, start.year=1950, present.year=2010, wpp.year=2010,
						countries=NULL, output.dir = file.path(getwd(), "bayesPop.output"),
						inputs=list(
							popM=NULL,
							popF=NULL,
							mxM=NULL,
							mxF=NULL,
							srb=NULL,
							pasfr=NULL,
							mig.type=NULL,
							migM=NULL,
							migF=NULL,	
							e0F.file=NULL, e0M.file=NULL, 
							tfr.file=NULL,
							e0F.sim.dir=NULL, e0M.sim.dir=NULL, 
							tfr.sim.dir=NULL	
						), nr.traj = 1000, keep.vital.events=FALSE,
						replace.output=FALSE, 
						verbose=TRUE) {
	prediction.exist <- FALSE
	ages=seq(0, by=5, length=27)
	if(is.null(countries)) inp <- load.inputs(inputs, start.year, present.year, end.year, wpp.year)
	else {
		if(has.pop.prediction(output.dir) && !replace.output) {
			pred <- get.pop.prediction(output.dir)
			inp <- pred$inputs
			nr.traj <- pred$nr.traj
			ages <- pred$ages
			prediction.exist <- TRUE
		} else inp <- load.inputs(inputs, start.year, present.year, end.year, wpp.year)
	}

	outdir <- file.path(output.dir, 'predictions')
	if(!prediction.exist) {
		if(!replace.output && has.pop.prediction(sim.dir=output.dir))
			stop('Prediction in ', outdir,
			' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		unlink(outdir, recursive=TRUE)
		.remove.cache.file(output.dir)
	} else pop.cleanup.cache(pred)
	
	#data('LOCATIONS', package='bayesPop')
	if(!is.null(countries) && is.na(countries[1])) { # all countries that are not included in the existing prediction
		all.countries <- unique(inp$POPm0[,'country_code'])
		country.codes <- if(!prediction.exist) all.countries
						else all.countries[!is.element(all.countries, pred$countries[,'code'])]
	} else {
		if(!is.null(countries)) {
			if (is.character(countries)) { # at least one of the codes is a character string
				for (icountry in 1:length(countries)) {
					if (is.character(countries[icountry])) {
						country.idx <- which(LOCATIONS[,'name'] == countries[icountry])
						if(length(country.idx) > 0)
							countries[icountry] <- LOCATIONS[country.idx,'country_code']
					}
				}
			}
			country.codes <- as.integer(countries)
		} else
			country.codes <- unique(inp$POPm0[,'country_code'])
	}
	do.pop.predict(country.codes, inp, outdir, nr.traj, ages, pred=if(prediction.exist) pred else NULL,
					keep.vital.events=keep.vital.events, verbose=verbose)
	invisible(get.pop.prediction(output.dir))
}

do.pop.predict <- function(country.codes, inp, outdir, nr.traj, ages, pred=NULL, keep.vital.events=FALSE, verbose=FALSE) {
	not.valid.countries.idx <- c()
	countries.idx <- rep(NA, length(country.codes))

	for(icountry in 1:length(country.codes)) {
		country.idx <- which(LOCATIONS[,'country_code'] == country.codes[icountry])
		if(length(country.idx) == 0) {
			not.valid.countries.idx <- c(not.valid.countries.idx, icountry)
			next
		}
		countries.idx[icountry] <- country.idx
	}

	if(length(not.valid.countries.idx) > 0) {
		warning('Countries ', paste(country.codes[not.valid.countries.idx], collapse=', '), 
					' not found in the LOCATIONS dataset.')
		country.codes <- country.codes[-not.valid.countries.idx]
		countries.idx <- countries.idx[-not.valid.countries.idx]
	}
	ncountries <- length(country.codes)
	nr_project <- length(inp$proj.years)
	nages <- length(ages)
	nest <- length(inp$estim.years)
	if(!file.exists(outdir)) 
		dir.create(outdir, recursive=TRUE)
	present.and.proj.years <- c(inp$estim.years[nest], inp$proj.years)
	prediction.file <- file.path(outdir, 'prediction.rda')	
	quantiles.to.keep <- get.quantiles.to.keep()
	nquant <- length(quantiles.to.keep)
	PIs_cqp <- quantM <- quantF <- array(NA, c(ncountries, nquant, nr_project+1),
						dimnames=list(country.codes, quantiles.to.keep, present.and.proj.years))
	quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(ncountries, nages, nquant, nr_project+1),
						dimnames=list(country.codes, ages, quantiles.to.keep, present.and.proj.years))
	mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(ncountries, 2, nr_project+1), 
						dimnames=list(country.codes, c('mean', 'sd'), present.and.proj.years))
	mx.ages <- c(0,1,ages[2:nages])
	status.for.gui <- paste('out of', ncountries, 'countries.')
	gui.options <- list()
	for(cidx in 1:ncountries) {
		if(getOption('bDem.Poppred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			gui.options$bDem.Poppred.status <- paste('finished', cidx, status.for.gui)
			unblock.gtk('bDem.Poppred', gui.options)
		}

		country <- country.codes[cidx]
		country.idx <- countries.idx[cidx]
		if(verbose)
			cat('\nProcessing country ', country, ' -- ', as.character(LOCATIONS[country.idx,'name']))
		# Extract the country-specific stuff from the inputs
		inpc <- get.country.inputs(country, inp, nr.traj, LOCATIONS[country.idx,'name'])
		if(is.null(inpc)) next
		nr.traj <- min(ncol(inpc$TFRpred), nr.traj)		
		if(verbose)
			cat(' (', nr.traj, ' trajectories )')
		
		npred <- min(nrow(inpc$TFRpred), nr_project)
		npredplus1 <- npred+1
		totp <- matrix(NA, nrow=npredplus1, ncol=nr.traj, 
					dimnames=list(present.and.proj.years, NULL))
		totpm <- totpf <- array(NA, dim=c(27, npredplus1, nr.traj), 
							dimnames=list(ages, present.and.proj.years, NULL))
		nvariants <- nrow(inpc$TFRhalfchild)
		totp.hch <- matrix(NA, nrow=npredplus1, ncol=nvariants,
					dimnames=list(present.and.proj.years, NULL))
		totpm.hch <- totpf.hch <- array(NA, dim=c(27, npredplus1, nvariants), 
					dimnames=list(ages, present.and.proj.years, NULL))
		if(keep.vital.events) {
			btm <- btf <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm <- deathsf <- array(0, dim=c(27, npredplus1, nr.traj), dimnames=list(ages, present.and.proj.years, NULL))
			asfert <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			btm.hch <- btf.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years, NULL))
			asfert.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			mxm <- mxf <- array(0, dim=c(28, npredplus1, nr.traj), dimnames=list(mx.ages, present.and.proj.years, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(28, npredplus1, nvariants), dimnames=list(mx.ages, present.and.proj.years, NULL))
		}
		debug <- FALSE
		MxKan <- runKannisto(nest, inpc)
		npasfr <- nrow(inpc$PASFR)
		if(keep.vital.events) observed <- compute.observedVE(inpc, inp$pop.matrix, inpc$MIGtype, MxKan, country, inp$estim.years)
		for(itraj in 1:nr.traj) {
			if(any(is.na(inpc$TFRpred[,itraj]))) next
			asfr <- inpc$PASFR/100.
			for(i in 1:npasfr) asfr[i,] <- inpc$TFRpred[,itraj] * asfr[i,]
			LTres <- modifiedLC(npred, MxKan, inpc$e0Mpred[,itraj], 
									inpc$e0Fpred[,itraj], verbose=verbose, debug=debug)
			popres <- StoPopProj(npred, inpc, LTres, asfr, inpc$MIGtype, LOCATIONS[country.idx,'name'],
									keep.vital.events=keep.vital.events)
			totp[,itraj] <- popres$totpop
			totpm[,,itraj] <- popres$mpop
			totpf[,,itraj] <- popres$fpop
			if(keep.vital.events) {
				btm[,2:npredplus1,itraj] <- popres$mbt
				btm[1:dim(observed$btm)[1],1,itraj] <- observed$btm[,dim(observed$btm)[2],]
				btf[,2:npredplus1,itraj] <- popres$fbt
				btf[1:dim(observed$btf)[1],1,itraj] <- observed$btf[,dim(observed$btf)[2],]
				deathsm[,2:npredplus1,itraj] <- popres$mdeaths
				deathsm[1:dim(observed$deathsm)[1],1,itraj] <- observed$deathsm[,dim(observed$deathsm)[2],]
				deathsf[,2:npredplus1,itraj] <- popres$fdeaths
				deathsf[1:dim(observed$deathsf)[1],1,itraj] <- observed$deathsf[,dim(observed$deathsf)[2],]
				asfert[,2:npredplus1,itraj] <- asfr
				asfert[1:dim(observed$asfert)[1],1,itraj] <- observed$asfert[,dim(observed$asfert)[2],]
				mxm[,2:npredplus1,itraj] <- LTres$mx[[1]]
				mxm[1:dim(MxKan[[1]]$mx)[1],1,itraj] <- MxKan[[1]]$mx[,dim(MxKan[[1]]$mx)[2]]
				mxf[,2:npredplus1,itraj] <- LTres$mx[[2]]
				mxf[1:dim(MxKan[[2]]$mx)[1],1,itraj] <- MxKan[[2]]$mx[,dim(MxKan[[2]]$mx)[2]]
			}
		}
		for (variant in 1:nvariants) { # compute the two half child variants
			asfr <- inpc$PASFR/100.
			for(i in 1:npasfr) asfr[i,] <- inpc$TFRhalfchild[variant,] * asfr[i,]
			LTres <- modifiedLC(npred, MxKan, inpc$e0Mmedian, 
									inpc$e0Fmedian, verbose=verbose, debug=debug)
			popres <- StoPopProj(npred, inpc, LTres, asfr, inpc$MIGtype, LOCATIONS[country.idx,'name'], 
							keep.vital.events=keep.vital.events)
			totp.hch[,variant] <- popres$totpop
			totpm.hch[,,variant] <- popres$mpop
			totpf.hch[,,variant] <- popres$fpop
			if(keep.vital.events) {
				btm.hch[,2:npredplus1,variant] <- popres$mbt
				btm.hch[,1,variant] <- btm[,1,1]
				btf.hch[,2:npredplus1,variant] <- popres$fbt
				btf.hch[,1,variant] <- btf[,1,1]
				deathsm.hch[,2:npredplus1,variant] <- popres$mdeaths
				deathsm.hch[,1,variant] <- deathsm[,1,1]
				deathsf.hch[,2:npredplus1,variant] <- popres$fdeaths
				deathsf.hch[,1,variant] <- deathsf[,1,1]
				asfert.hch[,2:npredplus1,variant] <- asfr
				asfert.hch[,1,variant] <- asfert[,1,1]
				mxm.hch[,2:npredplus1,variant] <- LTres$mx[[1]]
				mxf.hch[,2:npredplus1,variant] <- LTres$mx[[2]]
				mxm.hch[,1,variant] <- mxm[,1,1]
				mxf.hch[,1,variant] <- mxf[,1,1]
			}
		}
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch,
			 file = file.path(outdir, paste('totpop_country', country, '.rda', sep='')))
		if(keep.vital.events) 
			save(btm, btf, deathsm, deathsf, asfert, mxm, mxf,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, 
				mxm.hch, mxf.hch, 
				observed,
					file=file.path(outdir, paste('vital_events_country', country, '.rda', sep='')))
		PIs_cqp[cidx,,] = apply(totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sd[cidx,1,] <- apply(totp, 1, mean, na.rm = TRUE)
		mean_sd[cidx,2,] = apply(totp, 1, sd, na.rm = TRUE)
		for (i in 1:nages) {
			if(nr.traj == 1) {
				quantMage[cidx,i,,] <- matrix(rep(totpm[i,,1],nquant) , nrow=nquant, byrow=TRUE)
				quantFage[cidx,i,,] <- matrix(rep(totpf[i,,1],nquant) , nrow=nquant, byrow=TRUE)
				quantPropMage[cidx,i,,] <- matrix(rep(totpm[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
				quantPropFage[cidx,i,,] <- matrix(rep(totpf[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
			} else {
				quantMage[cidx,i,,] <- apply(totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantFage[cidx,i,,] <- apply(totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropMage[cidx,i,,] <- apply(totpm[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropFage[cidx,i,,] <- apply(totpf[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			}
		}
		stotpm <- colSums(totpm, na.rm=TRUE)
		quantM[cidx,,] = apply(stotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdM[cidx,1,] <- apply(stotpm, 1, mean, na.rm = TRUE)
		mean_sdM[cidx,2,] = apply(stotpm, 1, sd, na.rm = TRUE)
		stotpf <- colSums(totpf, na.rm=TRUE)
		quantF[cidx,,] = apply(stotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdF[cidx,1,] <- apply(stotpf, 1, mean, na.rm = TRUE)
		mean_sdF[cidx,2,] = apply(stotpf, 1, sd, na.rm = TRUE)
		
		#save updated meta file
		country.row <- LOCATIONS[country.idx,c('country_code', 'name')]
		colnames(country.row) <- c('code', 'name')
		if(!exists('bayesPop.prediction')) { # first pass
			bayesPop.prediction <- if(!is.null(pred)) pred 
					else structure(list(output.directory=outdir,
							nr.traj = nr.traj,	
							# assign empty arrays
							quantiles = PIs_cqp[c(),,,drop=FALSE],
               				traj.mean.sd = mean_sd[c(),,,drop=FALSE],
               				quantilesM = quantM[c(),,,drop=FALSE], 
               				traj.mean.sdM = mean_sdM[c(),,,drop=FALSE],
               				quantilesF = quantF[c(),,,drop=FALSE], 
               				traj.mean.sdF = mean_sdF[c(),,,drop=FALSE],
               				quantilesMage = quantMage[c(),,,,drop=FALSE], 
               				quantilesFage = quantFage[c(),,,,drop=FALSE], 
               				quantilesPropMage = quantPropMage[c(),,,,drop=FALSE], 
               				quantilesPropFage = quantPropFage[c(),,,,drop=FALSE],
               				estim.years=inp$estim.years, 
               				proj.years=c(inp$estim.years[length(inp$estim.years)], 
               									inp$proj.years), # includes present period
			   				inputs = inp,
			   				countries=as.data.frame(matrix(NA, nrow=0, ncol=2, dimnames=list(NULL, c('code', 'name')))),
			   				ages=ages), class='bayesPop.prediction')
		}	
		idx.in.pred.overwrite <- which(bayesPop.prediction$countries[,'code'] == country)
		if(length(idx.in.pred.overwrite)>0) {
			bayesPop.prediction$quantiles[idx.in.pred.overwrite,,] <- PIs_cqp[cidx,,,drop=FALSE]
			bayesPop.prediction$traj.mean.sd[idx.in.pred.overwrite,,] <- mean_sd[cidx,,,drop=FALSE]
			bayesPop.prediction$traj.mean.sdM[idx.in.pred.overwrite,,] <- mean_sdM[cidx,,,drop=FALSE]
			bayesPop.prediction$traj.mean.sdF[idx.in.pred.overwrite,,] <- mean_sdF[cidx,,,drop=FALSE]
			bayesPop.prediction$quantilesM[idx.in.pred.overwrite,,] <- quantM[cidx,,,drop=FALSE]
			bayesPop.prediction$quantilesF[idx.in.pred.overwrite,,] <- quantF[cidx,,,drop=FALSE]
			bayesPop.prediction$quantilesMage[idx.in.pred.overwrite,,,] <- quantMage[cidx,,,,drop=FALSE]
			bayesPop.prediction$quantilesFage[idx.in.pred.overwrite,,,] <- quantFage[cidx,,,,drop=FALSE]
			bayesPop.prediction$quantilesPropMage[idx.in.pred.overwrite,,,] <- quantPropMage[cidx,,,,drop=FALSE]
			bayesPop.prediction$quantilesPropFage[idx.in.pred.overwrite,,,] <- quantPropFage[cidx,,,,drop=FALSE]
		} else { 
			bayesPop.prediction$quantiles <- abind(bayesPop.prediction$quantiles, 
												PIs_cqp[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$traj.mean.sd <- abind(bayesPop.prediction$traj.mean.sd, 
												mean_sd[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$traj.mean.sdM <- abind(bayesPop.prediction$traj.mean.sdM, 
												mean_sdM[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$traj.mean.sdF <- abind(bayesPop.prediction$traj.mean.sdF, 
												mean_sdF[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesM <- abind(bayesPop.prediction$quantilesM, 
												quantM[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesF <- abind(bayesPop.prediction$quantilesF, 
												quantF[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesMage <- abind(bayesPop.prediction$quantilesMage, 
												quantMage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesFage <- abind(bayesPop.prediction$quantilesFage, 
												quantFage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesPropMage <- abind(bayesPop.prediction$quantilesPropMage, 
												quantPropMage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesPropFage <- abind(bayesPop.prediction$quantilesPropFage, 
												quantPropFage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$countries <- rbind(bayesPop.prediction$countries, country.row)
		}
		save(bayesPop.prediction, file=prediction.file)
	} 
	cat('\nPrediction stored into', outdir, '\n')
	return(bayesPop.prediction)
}

read.pop.file <- function(file) 
	return(read.delim(file=file, comment.char='#', check.names=FALSE))
	
read.bayesPop.file <- function(file)
	return(get(do.call(data, list(strsplit(file, '.', fixed=TRUE)[[1]][-2]))))
	#return(read.pop.file(file.path(find.package("bayesPop"), "data", file)))

load.inputs <- function(inputs, start.year, present.year, end.year, wpp.year) {
	observed <- list()
	# Get initial population counts
	if(is.null(inputs$popM)) 
		POPm0 <- read.bayesPop.file(paste('PopByAgeMale', wpp.year, '.txt', sep=''))	
	else POPm0 <- read.pop.file(inputs$popM)
	num.columns <- grep('^[0-9]{4}$', colnames(POPm0), value=TRUE) # values of year-columns
	if(!is.element(as.character(present.year), colnames(POPm0))) {
		stop('Wrong present.year. ', present.year, ' not available in the popM file.\nAvailable years: ',
				num.columns)
	}
	popm.matrix <- POPm0[,num.columns]
	dimnames(popm.matrix) <- list(paste(POPm0[,'country_code'], POPm0[,'age'], sep='_'), 
									as.character(as.integer(num.columns)-2))
	POPm0 <- POPm0[,c('country_code', 'age', present.year)]
	
	if(is.null(inputs$popF)) 
		POPf0 <- read.bayesPop.file(paste('PopByAgeFemale', wpp.year, '.txt', sep=''))
	else POPf0 <- read.pop.file(inputs$popF)
	num.columns <- grep('^[0-9]{4}$', colnames(POPf0), value=TRUE) # values of year-columns
	if(!is.element(as.character(present.year), colnames(POPf0))) {
		stop('Wrong present.year. ', present.year, ' not available in the popF file.\nAvailable years: ',
				num.columns)
	}
	popf.matrix <- POPf0[, num.columns]
	dimnames(popf.matrix) <- list(paste(POPf0[,'country_code'], POPf0[,'age'], sep='_'), 
									as.character(as.integer(num.columns)-2))
	POPf0 <- POPf0[,c('country_code', 'age', present.year)]
	# Get death rates
	if(is.null(inputs$mxM)) 
		MXm <- read.bayesPop.file(paste('ASMRMale', wpp.year, '.txt', sep=''))
	else MXm <- read.pop.file(inputs$mxM)
	names.MXm.data <- names(MXm)
	num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.MXm.data) # index of year-columns
	cols.starty <- as.integer(substr(names.MXm.data[num.columns], 1,4))
    cols.endy <- as.integer(substr(names.MXm.data[num.columns], 6,9))
	start.index <- which((cols.starty <= start.year) & (cols.endy > start.year))
	present.index <- which((cols.endy >= present.year) & (cols.starty < present.year))
	estim.periods <- names.MXm.data[num.columns[start.index:present.index]]
	MXm <- MXm[,c('country_code', 'age', estim.periods)]
	mid.est.years <- cols.starty[start.index:present.index] + 3
	if(is.null(inputs$mxF)) 
		MXf <- read.bayesPop.file(paste('ASMRFemale', wpp.year, '.txt', sep=''))
	else MXf <- read.pop.file(inputs$mxF)
	MXf <- MXf[,c('country_code', 'age', estim.periods)]
	
	# Get sex ratio at birth
	if(is.null(inputs$srb)) 
		SRB <- read.bayesPop.file(paste('SexRatioAtBirth', wpp.year, '.txt', sep=''))
	else SRB <- read.pop.file(inputs$srb)
	names.SRB.data <- names(SRB)
	num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.SRB.data) # index of year-columns
	cols.starty <- as.integer(substr(names.SRB.data[num.columns], 1,4))
    cols.endy <- as.integer(substr(names.SRB.data[num.columns], 6,9))
	start.index <- which((cols.starty <= present.year) & (cols.endy > present.year))
	end.index <- which((cols.endy >= end.year) & (cols.starty < end.year))
	if(length(end.index) == 0) {
		end.index <- length(num.columns)
		warning('Data for SexRatioAtBirth not available for all projection periods.\nLast projection period set to ', 
					names.SRB.data[num.columns[end.index]])
	}
	proj.periods <- names.SRB.data[num.columns[start.index:end.index]]
	obs.periods <- NULL
	if(start.index > 1) {
		obs.periods <- names.SRB.data[num.columns[1:(start.index-1)]]
		observed$SRB <- SRB[,c('country_code', obs.periods)]
	}
	SRB <- SRB[,c('country_code', proj.periods)]
	mid.proj.years <- cols.starty[start.index:end.index] + 3
	
	# Get percentage age-specific fertility rate
	if(is.null(inputs$pasfr)) 
		PASFR <- read.bayesPop.file(paste('PercentASFR', wpp.year, '.txt', sep=''))
	else PASFR <- read.pop.file(inputs$pasfr)
	if(!is.null(obs.periods)) {
		avail.obs.periods <- is.element(obs.periods, colnames(PASFR))
		observed$PASFR <- PASFR[,c('country_code', 'age', obs.periods[avail.obs.periods])]
	}
	PASFR <- PASFR[,c('country_code', 'age', proj.periods)]
	
	# Get migration type and base year
	if(is.null(inputs$mig.type)) 
		MIGtype <- read.bayesPop.file(paste('vwBaseYear', wpp.year, '.txt', sep=''))
	else MIGtype <- read.pop.file(inputs$mig.type)
	MIGtype <- MIGtype[,c('country_code', 'ProjFirstYear', 'MigCode')]
	# Get age-specific migration
	if(is.null(inputs$migM)) 
		MIGm <- read.bayesPop.file(paste('MigByAgeMale', wpp.year, '.txt', sep=''))
	else MIGm <- read.pop.file(inputs$migM)
	if(is.null(inputs$migF))
		MIGf <- read.bayesPop.file(paste('MigByAgeFemale', wpp.year, '.txt', sep=''))
	else MIGf <- read.pop.file(inputs$migF)
	if(!is.null(obs.periods)) {
		avail.obs.periods <- is.element(obs.periods, colnames(MIGm))
		observed$MIGm <- MIGm[,c('country_code', 'age', obs.periods[avail.obs.periods])]
		observed$MIGf <- MIGf[,c('country_code', 'age', obs.periods[avail.obs.periods])]
	}
	MIGm <- MIGm[,c('country_code', 'age', proj.periods)]
	MIGf <- MIGf[,c('country_code', 'age', proj.periods)]
	
	# Get life expectancy
	if(!is.null(inputs$e0F.sim.dir))  # female
		e0Fpred <- get.e0.prediction(inputs$e0F.sim.dir, mcmc.dir=NA)
	else {
		file.name <- if(!is.null(inputs$e0F.file)) inputs$e0F.file
					else file.path(find.package("bayesPop"), "ex-data", 
							paste('e0Fwpp', wpp.year, '.csv', sep=''))
		if(!file.exists(file.name))
			stop('File ', file.name, 
				' does not exist.\nSet e0F.sim.dir, e0F.file or change WPP year.')
		e0Fpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
		e0Fpred <- e0Fpred[,c('LocID', 'Year', 'e0')]
		colnames(e0Fpred) <- c('country_code', 'year', 'value')
		#observed$e0Fpred <- e0Fpred[e0Fpred[,'year']+2<present.year,]
		#e0Fpred <- e0Fpred[e0Fpred[,'year']+2>=present.year,]
	} 

	if(!is.null(inputs$e0M.sim.dir)) { # male
		if(inputs$e0M.sim.dir == 'joint_') {
			if(!has.e0.jmale.prediction(e0Fpred))
				stop('No joint prediction for female and male available. Correct the e0M.sim.dir argument.' )
			e0Mpred <- get.e0.jmale.prediction(e0Fpred)
		} else e0Mpred <- get.e0.prediction(inputs$e0M.sim.dir, mcmc.dir=NA)
	} else {
		file.name <- if(!is.null(inputs$e0M.file)) inputs$e0M.file 
					else file.path(find.package("bayesPop"), "ex-data", 
							paste('e0Mwpp', wpp.year, '.csv', sep=''))
		if(!file.exists(file.name))
			stop('File ', file.name, 
				' does not exist.\nSet e0M.sim.dir, e0M.file or change WPP year.')
		e0Mpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
		e0Mpred <- e0Mpred[,c('LocID', 'Year', 'e0')]
		colnames(e0Mpred) <- c('country_code', 'year', 'value')
		#observed$e0Mpred <- e0Mpred[e0Mpred[,'year']+2<present.year,]
		#e0Mpred <- e0Mpred[e0Mpred[,'year']+2>=present.year,]

	} 
		
	# Get TFR
	if(!is.null(inputs$tfr.sim.dir)) 
		TFRpred <- get.tfr.prediction(inputs$tfr.sim.dir, mcmc.dir=NA)
	else {
		file.name <- if(!is.null(inputs$tfr.file)) inputs$tfr.file
					else file.path(find.package("bayesPop"), "ex-data", 
							paste('TFRwpp', wpp.year, '.csv', sep=''))
		if(!file.exists(file.name))
			stop('File ', file.name, 
				' does not exist.\nSet tfr.sim.dir, tfr.file or change WPP year.')
		TFRpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
		TFRpred <- TFRpred[,c('LocID', 'Year', 'TF')]
		colnames(TFRpred) <- c('country_code', 'year', 'value')
		#observed$TFRpred <- TFRpred[TFRpred[,'year']+2<present.year,]
		#TFRpred <- TFRpred[TFRpred[,'year']+2>=present.year,]
	} 
	return(list(POPm0=POPm0, POPf0=POPf0, MXm=MXm, MXf=MXf, SRB=SRB,
				PASFR=PASFR, MIGtype=MIGtype, MIGm=MIGm, MIGf=MIGf,
				e0Mpred=e0Mpred, e0Fpred=e0Fpred, TFRpred=TFRpred, 
				estim.years=mid.est.years, proj.years=mid.proj.years,
				pop.matrix=list(male=popm.matrix, female=popf.matrix),
				observed=observed))
}

.get.par.from.inputs <- function(par, inputs, country) {
	if(is.null(inputs[[par]])) return (NULL)
	idx <- inputs[[par]][,'country_code'] == country
	res <- inputs[[par]][idx,,drop=FALSE]
	return (as.matrix(res[, !is.element(colnames(res), c('country_code', 'age')),drop=FALSE]))
}

get.country.inputs <- function(country, inputs, nr.traj, country.name) {
	inpc <- list()
	obs <- list()
	for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'SRB',
						'PASFR', 'MIGtype', 'MIGm', 'MIGf')) {
		inpc[[par]] <- .get.par.from.inputs(par, inputs, country)
		obs[[par]] <- .get.par.from.inputs(par, inputs$observed, country)
	}
	inpc[['MIGBaseYear']] <- inpc[['MIGtype']][,'ProjFirstYear']
	inpc[['MIGtype']] <- inpc[['MIGtype']][,'MigCode']
	what.traj <- list(TFRpred='TFR', e0Mpred='male e0', e0Fpred='female e0')
	medians <- list()
	for(par in names(what.traj)) {
		if (is.data.frame(inputs[[par]])) {
			cidx <- inputs[[par]][,'country_code'] == country 
			idx <- cidx & is.element(inputs[[par]][,'year'], inputs$proj.years)
			if(sum(idx) == 0) {
				warning('No ', what.traj[[par]], ' trajectories for ', country.name, 
					'. No population projection generated.')
				return(NULL)
			}
			inpc[[par]] <- inputs[[par]][idx,]
			inpc[[par]] <- matrix(inpc[[par]][, 'value'], nrow=length(inputs$proj.years))
			rownames(inpc[[par]]) <- as.character(inputs$proj.years)
			medians[[par]] <- apply(inpc[[par]], 1, quantile, 0.5, na.rm = TRUE)
			obsidx <- cidx & inputs[[par]][,'year'] < min(inputs$proj.years)
			if(sum(obsidx) > 0) {
				obs[[par]] <- inputs[[par]][obsidx,]
				obs.years <- unique(obs[[par]][, 'year'])
				obs[[par]] <- matrix(obs[[par]][, 'value'], nrow=length(obs.years))[,1]
				names(obs[[par]]) <- as.character(obs.years)
			}
		}
	}
	if(is.null(inpc$TFRpred)) {
		inpc$TFRpred <- get.tfr.trajectories(inputs$TFRpred, country)
		if(is.null(inpc$TFRpred)) {
			warning('No TFR trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)	
		}
		country.obj <- get.country.object(country, inputs$TFRpred$mcmc.set$meta)
		medians$TFRpred <- bayesTFR:::get.median.from.prediction(inputs$TFRpred, country.obj$index, country.obj$code)[-1]
		obs.tfr <- bayesTFR:::get.tfr.reconstructed(inputs$TFRpred$tfr_matrix_reconstructed, inputs$TFRpred$mcmc.set$meta)
		obs$TFRpred <- obs.tfr[1:if(!is.null(inputs$TFRpred$present.year.index)) inputs$TFRpred$present.year.index else nrow(obs.tfr),country.obj$index]
	} 
	inpc$TFRhalfchild <- bayesTFR:::get.half.child.variant(median=medians$TFRpred, increment=c(0.25, 0.4, 0.5))
	if(!all(is.element(inputs$proj.years, colnames(inpc$TFRhalfchild))))
		stop('Mismatch in projection periods of TFR and the target projection years.')
	inpc$TFRhalfchild <- inpc$TFRhalfchild[,as.character(inputs$proj.years)]
	
	if(is.null(inpc$e0Mpred)) {
		inpc$e0Mpred <- get.e0.trajectories(inputs$e0Mpred, country)
		if(is.null(inpc$e0Mpred)) {
			warning('No male e0 trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)	
		}
		country.obj <- get.country.object(country, inputs$e0Mpred$mcmc.set$meta)
		medians$e0Mpred <- bayesTFR:::get.median.from.prediction(inputs$e0Mpred, country.obj$index, country.obj$code)
		obs.e0M <- bayesLife:::get.e0.reconstructed(inputs$e0Mpred$e0.matrix.reconstructed, inputs$e0Mpred$mcmc.set$meta)
		obs$e0Mpred <- obs.e0M[1:if(!is.null(inputs$e0Mpred$present.year.index)) inputs$e0Mpred$present.year.index else nrow(obs.e0M),country.obj$index]
	}
	if(!all(is.element(inputs$proj.years, names(medians$e0Mpred))))
		stop('Mismatch in projection periods of male e0 and the target projection years.')
	inpc$e0Mmedian <- medians$e0Mpred[as.character(inputs$proj.years)]
	
	if(is.null(inpc$e0Fpred)) {
		inpc$e0Fpred <- get.e0.trajectories(inputs$e0Fpred, country)
		if(is.null(inpc$e0Fpred)) {
			warning('No female e0 trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)
		}
		country.obj <- get.country.object(country, inputs$e0Fpred$mcmc.set$meta)
		medians$e0Fpred <- bayesTFR:::get.median.from.prediction(inputs$e0Fpred, country.obj$index, country.obj$code)
		obs.e0F <- bayesLife:::get.e0.reconstructed(inputs$e0Fpred$e0.matrix.reconstructed, inputs$e0Fpred$mcmc.set$meta)
		obs$e0Fpred <- obs.e0F[1:if(!is.null(inputs$e0Fpred$present.year.index)) inputs$e0Fpred$present.year.index else nrow(obs.e0F),country.obj$index]

	}
	if(!all(is.element(inputs$proj.years, names(medians$e0Fpred))))
		stop('Mismatch in projection periods of female e0 and the target projection years.')
	inpc$e0Fmedian <- medians$e0Fpred[as.character(inputs$proj.years)]
	
	indices <- list()
	indices$e0Mpred <- get.traj.index(nr.traj, inpc$e0Mpred)
	indices$e0Fpred <- if (ncol(inpc$e0Mpred) != ncol(inpc$e0Fpred)) get.traj.index(nr.traj, inpc$e0Fpred)
					else indices$e0Mpred
	indices$TFRpred <- get.traj.index(nr.traj, inpc$TFRpred)
	nr.traj <- min(length(indices$e0Mpred), length(indices$e0Fpred), length(indices$TFRpred))
	for (par in c('TFRpred', 'e0Mpred', 'e0Fpred')) {
		if(length(indices[[par]]) != nr.traj)
			indices[[par]] <- get.traj.index(nr.traj, inpc[[par]])
	}
			
	if(length(indices$TFRpred) != length(indices$e0Mpred) || length(indices$e0Mpred) != length(indices$e0Fpred)) {
		warning('Mismatch in number of trajectories of TFR and/or life expactancy for ', 
						country.name, '. No population projection generated.')
		return(NULL)	
	}
	for (par in c('TFRpred', 'e0Mpred', 'e0Fpred')) {
		inpc[[par]] <- inpc[[par]][,indices[[par]], drop=FALSE]
		inpc[[par]] <- inpc[[par]][as.character(inputs$proj.years),, drop=FALSE]
	}
	inpc$observed <- obs
	return(inpc)
}

get.traj.index <- function(nr.traj, traj) {
	ncoltraj <- ncol(traj)
	if(nr.traj >= ncoltraj) return(1:ncoltraj)
	return(seq(1, ncoltraj, length=nr.traj))
}

runKannisto <- function(nest, inputs) {
	# extend mx, get LC ax,bx,k1
	#mxFKan <- c(KannistoAxBx(nest, inputs$MXf, inputs$MIGBaseYear), sex=2)
	#mxMKan <- c(KannistoAxBx(nest, inputs$MXm, inputs$MIGBaseYear), sex=1)
	Kan <- KannistoAxBx.joint(nest, inputs$MXm, inputs$MXf, inputs$MIGBaseYear)
	mxMKan <- c(Kan$male, sex=1)
	mxFKan <- c(Kan$female, sex=2)
	#stop('')
	#mxFKan$ax.orig <- mxFKan$ax
	#mxFKan$ax[(mxMKan$ax - mxFKan$ax)<0.2] <- mxMKan$ax[(mxMKan$ax - mxFKan$ax)<0.2]-0.2
	#mxFKan$bx.orig <- mxFKan$bx
	#mxFKan$bx <- mxMKan$bx
	bx <- 0.5 * (mxMKan$bx + mxFKan$bx)
	bx.lim=c(min(bx[bx>0]), max(bx[bx>0]))
	lmin <- -12
	lmax <- 0
	machine.max <- log(.Machine$double.xmax)
	machine.min <- log(.Machine$double.xmin)
	klM <- min((lmax - min(mxMKan$ax))/bx.lim[2], machine.max)
	klF <- min((lmax - min(mxFKan$ax))/bx.lim[2], machine.max)
	kuM <- max((lmin - max(mxMKan$ax))/bx.lim[1], machine.min)
	kuF <- max((lmin - max(mxFKan$ax))/bx.lim[1], machine.min)
	return(list(male=mxMKan, female=mxFKan, 
				bx=bx, kl=c(klM, klF), ku=c(kuM, kuF)))
}

modifiedLC <- function (npred, mxKan, eopm, eopf, verbose=FALSE, debug=FALSE) {
	eop  <- list(eopm, eopf)
    # Using combined bx, This differs from ModifiedLC0!    
    sr <- LLm <- list(matrix(0, nrow=27, ncol=npred), matrix(0, nrow=27, ncol=npred))
    Mx <- lx <- list(matrix(0, nrow=28, ncol=npred), matrix(0, nrow=28, ncol=npred))
    #sr1 <- list(matrix(0, nrow=27, ncol=npred), matrix(0, nrow=27, ncol=npred))
    #if(debug) print('Start check ===========================')
    #Get the projected kt from eo, and make projection of Mx
    #stop('')
    for (mxYKan in list(mxKan$female, mxKan$male)) { # iterate over male and female
    	#print(c('sex: ', mxYKan$sex))
    	res <- .C("LC", as.integer(npred), as.integer(mxYKan$sex), as.numeric(mxYKan$ax), as.numeric(mxKan$bx), 
			as.numeric(eop[[mxYKan$sex]]), Kl=as.numeric(mxKan$kl[[mxYKan$sex]]), Ku=as.numeric(mxKan$ku[[mxYKan$sex]]), 
			constrain=as.integer(mxYKan$sex == 1), 
			#constrain=as.integer(0),
			FMx=as.numeric(Mx[[2]]), FEop=as.numeric(eop[[2]]),
			LLm = as.numeric(LLm[[mxYKan$sex]]), Sr=as.numeric(sr[[mxYKan$sex]]), 
			lx=as.numeric(lx[[mxYKan$sex]]), Mx=as.numeric(Mx[[mxYKan$sex]]))
		sr[[mxYKan$sex]] <- matrix(res$Sr, nrow=27)
		LLm[[mxYKan$sex]] <- matrix(res$LLm, nrow=27)
		Mx[[mxYKan$sex]] <- matrix(res$Mx, nrow=28)
		lx[[mxYKan$sex]] <- matrix(res$lx, nrow=28)
    }
    #stop('')
	return(list(sr=sr, LLm=LLm, mx=Mx, lx=lx))    
}

KannistoAxBx.joint <- function(ne, male.mx, female.mx, yb)  {
	# Extending mx to age 130 using Kannisto model and mx 80-99, OLS
	# ne is the number of projections
	Mxe.m <- as.matrix(male.mx)
	Mxe.m <- rbind(Mxe.m, 
			matrix(NA, nrow=28-nrow(male.mx), ncol=ncol(male.mx)))			
	Mxe.f <- as.matrix(female.mx)
	Mxe.f <- rbind(Mxe.f, 
			matrix(NA, nrow=28-nrow(female.mx), ncol=ncol(female.mx)))

	k <- 22:28
	npoints <- 4
	#h <- 100 + 5 * (0:6)
	#hminus80 <- h - 80
	age.group <- (21-npoints+1):21
	data <- data.frame(sex=c(rep(1,npoints), rep(0,npoints)), age=c(age.group, age.group))
	mxc <- rbind(male.mx[age.group, 1:ne], female.mx[age.group, 1:ne])
	logit.mxc <- log(mxc) - log(1-mxc)
	#lmxr <- log(mx[18:21, 1:ne] / (1 - mx[18:21, 1:ne]))
	#Xm1 <- apply(lmxr, 2, sum)
	#Xm2 <- apply((5 * (1:4) - 5) * lmxr, 2, sum)
	for(j in 1:ne) {
		data$lmx <- logit.mxc[,j]
		fit <- lm(lmx ~ sex + age, data=data)
		coefs <- coefficients(fit)
		aam.female <- exp(coefs[1]) # intercept
		aam.male <- exp(coefs[1] + coefs['sex'])
		bbm <- coefs['age']
		# Ages 100-105, ..., 130+
		expterm.m <- aam.male * exp(bbm * k)	
		Mxe.m[k, j] =  expterm.m / (1 + expterm.m)
		expterm.f <- aam.female * exp(bbm * k)	
		Mxe.f[k, j] =  expterm.f / (1 + expterm.f)
	}
    
	#Get Lee-Cater Ax and Bx
	if(yb > 1980) yb <- 1980
	ns <- (yb - 1950) / 5 + 1
    
    result <- list(male=list(mx=Mxe.m), female=list(mx=Mxe.f))
    for(sex in c('male', 'female')) {
    	lMxe <- log(result[[sex]]$mx)
    	x1 <- apply(lMxe[,ns:ne], 1, sum)
    	ax <- x1 / (ne - ns + 1)
		kt <- rep(NA, ne)
		kt[ns:ne] = apply(lMxe[,ns:ne], 2, sum) - sum(ax)
    
		x2 <- sum(kt[ns:ne]*kt[ns:ne])
		x1 <- rep(NA, nrow(lMxe))
		for (i in 1:nrow(lMxe)) 
			x1[i] <- sum((lMxe[i,ns:ne]-ax[i])*kt[ns:ne])
		bx <- x1/x2
		negbx <- which(bx <= 0)
		lnegbx <- length(negbx)
		if(lnegbx > 0 && negbx[1] == 1) {
			bx[1] <- 0
			negbx <- if(lnegbx > 1) negbx[2:lnegbx] else c()
			lnegbx <- length(negbx)
		}
		while(lnegbx > 0) {
			bx[negbx] <- 0.5 * bx[negbx-1]
			negbx <- which(bx < 0)
			lnegbx <- length(negbx)
		}      
		for (i in 1:27) { # ?
			if (bx[29 - i] == 0) bx[29 - i] <- bx[29 - i - 1]
		}
    	bx[19:28] <- bx[18]   # bx(80) is used for all older ages to fit large Eo
		result[[sex]]$ax <- ax
		result[[sex]]$bx <- bx
		result[[sex]]$k0 <- kt[ne]
		result[[sex]]$d1 <- (kt[ne] - kt[ns]) / (ne - ns + 1)
	}
	return(result)
}


KannistoAxBx <- function(ne, mx, yb, female.mx=NULL, male.mx=NULL)  {
	# Extending mx to age 130 using Kannisto model and mx 80-99, OLS
	# ne is the number of projections
	Mxe <- as.matrix(mx)
	Mxe <- rbind(Mxe, 
			matrix(NA, nrow=28-nrow(mx), ncol=ncol(mx)))

	k <- 22:28
	h <- 100 + 5 * (0:6)
	hminus80 <- h - 80
	lmxr <- log(mx[18:21, 1:ne] / (1 - mx[18:21, 1:ne]))
	Xm1 <- apply(lmxr, 2, sum)
	Xm2 <- apply((5 * (1:4) - 5) * lmxr, 2, sum)
	for(j in 1:ne) {
		aam = exp((350 * Xm1[j] - 30 * Xm2[j]) / 500) # 500 = 4 * 350 - 30^2
		bbm = (Xm1[j] - 4 * log(aam)) / 30
		# Ages 100-105, ..., 130+
		expterm <- aam * exp(bbm * (hminus80))	
		Mxe[k, j] =  expterm / (1 + expterm)
		#if(!is.null(female.mx) && any(Mxe[k, j] < female.mx[k,j])) {
		#	smaller.f <- which(Mxe[k, j] < female.mx[k,j])
		#	Mxe[k[smaller.f],j] <- female.mx[k[smaller.f],j]
		#}
	}
    lMxe <- log(Mxe)
    
	#Get Lee-Cater Ax and Bx
	if(yb > 1980) yb <- 1980
	ns <- (yb - 1950) / 5 + 1
    
    x1 <- apply(lMxe[,ns:ne], 1, sum)
    ax <- x1 / (ne - ns + 1)
    #if(!is.null(male.mx))
	#	ax[(male.mx$ax - ax)<0.2] <- male.mx$ax[(male.mx$ax - ax)<0.2]-0.2
	kt <- rep(NA, ne)
	kt[ns:ne] = apply(lMxe[,ns:ne], 2, sum) - sum(ax)
    
	x2 <- sum(kt[ns:ne]*kt[ns:ne])
	x1 <- rep(NA, nrow(Mxe))
	for (i in 1:nrow(Mxe)) 
		x1[i] <- sum((lMxe[i,ns:ne]-ax[i])*kt[ns:ne])
	bx <- x1/x2
	negbx <- which(bx <= 0)
	lnegbx <- length(negbx)
	if(lnegbx > 0 && negbx[1] == 1) {
		#warning('First value of bx is negative.', immediate. = TRUE)
		#stop('')
		#posbx <- which(bx > 0)[1]
		bx[1] <- 0
		negbx <- if(lnegbx > 1) negbx[2:lnegbx] else c()
		lnegbx <- length(negbx)
	}
	while(lnegbx > 0) {
		bx[negbx] <- 0.5 * bx[negbx-1]
		negbx <- which(bx < 0)
		lnegbx <- length(negbx)
	}
       
	for (i in 1:27) { # ?
		if (bx[29 - i] == 0) bx[29 - i] <- bx[29 - i - 1]
	}
    bx[19:28] <- bx[18]   # bx(80) is used for all older ages to fit large Eo

	return(list(mx=Mxe, ax=ax, bx=bx, k0=kt[ne], d1=(kt[ne] - kt[ns]) / (ne - ns + 1)))
}


StoPopProj <- function(npred, inputs, LT, asfr, mig.type, country.name, keep.vital.events=FALSE) {
	popm <- popf <- matrix(0, nrow=27, ncol=npred+1)
	popm[,1] <- c(inputs$POPm0, rep(0, 6))
	popf[,1] <- c(inputs$POPf0, rep(0, 6))
	totp <- c(sum(popm[,1]+popf[,1]), rep(0, npred))
	btageM <- btageF <- matrix(0, nrow=7, ncol=npred) # births by age of mother and sex of child
	deathsM <- deathsF <- matrix(0, nrow=27, ncol=npred)
	#stop('')
	res <- .C("TotalPopProj", npred, as.numeric(as.matrix(inputs$MIGm)), 
		as.numeric(as.matrix(inputs$MIGf)), nrow(inputs$MIGm), ncol(inputs$MIGm), mig.type,
		srm=LT$sr[[1]], srf=LT$sr[[2]], asfr=as.numeric(as.matrix(asfr)), 
		srb=as.numeric(as.matrix(inputs$SRB)), popm=popm, popf=popf, totp=totp,
		Lm=as.numeric(LT$LLm[[1]]), Lf=as.numeric(LT$LLm[[2]]), 
		lxm=as.numeric(LT$lx[[1]]), lxf=as.numeric(LT$lx[[2]]), 
		btagem=as.numeric(btageM), btagef=as.numeric(btageF), 
		deathsm=as.numeric(deathsM), deathsf=as.numeric(deathsF)
		)
	
 	if(any(res$popm < 0)) warnings('Negative male population for ', country.name)
	if(any(res$popf < 0)) warnings('Negative female population for ', country.name)

	vital.events <- list()
	if(keep.vital.events) {
		vital.events$mbt <- res$btagem
		vital.events$fbt <- res$btagef
		vital.events$mdeaths <- res$deathsm
		vital.events$fdeaths <- res$deathsf
	}
	return(c(list(totpop=res$totp, mpop=res$popm, fpop=res$popf), vital.events))
}

compute.observedVE <- function(inputs, pop.matrix, mig.type, mxKan, country.code, estim.years) {
	obs <- inputs$observed
	if(is.null(obs$PASFR)) return(NULL)
	npasfr <- nrow(obs$PASFR)
	nest <- min(length(obs$TFRpred), ncol(obs$PASFR))
	estim.years <- estim.years[(length(estim.years)-nest+1):length(estim.years)]
	asfr <- obs$PASFR[,(ncol(obs$PASFR)-nest+1):ncol(obs$PASFR), drop=FALSE]/100.
	tfr <- obs$TFRpred[(length(obs$TFRpred)-nest+1):length(obs$TFRpred)]
	mig.data <- list(as.matrix(obs$MIGm[,(ncol(obs$MIGm)-nest+1):ncol(obs$MIGm)]), 
					as.matrix(obs$MIGf[,(ncol(obs$MIGf)-nest+1):ncol(obs$MIGf)]))
	for(i in 1:npasfr) asfr[i,] <- tfr * asfr[i,]
	
	pop <- D10 <- list()
	nmx <- ncol(inputs$MXm)
	mx <-  list(inputs$MXm[,(nmx-nest+1):nmx, drop=FALSE], inputs$MXf[,(nmx-nest+1):nmx, drop=FALSE])
	srb <- obs$SRB[(length(obs$SRB)-nest+1):length(obs$SRB)]
	srb.ratio <- srb / (1 + srb)
	sr <- list(matrix(0, nrow=21, ncol=nest), matrix(0, nrow=21, ncol=nest))
	births <- deaths <- list()
	for(sex in 1:2) {
		tpop <- get.pop.observed.with.age(NULL, country.code, sex=c('male', 'female')[sex], 
						data=pop.matrix)
		pop[[sex]] <- tpop$data[tpop$age.idx,(ncol(tpop$data)-nest):ncol(tpop$data)]
	}
	bt <- (pop[[2]][4:10,2:ncol(pop[[2]])] + pop[[2]][4:10,1:(ncol(pop[[2]])-1)]) * asfr * 0.5
	births[[1]] <- bt * srb.ratio
	births[[2]] <- bt - births[[1]]
	for(sex in 1:2) {		
		sr[[sex]] <- get.survival(abind(mx[[sex]],along=3), sex=c("M","F")[sex], age05=c(TRUE, TRUE, FALSE))[,,1]
		deaths[[sex]] <- matrix(0, nrow=21, ncol=nest)
		res <- .C("get_deaths_from_sr", as.numeric(sr[[sex]]), nest, as.numeric(as.matrix(pop[[sex]])), 
					as.numeric(mig.data[[sex]]), mig.type,
					as.numeric(colSums(as.matrix(births[[sex]]))), Deaths=as.numeric(deaths[[sex]]))
		deaths[[sex]] <- matrix(res$Deaths, nrow=21)
		colnames(deaths[[sex]]) <- estim.years
		rownames(deaths[[sex]]) <- rownames(pop[[sex]])
	}	
	colnames(asfr) <- estim.years
	rownames(asfr) <- rownames(births[[1]])
	res <- list(btm=births[[1]], btf=births[[2]], 
				deathsm=deaths[[1]], deathsf=deaths[[2]], asfert=asfr, mxm=mx[[1]], mxf=mx[[2]])
	for(par in names(res))
		res[[par]] <- abind(res[[par]], NULL, along=3) # add trajectory dimension
	return(res)
}

write.pop.projection.summary <- function(pop.pred, what=NULL, expression=NULL, output.dir=NULL, ...) {
	if (is.null(output.dir)) output.dir <- pop.pred$output.directory
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.what <- c('pop', 'popsex', 'popsexage', 'popage', 'births', 'birthssex', 'birthsage', 'birthssexage', 
			'deaths', 'deathssex', 'deathsage', 'deathssexage', 'srsexage', 'fertility', 'fertilityage')
	what <- if(is.null(what)) all.what else match.arg(what, all.what, several.ok=TRUE)
	params <- list()
	if(!is.null(expression)) {
		what <- 'expression'
		params <- list(expression=expression)
	}
	for(summary.type in what) {
		if(is.element(summary.type, what))
			do.call(paste('write.', summary.type, sep=''), c(list(pop.pred, output.dir=output.dir), params, ...))
	}
}

write.pop <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE)
	
write.popsex <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, what.log='population')
	
write.popsexage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, what.log='population')
	
write.popage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, what.log='population')
	
write.births <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='births', 
			file.suffix='births', what.log='total births')
	
write.birthssex <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, vital.event='births', 
			file.suffix='births', what.log='births')

write.birthsage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='births', 
			file.suffix='births', what.log='births')

write.birthssexage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='births', 
			file.suffix='births', what.log='births')

write.deaths <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='deaths', 
			file.suffix='deaths', what.log='total deaths')
	
write.deathssex <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths')

write.deathsage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths')
	
write.deathssexage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths')
			
write.srsexage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='survival', 
			file.suffix='sr', what.log='survival ratio')
			
write.fertility <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='fertility', 
			file.suffix='asfr', what.log='fertility rate')
			
write.fertilityage <- function(pop.pred, output.dir) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='fertility', 
			file.suffix='asfr', what.log='fertility rate')	
	
write.expression <- function(pop.pred, expression, output.dir, file.suffix='expression', expression.label=expression,  
								include.observed=FALSE, decimal=NULL) {
	cat('Creating summary file for expression ', expression, ' ...\n')
	header <- list(country.name='country_name',  country.code='country_code', variant='variant')
	variant.names <- c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95')
	nr.var <- length(variant.names)
	pred.period <- get.pop.prediction.periods(pop.pred)
	nr.proj <- length(pred.period)
	nr.obs <- 0
	if(include.observed) {
		obs.period <- get.pop.observed.periods(pop.pred)
		nr.obs <- length(obs.period)-1
		obs.period <- obs.period[-(nr.obs+1)] # remove the last one because the same as the first projection period
		for(iyear in 1:nr.obs) header[[paste('year', iyear, sep='')]] <- obs.period[iyear]
	}
	for(iyear in 1:nr.proj) header[[paste('year', iyear+nr.obs, sep='')]] <- pred.period[iyear]
	col.names <- grep('year', names(header), value=TRUE)
	result <- NULL
	if(include.observed) {
		res <- get.pop.observed.from.expression.all.countries(expression, pop.pred, time.index=1:nr.obs)
		#copy the same data into the variant rows 
		result <- matrix(NA, nrow=nrow(res)*5, ncol=ncol(res))
		for(i in 1:5) result[seq(i,by=5, length=nrow(res)),] <- res
	}	
	for(iyear in 1:nr.proj) {	
		result <- cbind(result, as.vector(t(get.pop.from.expression.all.countries(expression, pop.pred, 
						quantiles=c(0.5, 0.1, 0.9, 0.025, 0.975), projection.index=iyear))))
	}
	if(!is.null(decimal)) result <- round(result, decimal)
	colnames(result) <- col.names
	result <- cbind(country.name=rep(as.character(pop.pred$countries$name), each=nr.var), 
					country.code=rep(pop.pred$countries$code, each=nr.var),
					variant=rep(variant.names, length(pop.pred$countries$code)), result)
	colnames(result)[colnames(result)==names(header)] <- header
	file <- file.path(output.dir, paste('projection_summary_', file.suffix, '.csv', sep=''))
	write(paste('#', expression.label), file)
	warn <- getOption('warn')
	options(warn=-1) # disable warning messages (it doesn't like that col.names is TRUE and append=TRUE)
	write.table(result, file=file, sep=',', row.names=FALSE, col.names=TRUE, append=TRUE,
				quote=which(is.element(colnames(result), c('country_name', 'variant'))))
	options(warn=warn)
	cat('Stored into: ', file, '\n')
}	
	
.write.pop <- function(pop.pred, output.dir, bysex=FALSE, byage=FALSE, vital.event=NULL, file.suffix='tpop', 
							what.log='total population') {
	cat('Creating summary file of ', what.log, ' ')
	if(bysex) cat('by sex ')
	if(byage) cat('by age ')
	cat ('...\n')
	header <- list(country.name='country_name',  country.code='country_code', variant='variant')
	if(bysex) header[['sex']] <- 'sex'
	if(byage) header[['age']] <- 'age'
	variant.names <- c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95')
	nr.var <- length(variant.names)
	pred.period <- get.pop.prediction.periods(pop.pred)
	if(!is.null(vital.event)) pred.period <- pred.period[2:length(pred.period)]
	nr.proj <- length(pred.period)
	for (i in 1:nr.proj) 
		header[[paste('year', i, sep='')]] <- pred.period[i]
	col.names <- grep('year', names(header), value=TRUE)
	result <- NULL
	sex.index <- c(TRUE, FALSE, FALSE)
	if(bysex) sex.index <- !sex.index
	age.index <- c(TRUE, rep(FALSE, length(pop.pred$ages)))
	if(byage) age.index <- !age.index
	ages <- 1:length(pop.pred$ages)
	for (country in 1:nrow(pop.pred$countries)) {
		country.obj <- get.country.object(country, country.table=pop.pred$countries, index=TRUE)
		for(sex in c('both', 'male', 'female')[sex.index]) {
			if(!is.null(vital.event)) {
			 	sum.over.ages <- age.index[1]
				traj.and.quantiles <- get.popVE.trajectories.and.quantiles(pop.pred, country.obj$code, event=vital.event, 
										sex=sex, age='all', sum.over.ages=sum.over.ages)
				if(is.null(traj.and.quantiles$trajectories)) {
					warning('Problem with loading ', vital.event, '. Possibly no vital events stored during prediction.')
					return()
				}
				if(!sum.over.ages) { # This is because births have only subset of ages
					ages <- traj.and.quantiles$age.idx.raw
					age.index <- age.index[1:(length(ages)+1)]
					subtract.from.age <- traj.and.quantiles$age.idx.raw[1]-traj.and.quantiles$age.idx[1]
				}
				
			}
			for(age in c('all', ages)[age.index]) {
				this.result <- cbind(
							country.name=rep(country.obj$name, nr.var), 
							country.code=rep(country.obj$code, nr.var),
							variant=variant.names)
				if(sex != 'both')
					this.result <- cbind(this.result, sex=rep(sex, nr.var))
				if(age != 'all') {
					age <- as.integer(age)
					this.result <- cbind(this.result, age=rep(get.age.labels(pop.pred$ages)[age], nr.var))
				}
				if(is.null(vital.event)) {
					quant <- get.pop.trajectories(pop.pred, country.obj$code, nr.traj=0, sex=sex, age=age)$quantiles
					traj <- NULL
					reload <- TRUE
				} else { # vital event
					quant <- traj.and.quantiles$quantiles
					traj <- traj.and.quantiles$trajectories
					if(!sum.over.ages) {
						quant <- quant[age-subtract.from.age,,]
						traj <- traj[age-subtract.from.age,,]
					}
					reload <- FALSE
				}
				proj.result <- round(rbind(
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, q=0.5, 
											trajectories=traj, reload=reload), 
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, pi=80, 
											trajectories=traj, reload=reload),
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, pi=95, 
											trajectories=traj, reload=reload)),
					0)
				colnames(proj.result) <- col.names
				this.result <- cbind(this.result, proj.result)
				result <- rbind(result, this.result)
			}
		}
	}
	colnames(result)[colnames(result)==names(header)] <- header
	suffix <- paste(file.suffix, paste(c('sex', 'age')[c(bysex, byage)], collapse=''), sep='')
	file <- paste('projection_summary_', suffix, '.csv', sep='')
	write.table(result, file=file.path(output.dir, file), sep=',', row.names=FALSE, col.names=TRUE, 
				quote=which(is.element(colnames(result), c('country_name', 'variant', 'sex', 'age'))))
	cat('Stored into: ', file.path(output.dir, file), '\n')
}

LifeTableMxCol <- function(mx, colname=c('Lx', 'lx', 'qx', 'mx'), ...){
	colname <- match.arg(colname)
	if(is.null(dim(mx))) return(.doLifeTableMxCol(mx, colname, ...))
	return(apply(mx, 2, .doLifeTableMxCol, colname=colname, ...))
}

.collapse.Lx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	Lx.start <- c(LT$Lx[1:2], LT$Lx[1] + LT$Lx[2])[age05]
	return(c(Lx.start, LT$Lx[-(1:2)]))
}

.collapse.lx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	lx.start <- LT$lx[c(1,2,2)][age05]
	return(c(lx.start, LT$lx[-(1:2)]))
}

.collapse.qx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	qx.start <- c(LT$qx[1:2], 1-(LT$lx[3]/LT$lx[1]))[age05]
	return(c(qx.start, LT$qx[-(1:2)]))
}

.collapse.mx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	mx.start <- c(LT$mx[1:2], (LT$lx[1] - LT$lx[3])/(LT$Lx[1] + LT$Lx[2]))[age05]
	return(c(mx.start, LT$mx[-(1:2)]))
}

.doLifeTableMxCol <- function(mx, colname, age05=c(FALSE, FALSE, TRUE), ...) {
	# age05 determines the inclusion of ages 0-1, 1-4, 0-4
	LT <- LifeTableMx(mx, ...)
	return(do.call(paste('.collapse', colname, sep='.'), list(LT, age05=age05)))
}


LifeTableMx <- function(mx, sex=c('Male', 'Female')){
	
	sex <- match.arg(sex)
	sex <- list(Male=1, Female=2)[[sex]]
	nage <- length(mx)
	Lx <- lx <- qx <- rep(0, nage)
	ax <- rep(0, 21)
	nagem1 <- nage-1
	LTC <- .C("LifeTable", as.integer(sex), as.integer(nagem1), as.numeric(mx), 
			Lx=Lx, lx=lx, qx=qx, ax=ax)
	LT <- data.frame(age=c(0,1, seq(5, by=5, length=nage-2)), 
					Lx=LTC$Lx, lx=LTC$lx, qx=LTC$qx, ax=rep(NA,nage), mx=mx)
	l <- min(length(LTC$ax), nrow(LT))
	LT$ax[1:l] <- LTC$ax[1:l]
	return(LT)
}

unblock.gtk <- function(...) bayesTFR:::unblock.gtk(...)