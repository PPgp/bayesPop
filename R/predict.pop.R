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
						), nr.traj = 1000, 
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
	} 
	
	data(LOCATIONS)
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
			country.codes <- countries
		} else
			country.codes <- unique(inp$POPm0[,'country_code'])
	}
	invisible(do.pop.predict(country.codes, inp, outdir, nr.traj, ages, pred=if(prediction.exist) pred else NULL,
					verbose=verbose))
}

do.pop.predict <- function(country.codes, inp, outdir, nr.traj, ages, pred=NULL, verbose=FALSE) {
	ncountries <- length(country.codes)
	nr_project <- length(inp$proj.years)
	nages <- length(ages)
	nest <- length(inp$estim.years)
	if(!file.exists(outdir)) 
		dir.create(outdir, recursive=TRUE)
	prediction.file <- file.path(outdir, 'prediction.rda')	
	quantiles.to.keep <- c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1)
	PIs_cqp <- quantM <- quantF <- array(NA, c(ncountries, length(quantiles.to.keep), nr_project+1),
						dimnames=list(NULL, quantiles.to.keep, c(inp$estim.years[nest], inp$proj.years)))
	quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(ncountries, nages, length(quantiles.to.keep), nr_project+1),
						dimnames=list(NULL, ages, quantiles.to.keep, c(inp$estim.years[nest], inp$proj.years)))
	mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(ncountries, 2, nr_project+1))
	country.idx.proj <- c()
	cidx <- 0
	for(country in country.codes) {
		country.idx <- which(LOCATIONS[,'country_code'] == country)
		if(length(country.idx) == 0) {
			warning('Country ', country, ' not found in the LOCATIONS dataset.')
			next
		}
		if(verbose)
			cat('\nProcessing country ', country, ' -- ', as.character(LOCATIONS[country.idx,'name']))
		# Extract the country-specific stuff from the inputs
		inpc <- get.country.inputs(country, inp, nr.traj, LOCATIONS[country.idx,'name'])
		if(is.null(inpc)) next
		nr.traj <- min(ncol(inpc$TFRpred), nr.traj)		
		if(verbose)
			cat(' (', nr.traj, ' trajectories )')
		
		npred <- min(nrow(inpc$TFRpred), nr_project)
		#npred <- 19
		totp <- matrix(NA, nrow=npred+1, ncol=nr.traj)
		totpm <- totpf <- array(NA, dim=c(27, npred+1, nr.traj))
		debug <- FALSE
		MxKan <- runKannisto(nest, inpc)
		npasfr <- nrow(inpc$PASFR)			
		for(itraj in 1:nr.traj) {
			asfr <- inpc$PASFR/100.
			for(i in 1:npasfr) asfr[i,] <- inpc$TFRpred[,itraj] * asfr[i,]
			#debug <- country == 478
			#if(itraj == 2) debug <- TRUE
			sr <- modifiedLC(nest, npred, MxKan, inpc$e0Mpred[,itraj], 
									inpc$e0Fpred[,itraj], verbose=verbose, debug=debug)
			popres <- StoPopProj(npred, inpc, sr, asfr, inpc$MIGtype, LOCATIONS[country.idx,'name'])
			totp[,itraj] <- popres$totpop
			totpm[,,itraj] <- popres$mpop
			totpf[,,itraj] <- popres$fpop
		}
		totp.hch <- matrix(NA, nrow=npred+1, ncol=2)
		totpm.hch <- totpf.hch <- array(NA, dim=c(27, npred+1, 2))
		for (variant in 1:nrow(inpc$TFRhalfchild)) { # compute the two half child variants
			asfr <- inpc$PASFR/100.
			for(i in 1:npasfr) asfr[i,] <- inpc$TFRhalfchild[variant,] * asfr[i,]
			sr <- modifiedLC(nest, npred, MxKan, inpc$e0Mmedian, 
									inpc$e0Fmedian, verbose=verbose, debug=debug)
			popres <- StoPopProj(npred, inpc, sr, asfr, inpc$MIGtype, LOCATIONS[country.idx,'name'])
			totp.hch[,variant] <- popres$totpop
			totpm.hch[,,variant] <- popres$mpop
			totpf.hch[,,variant] <- popres$fpop
		}
		country.idx.proj <- c(country.idx.proj, country.idx)
		cidx <- cidx + 1
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch,
			 file = file.path(outdir, paste('totpop_country', country, '.rda', sep='')))
		PIs_cqp[cidx,,] = apply(totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sd[cidx,1,] <- apply(totp, 1, mean, na.rm = TRUE)
		mean_sd[cidx,2,] = apply(totp, 1, sd, na.rm = TRUE)
		for (i in 1:nages) {
			quantMage[cidx,i,,] <- apply(totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
			quantFage[cidx,i,,] = apply(totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
			quantPropMage[cidx,i,,] <- apply(totpm[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			quantPropFage[cidx,i,,] <- apply(totpf[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		}
		stotpm <- colSums(totpm)
		quantM[cidx,,] = apply(stotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdM[cidx,1,] <- apply(stotpm, 1, mean, na.rm = TRUE)
		mean_sdM[cidx,2,] = apply(stotpm, 1, sd, na.rm = TRUE)
		stotpf <- colSums(totpf)
		quantF[cidx,,] = apply(stotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdF[cidx,1,] <- apply(stotpf, 1, mean, na.rm = TRUE)
		mean_sdF[cidx,2,] = apply(stotpf, 1, sd, na.rm = TRUE)
		
		#save updated meta file
		country.row <- LOCATIONS[country.idx,c('country_code', 'name')]
		colnames(country.row) <- c('code', 'name')
		if(cidx == 1) { # first pass
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
	return(read.pop.file(file.path(.find.package("bayesPop"), "data", file)))

load.inputs <- function(inputs, start.year, present.year, end.year, wpp.year) {
	# Get initial population counts
	if(is.null(inputs$popM)) 
		POPm0 <- read.bayesPop.file(paste('PopByAgeMale', wpp.year, '.txt', sep=''))	else POPm0 <- read.pop.file(inputs$popM)
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
	SRB <- SRB[,c('country_code', proj.periods)]
	mid.proj.years <- cols.starty[start.index:end.index] + 3
	
	# Get percentage age-specific fertility rate
	if(is.null(inputs$pasfr)) 
		PASFR <- read.bayesPop.file(paste('PercentASFR', wpp.year, '.txt', sep=''))
	else PASFR <- read.pop.file(inputs$pasfr)
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
	MIGm <- MIGm[,c('country_code', 'age', proj.periods)]
	if(is.null(inputs$migF))
		MIGf <- read.bayesPop.file(paste('MigByAgeFemale', wpp.year, '.txt', sep=''))
	else MIGf <- read.pop.file(inputs$migF)
	MIGf <- MIGf[,c('country_code', 'age', proj.periods)]
	
	# Get life expectancy
	if(!is.null(inputs$e0F.sim.dir))  # female
		e0Fpred <- get.e0.prediction(inputs$e0F.sim.dir, mcmc.dir=NA)
	else {
		file.name <- if(!is.null(inputs$e0F.file)) inputs$e0F.file
					else file.path(.find.package("bayesPop"), "ex-data", 
							paste('e0Fwpp', wpp.year, '.csv', sep=''))
		if(!file.exists(file.name))
			stop('File ', file.name, 
				' does not exist.\nSet e0F.sim.dir, e0F.file or change WPP year.')
		e0Fpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
		e0Fpred <- e0Fpred[,c('LocID', 'Year', 'e0')]
		colnames(e0Fpred) <- c('country_code', 'year', 'value')
	} 

	if(!is.null(inputs$e0M.sim.dir)) { # male
		if(inputs$e0M.sim.dir == 'joint_') {
			if(!has.e0.jmale.prediction(e0Fpred))
				stop('No joint prediction for female and male available. Correct the e0M.sim.dir argument.' )
			e0Mpred <- get.e0.jmale.prediction(e0Fpred)
		} else e0Mpred <- get.e0.prediction(inputs$e0M.sim.dir, mcmc.dir=NA)
	} else {
		file.name <- if(!is.null(inputs$e0M.file)) inputs$e0M.file 
					else file.path(.find.package("bayesPop"), "ex-data", 
							paste('e0Mwpp', wpp.year, '.csv', sep=''))
		if(!file.exists(file.name))
			stop('File ', file.name, 
				' does not exist.\nSet e0M.sim.dir, e0M.file or change WPP year.')
		e0Mpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
		e0Mpred <- e0Mpred[,c('LocID', 'Year', 'e0')]
		colnames(e0Mpred) <- c('country_code', 'year', 'value')
	} 
		
	# Get TFR
	if(!is.null(inputs$tfr.sim.dir)) 
		TFRpred <- get.tfr.prediction(inputs$tfr.sim.dir, mcmc.dir=NA)
	else {
		file.name <- if(!is.null(inputs$tfr.file)) inputs$tfr.file
					else file.path(.find.package("bayesPop"), "ex-data", 
							paste('TFRwpp', wpp.year, '.csv', sep=''))
		if(!file.exists(file.name))
			stop('File ', file.name, 
				' does not exist.\nSet tfr.sim.dir, tfr.file or change WPP year.')
		TFRpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
		TFRpred <- TFRpred[,c('LocID', 'Year', 'TF')]
		colnames(TFRpred) <- c('country_code', 'year', 'value')
	} 
	return(list(POPm0=POPm0, POPf0=POPf0, MXm=MXm, MXf=MXf, SRB=SRB,
				PASFR=PASFR, MIGtype=MIGtype, MIGm=MIGm, MIGf=MIGf,
				e0Mpred=e0Mpred, e0Fpred=e0Fpred, TFRpred=TFRpred, 
				estim.years=mid.est.years, proj.years=mid.proj.years,
				pop.matrix=list(male=popm.matrix, female=popf.matrix)))
}

get.country.inputs <- function(country, inputs, nr.traj, country.name) {
	inpc <- list()
	for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'SRB',
						'PASFR', 'MIGtype', 'MIGm', 'MIGf')) {
		idx <- inputs[[par]][,'country_code'] == country
		inpc[[par]] <- inputs[[par]][idx,,drop=FALSE]
		inpc[[par]] <- as.matrix(inpc[[par]][, !is.element(colnames(inpc[[par]]), 
							c('country_code', 'age')),drop=FALSE])
	}
	inpc[['MIGBaseYear']] <- inpc[['MIGtype']][,'ProjFirstYear']
	inpc[['MIGtype']] <- inpc[['MIGtype']][,'MigCode']
	what.traj <- list(TFRpred='TFR', e0Mpred='male e0', e0Fpred='female e0')
	medians <- list()
	for(par in names(what.traj)) {
		if (is.data.frame(inputs[[par]])) {
			idx <- inputs[[par]][,'country_code'] == country & is.element(inputs[[par]][,'year'], inputs$proj.years)
			if(sum(idx) == 0) {
				warning('No ', what.traj[[par]], ' trajectories for ', country.name, 
					'. No population projection generated.')
				return(NULL)
			}
			inpc[[par]] <- inputs[[par]][idx,]
			inpc[[par]] <- matrix(inpc[[par]][, 'value'], nrow=length(inputs$proj.years))
			rownames(inpc[[par]]) <- as.character(inputs$proj.years)
			medians[[par]] <- apply(inpc[[par]], 1, quantile, 0.5, na.rm = TRUE)
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
	} 
	inpc$TFRhalfchild <- bayesTFR:::get.half.child.variant(median=medians$TFRpred, increment=c(0.25, 0.4, 0.5))
	
	if(is.null(inpc$e0Mpred)) {
		inpc$e0Mpred <- get.e0.trajectories(inputs$e0Mpred, country)
		if(is.null(inpc$e0Mpred)) {
			warning('No male e0 trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)	
		}
		country.obj <- get.country.object(country, inputs$e0Mpred$mcmc.set$meta)
		medians$e0Mpred <- bayesTFR:::get.median.from.prediction(inputs$e0Mpred, country.obj$index, country.obj$code)
	}
	inpc$e0Mmedian <- medians$e0Mpred
	
	if(is.null(inpc$e0Fpred)) {
		inpc$e0Fpred <- get.e0.trajectories(inputs$e0Fpred, country)
		if(is.null(inpc$e0Fpred)) {
			warning('No female e0 trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)
		}
		country.obj <- get.country.object(country, inputs$e0Fpred$mcmc.set$meta)
		medians$e0Fpred <- bayesTFR:::get.median.from.prediction(inputs$e0Fpred, country.obj$index, country.obj$code)
	}
	inpc$e0Fmedian <- medians$e0Fpred
	
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
		inpc[[par]] <- inpc[[par]][,indices[[par]]]
		inpc[[par]] <- inpc[[par]][as.character(inputs$proj.years),]
	}
	return(inpc)
}

get.traj.index <- function(nr.traj, traj) {
	ncoltraj <- ncol(traj)
	if(nr.traj >= ncoltraj) return(1:ncoltraj)
	return(seq(1, ncoltraj, length=nr.traj))
}

runKannisto <- function(nest, inputs) {
	# extend mx, get LC ax,bx,k1
	mxMKan <- c(KannistoAxBx(nest, inputs$MXm, inputs$MIGBaseYear), sex=1)
	mxFKan <- c(KannistoAxBx(nest, inputs$MXf, inputs$MIGBaseYear), sex=2)
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

modifiedLC <- function (nest, npred, mxKan, eopm, eopf, verbose=FALSE, debug=FALSE) {
	eop  <- list(eopm, eopf)
    # Using combined bx, This differs from ModifiedLC0!    
    sr <- list(matrix(0, nrow=27, ncol=npred), matrix(0, nrow=27, ncol=npred))
    sr1 <- list(matrix(0, nrow=27, ncol=npred), matrix(0, nrow=27, ncol=npred))
    #if(debug) print('Start check ===========================')
    #Get the projected kt from eo, and make projection of Mx
    for (mxYKan in list(mxKan$male, mxKan$female)) { # iterate over male and female
    	res <- .C("LC", as.integer(npred), as.integer(mxYKan$sex), as.numeric(mxYKan$ax), as.numeric(mxYKan$bx), 
			as.numeric(eop[[mxYKan$sex]]), Kl=as.numeric(mxKan$kl[[mxYKan$sex]]), Ku=as.numeric(mxKan$ku[[mxYKan$sex]]), 
			LLm = rep(0, 27), Sr=as.numeric(sr[[mxYKan$sex]]))
		sr[[mxYKan$sex]] <- matrix(res$Sr, nrow=27)
    }
	return(sr)    
}

KannistoAxBx <- function(ne, mx, yb)  {
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
	}
    lMxe <- log(Mxe)
    
	#Get Lee-Cater Ax and Bx
	if(yb > 1980) yb <- 1980
	ns <- (yb - 1950) / 5 + 1
    
    x1 <- apply(lMxe[,ns:ne], 1, sum)
    ax <- x1 / (ne - ns + 1)

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

	return(list(mx=mx, ax=ax, bx=bx, k0=kt[ne], d1=(kt[ne] - kt[ns]) / (ne - ns + 1)))
}


StoPopProj <- function(npred, inputs, sr, asfr, mig.type, country.name) {
	popm <- popf <- matrix(0, nrow=27, ncol=npred+1)	
	popm[,1] <- c(inputs$POPm0, rep(0, 6))
	popf[,1] <- c(inputs$POPf0, rep(0, 6))
	totp <- c(sum(popm[,1]+popf[,1]), rep(0, npred))
	
	res <- .C("TotalPopProj", npred, as.numeric(as.matrix(inputs$MIGm)), 
		as.numeric(as.matrix(inputs$MIGf)), nrow(inputs$MIGm), ncol(inputs$MIGm), mig.type,
		srm=sr[[1]], srf=sr[[2]], asfr=as.numeric(as.matrix(asfr)), 
		srb=as.numeric(as.matrix(inputs$SRB)), popm=popm, popf=popf, totp=totp
		)
	
 	if(any(res$popm < 0)) warnings('Negative male population for ', country.name)
	if(any(res$popf < 0)) warnings('Negative female population for ', country.name)

	return(list(totpop=res$totp, mpop=res$popm, fpop=res$popf))
}
