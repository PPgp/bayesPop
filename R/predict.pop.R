if(getRversion() >= "2.15.1") utils::globalVariables(c("UNlocations", "MLTbx"))

pop.predict <- function(end.year=2100, start.year=1950, present.year=2020, wpp.year=2019,
						countries=NULL, output.dir = file.path(getwd(), "bayesPop.output"),
						annual = FALSE,
						inputs=list(
							popM=NULL,
							popF=NULL,
							mxM=NULL,
							mxF=NULL,
							srb=NULL,
							pasfr=NULL,
							patterns=NULL,
							migM=NULL, migF=NULL,
							migMt = NULL, migFt = NULL, mig = NULL,
							e0F.file=NULL, e0M.file=NULL, 
							tfr.file=NULL, 
							e0F.sim.dir=NULL, e0M.sim.dir=NULL, 
							tfr.sim.dir=NULL,
							migMtraj=NULL, migFtraj=NULL, migtraj = NULL,
							GQpopM = NULL, GQpopF = NULL, average.annual = NULL
						), nr.traj = 1000, keep.vital.events=FALSE,
						fixed.mx=FALSE, fixed.pasfr=FALSE, lc.for.hiv = TRUE, lc.for.all = TRUE,
						my.locations.file = NULL, 
						replace.output=FALSE, verbose=TRUE, ...) {
	prediction.exist <- FALSE
	ages <- all.ages(annual, observed = FALSE)
	unblock.gtk.if.needed('reading inputs')
	if(!is.null(my.locations.file)) {
		UNlocations <- NULL # needed for R check not to complain
		UNlocations <<- read.delim(file=my.locations.file, comment.char='#', check.names=FALSE)
		if(verbose) cat('Loading ', my.locations.file, '.\n')
	} else bayesTFR:::load.bdem.dataset('UNlocations', wpp.year, envir=globalenv(), verbose=verbose)
	if(is.null(countries)) {
		if(!replace.output && has.pop.prediction(sim.dir=output.dir))
			stop('Prediction in ', output.dir,
				' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		inp <- load.inputs(inputs, start.year, present.year, end.year, wpp.year, fixed.mx=fixed.mx, fixed.pasfr=fixed.pasfr, 
		                   lc.for.hiv = lc.for.hiv, lc.for.all = lc.for.all, annual = annual, verbose=verbose)
	}else {
		if(has.pop.prediction(output.dir) && !replace.output) {
			pred <- get.pop.prediction(output.dir)
			inp <- load.inputs(pred$function.inputs, pred$inputs$start.year, pred$inputs$present.year, pred$inputs$end.year, 
								pred$wpp.year, fixed.mx=pred$inputs$fixed.mx, fixed.pasfr=pred$inputs$fixed.pasfr, all.countries=FALSE, 
								existing.mig=list(MIGm=pred$inputs$MIGm, MIGf=pred$inputs$MIGf, 
												obsMIGm=pred$inputs$observed$MIGm, obsMIGf=pred$inputs$observed$MIGf),
								lc.for.hiv = pred$inputs$lc.for.hiv, lc.for.all = pred$inputs$lc.for.all,
								annual = pred$inputs$annual, verbose=verbose)
			if(!missing(inputs)) 
				warning('Projection already exists. Using inputs from existing projection. Use replace.output=TRUE for updating inputs.')
			nr.traj <- pred$nr.traj
			ages <- pred$ages
			prediction.exist <- TRUE
		} else inp <- load.inputs(inputs, start.year, present.year, end.year, wpp.year, fixed.mx=fixed.mx, fixed.pasfr=fixed.pasfr,
		                          all.countries=FALSE, lc.for.hiv = lc.for.hiv, lc.for.all = lc.for.all, annual = annual, verbose=verbose)
	}

	outdir <- file.path(output.dir, 'predictions')
	if(!prediction.exist) {
		if(!replace.output && has.pop.prediction(sim.dir=output.dir))
			stop('Prediction in ', outdir,
			' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		unlink(outdir, recursive=TRUE)
		.remove.cache.file(outdir)
	} else pop.cleanup.cache(pred)
	if(!is.null(countries) && is.na(countries[1])) { # all countries that are not included in the existing prediction
		all.countries <- intersect(unique(inp$POPm0[,'country_code']), UNcountries())
		country.codes <- if(!prediction.exist) all.countries
						else all.countries[!is.element(all.countries, pred$countries[,'code'])]
	} else {
		if(!is.null(countries)) {
			if (is.character(countries)) { # at least one of the codes is a character string
				for (icountry in 1:length(countries)) {
					if (is.character(countries[icountry])) {
						country.idx <- which(UNlocations[,'name'] == countries[icountry])
						if(length(country.idx) > 0)
							countries[icountry] <- UNlocations[country.idx,'country_code']
					}
				}
			}
			country.codes <- as.integer(countries)
		} else
			country.codes <- intersect(unique(inp$POPm0[,'country_code']), UNcountries())
	}
	do.pop.predict(country.codes, inp, outdir, nr.traj, ages, pred=if(prediction.exist) pred else NULL,
					keep.vital.events=keep.vital.events, fixed.mx=inp$fixed.mx, fixed.pasfr=fixed.pasfr, 
					function.inputs=inputs, verbose=verbose, ...)
	invisible(get.pop.prediction(output.dir))
}

do.pop.predict <- function(country.codes, inp, outdir, nr.traj, ages, pred=NULL, keep.vital.events=FALSE, fixed.mx=FALSE, 
							fixed.pasfr=FALSE, function.inputs=NULL, verbose=FALSE, 
							parallel = FALSE, nr.nodes = NULL, ...) {
	not.valid.countries.idx <- c()
	countries.idx <- rep(NA, length(country.codes))

	for(icountry in 1:length(country.codes)) {
		country.idx <- which(UNlocations[,'country_code'] == country.codes[icountry])
		if(length(country.idx) == 0) {
			not.valid.countries.idx <- c(not.valid.countries.idx, icountry)
			next
		}
		countries.idx[icountry] <- country.idx
	}
	if(length(not.valid.countries.idx) > 0) {
		warning('Countries ', paste(country.codes[not.valid.countries.idx], collapse=', '), 
					' not found in the UNlocations dataset.')
		country.codes <- country.codes[-not.valid.countries.idx]
		countries.idx <- countries.idx[-not.valid.countries.idx]
	}
	ncountries <- length(country.codes)
	nr_project <- length(inp$proj.years)
	nages <- length(ages)
	if(!file.exists(outdir)) 
		dir.create(outdir, recursive=TRUE)
	present.and.proj.years <- c(inp$estim.years[length(inp$estim.years)], inp$proj.years)
	present.and.proj.years.pop <- present.and.proj.years
	if(!inp$annual) present.and.proj.years.pop <- present.and.proj.years.pop + 2
	prediction.file <- file.path(outdir, 'prediction.rda')	
	quantiles.to.keep <- get.quantiles.to.keep()
	nquant <- length(quantiles.to.keep)
	PIs_cqp <- quantM <- quantF <- array(NA, c(ncountries, nquant, nr_project+1),
						dimnames=list(country.codes, quantiles.to.keep, present.and.proj.years.pop))
	quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(ncountries, nages, nquant, nr_project+1),
						dimnames=list(country.codes, ages, quantiles.to.keep, present.and.proj.years.pop))
	mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(ncountries, 2, nr_project+1), 
						dimnames=list(country.codes, c('mean', 'sd'), present.and.proj.years.pop))
	mx.ages <- if(inp$annual) ages else c(0,1,ages[2:nages])
	status.for.gui <- paste('out of', ncountries, 'countries.')
	gui.options <- list()
	inp.to.save <- list()
	# remove big or redundant items from inputs to be saved
	for(item in ls(inp)[!grepl('^migMpred$|^migFpred$|^TFRpred$|^e0Fpred$|^e0Mpred$|^estim.years$|^proj.years$|^wpp.years$', ls(inp))]) 
		inp.to.save[[item]] <- get(item, inp)
		
	
	if(parallel) {
	    if(is.null(nr.nodes)) nr.nodes <- getOption("cl.cores", detectCores(logical = FALSE))
	}
	exporting.objects <- c("country.codes", "countries.idx", "UNlocations", "inp", "inp.to.save",
	                       "present.and.proj.years.pop", "present.and.proj.years", "keep.vital.events",
	                       "ages", "nages", "fixed.mx", "fixed.pasfr", "verbose", 
	                       "nquant", "quantiles.to.keep", "ncountries")

	
	# prediction function
    predict.one.country <- function(cidx, nr.traj, nr_project) {
		#unblock.gtk.if.needed(paste('finished', cidx, status.for.gui), gui.options)
		country <- country.codes[cidx]
		country.idx <- countries.idx[cidx]
		if(verbose)
			cat('\nProgress: ', round((cidx-1)/ncountries * 100), '%; now processing ', country, ' ', 
			    as.character(UNlocations[country.idx,'name']), ': ')
		# Extract the country-specific stuff from the inputs
		inpc <- get.country.inputs(country, inp, nr.traj, UNlocations[country.idx,'name'])
		if(is.null(inpc)) return(NULL)
		nr.traj <- min(ncol(inpc$TFRpred), nr.traj)		
		if(verbose)
			cat(nr.traj, ' trajectories')
		migr.modified <- .set.inp.migration.if.needed(inp, inpc, country)

		npred <- min(nrow(inpc$TFRpred), nr_project)
		npredplus1 <- npred+1
		nmortcat <- length(mx.ages)
		nfertcat <- fert.age.length(inp$annual)
		nmigcat <- all.age.length(inp$annual, observed = TRUE)
		fages <- fert.ages(inp$annual)
		totp <- matrix(NA, nrow=npredplus1, ncol=nr.traj, 
					dimnames=list(present.and.proj.years.pop, NULL))
		totpm <- totpf <- array(NA, dim=c(nages, npredplus1, nr.traj), 
							dimnames=list(ages, present.and.proj.years.pop, NULL))
		nvariants <- nrow(inpc$TFRhalfchild)
		totp.hch <- matrix(NA, nrow=npredplus1, ncol=nvariants,
					dimnames=list(present.and.proj.years.pop, NULL))
		totpm.hch <- totpf.hch <- array(NA, dim=c(nages, npredplus1, nvariants), 
					dimnames=list(ages, present.and.proj.years.pop, NULL))
		if(keep.vital.events) {
			btm <- btf <- array(0, dim=c(nfertcat, npredplus1, nr.traj), dimnames=list(fages, present.and.proj.years, NULL))
			deathsm <- deathsf <- array(0, dim=c(nages, npredplus1, nr.traj), dimnames=list(ages, present.and.proj.years, NULL))
			asfert <- array(0, dim=c(nfertcat, npredplus1, nr.traj), dimnames=list(fages, present.and.proj.years, NULL))
			pasfert <- array(0, dim=c(nfertcat, npredplus1, nr.traj), dimnames=list(fages, present.and.proj.years, NULL))
			btm.hch <- btf.hch <- array(0, dim=c(nfertcat, npredplus1, nvariants), dimnames=list(fages, present.and.proj.years, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(nages, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years, NULL))
			asfert.hch <- array(0, dim=c(nfertcat, npredplus1, nvariants), dimnames=list(fages, present.and.proj.years, NULL))
			pasfert.hch <- array(0, dim=c(nfertcat, npredplus1, nvariants), dimnames=list(fages, present.and.proj.years, NULL))
			mxm <- mxf <- array(0, dim=c(nmortcat, npredplus1, nr.traj), dimnames=list(mx.ages, present.and.proj.years, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(nmortcat, npredplus1, nvariants), dimnames=list(mx.ages, present.and.proj.years, NULL))
			migMntraj <- if(is.null(inpc[['migMpred']])) 1 else dim(inpc[['migMpred']])[2]
			migFntraj <- if(is.null(inpc[['migFpred']])) 1 else dim(inpc[['migFpred']])[2]  
			migm <- array(0, dim=c(nmigcat, npredplus1, migMntraj), dimnames=list(ages[1:nmigcat], present.and.proj.years, NULL))
			migf <- array(0, dim=c(nmigcat, npredplus1, migFntraj), dimnames=list(ages[1:nmigcat], present.and.proj.years, NULL))
			# extend current mx to 130 age groups
			Kan.present <- KannistoAxBx.joint(inpc$MXm[,ncol(inpc$MXm), drop=FALSE], inpc$MXf[,ncol(inpc$MXf), drop=FALSE], 
			                                  compute.AxBx=FALSE, annual = inp$annual)
		}
		debug <- FALSE
		#stop('')
		if(!fixed.mx) {
			MxKan <- runKannisto(inpc, inp$start.year, lc.for.all = inp$lc.for.all, npred=npred, annual = inp$annual) 
			mortcast.args <- .prepare.for.mortality.projection(pattern = inpc$MXpattern, mxKan = MxKan, 
			                                                   hiv.params = inpc$HIVparams, lc.for.all = inp$lc.for.all, 
			                                                   annual = inp$annual)
			if(verbose) cat(", mx via ", paste(c(mortcast.args$meth1, mortcast.args$meth2), collapse = ","))
		} else {
			MxKan <- runKannisto.noLC(inpc, annual = inp$annual)
			LTres <- survival.fromLT(npred, MxKan, annual = inp$annual, verbose=verbose, debug=debug)
		}
		#npasfr <- nrow(inpc$PASFR)
		if(keep.vital.events) 
		    observed <- compute.observedVE(inpc, inp$pop.matrix, inpc$MIGtype, MxKan, country, 
															estim.years=inp$estim.years, annual = inp$annual)
		tfr.med <- apply(inpc$TFRpred, 1, median, na.rm = TRUE)[nrow(inpc$TFRpred)]
		for(itraj in 1:nr.traj) {
			if(any(is.na(inpc$TFRpred[,itraj]))) next
			if(!fixed.pasfr) 
				pasfr <- kantorova.pasfr(c(inpc$observed$TFRpred, inpc$TFRpred[,itraj]), inpc, 
										norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med, annual = inp$annual)
			else pasfr <- inpc$PASFR/100.
			asfr <- pasfr
			for(i in 1:nrow(pasfr)) asfr[i,] <- inpc$TFRpred[,itraj] * asfr[i,]
			if(!fixed.mx) {
			    LTres <- project.mortality(inpc$e0Mpred[,itraj], inpc$e0Fpred[,itraj], npred, mortcast.args = mortcast.args,
									        annual = inp$annual, verbose = verbose, debug = debug)
			}
			migpred <- list(M=NULL, F=NULL)
			for(sex in c('M', 'F')) {
				par <- paste0('mig', sex, 'pred')
				migpred[[sex]] <- as.matrix(if(is.null(inpc[[par]])) inpc[[paste0('MIG', tolower(sex))]] else inpc[[par]][,itraj,])
			}
			popres <- StoPopProj(npred, inpc, LTres, asfr, migpred, inpc$MIGtype, country.name=UNlocations[country.idx,'name'],
									keep.vital.events=keep.vital.events, annual = inp$annual)
			totp[,itraj] <- popres$totpop
			totpm[,,itraj] <- popres$mpop
			totpf[,,itraj] <- popres$fpop
			if(keep.vital.events) {
			    migtrajm <- min(itraj, migMntraj)
			    migtrajf <- min(itraj, migFntraj)
			    if(!is.null(observed)) {
			        btm[1:dim(observed$btm)[1],1,itraj] <- observed$btm[,dim(observed$btm)[2],]
			        btf[1:dim(observed$btf)[1],1,itraj] <- observed$btf[,dim(observed$btf)[2],]
			        deathsm[1:dim(observed$deathsm)[1],1,itraj] <- observed$deathsm[,dim(observed$deathsm)[2],]
			        deathsf[1:dim(observed$deathsf)[1],1,itraj] <- observed$deathsf[,dim(observed$deathsf)[2],]
			        asfert[1:dim(observed$asfert)[1],1,itraj] <- observed$asfert[,dim(observed$asfert)[2],]
			        pasfert[1:dim(pasfr)[1],1,itraj] <- inpc$observed$PASFR[,dim(inpc$observed$PASFR)[2]]
			        migm[1:dim(inpc$observed$MIGm)[1],1,migtrajm] <- inpc$observed$MIGm[,dim(inpc$observed$MIGm)[2]]
			        migf[1:dim(inpc$observed$MIGf)[1],1,migtrajf] <- inpc$observed$MIGf[,dim(inpc$observed$MIGf)[2]]
			    }
				btm[,2:npredplus1,itraj] <- popres$mbt
				btf[,2:npredplus1,itraj] <- popres$fbt
				deathsm[,2:npredplus1,itraj] <- popres$mdeaths
				deathsf[,2:npredplus1,itraj] <- popres$fdeaths
				asfert[,2:npredplus1,itraj] <- asfr
				pasfert[,2:npredplus1,itraj] <- pasfr*100
				mxm[,2:npredplus1,itraj] <- LTres$mx[[1]]
				mxm[1:length(Kan.present$male$mx),1,itraj] <- Kan.present$male$mx
				mxf[,2:npredplus1,itraj] <- LTres$mx[[2]]
				mxf[1:length(Kan.present$female$mx),1,itraj] <- Kan.present$female$mx
				migm[,2:npredplus1,migtrajm] <- migpred[['M']]
				migf[,2:npredplus1,migtrajf] <- migpred[['F']]
			}
		}
		for (variant in 1:nvariants) { # compute the two half child variants
			if(!fixed.pasfr) 
				pasfr <- kantorova.pasfr(c(inpc$observed$TFRpred, inpc$TFRhalfchild[variant,]), inpc, 
										norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med, annual = inp$annual)
			else pasfr <- inpc$PASFR/100.
			asfr <- pasfr
			for(i in 1:nrow(pasfr)) asfr[i,] <- inpc$TFRhalfchild[variant,] * asfr[i,]
			if(!fixed.mx) LTres <- project.mortality(inpc$e0Mmedian, inpc$e0Fmedian, npred, mortcast.args = mortcast.args,
			                                         annual = inp$annual, verbose=verbose, debug=debug)
			migpred.hch <- list(M=NULL, F=NULL)
			for(sex in c('M', 'F')) {
				par <- paste0('mig', sex, 'median')
				migpred.hch[[sex]] <- as.matrix(if(is.null(inpc[[par]])) inpc[[paste0('MIG', tolower(sex))]] else inpc[[par]])
			}
			popres <- StoPopProj(npred, inpc, LTres, asfr, migpred.hch, inpc$MIGtype, 
							country.name=UNlocations[country.idx,'name'], 
							keep.vital.events=keep.vital.events, annual = inp$annual)
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
				pasfert.hch[,2:npredplus1,variant] <- pasfr*100
				pasfert.hch[,1,variant] <- pasfert[,1,1]
				mxm.hch[,2:npredplus1,variant] <- LTres$mx[[1]]
				mxf.hch[,2:npredplus1,variant] <- LTres$mx[[2]]
				mxm.hch[,1,variant] <- mxm[,1,1]
				mxf.hch[,1,variant] <- mxf[,1,1]
			}
		}
		trajectory.indices <- inpc$trajectory.indices
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch, trajectory.indices,
			 file = file.path(outdir, paste0('totpop_country', country, '.rda')))
		if(keep.vital.events) 
			save(btm, btf, deathsm, deathsf, asfert, pasfert, mxm, mxf, migm, migf,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, pasfert.hch, 
				mxm.hch, mxf.hch, 
				observed,
					file=file.path(outdir, paste0('vital_events_country', country, '.rda')))
		
		res <- list()
		within(res, {
		    PIs <- apply(totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		    dimnames(PIs)[1:2] <- list(quantiles.to.keep, present.and.proj.years.pop)
		    means <- abind(apply(totp, 1, mean, na.rm = TRUE),
		                    apply(totp, 1, sd, na.rm = TRUE), along = 0, 
		                   new.names = c("mean", "sd"))
		    qMage <- qFage <- qPropMage <- qPropFage <- array(NA, c(nages, nquant, nr_project+1),
		          dimnames=list(ages, quantiles.to.keep, present.and.proj.years.pop))
		    for (i in 1:nages) {
		        if(nr.traj == 1) {
		            qMage[i,,] <- matrix(rep(totpm[i,,1],nquant) , nrow=nquant, byrow=TRUE)
		            qFage[i,,] <- matrix(rep(totpf[i,,1],nquant) , nrow=nquant, byrow=TRUE)
		            qPropMage[i,,] <- matrix(rep(totpm[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
		            qPropFage[i,,] <- matrix(rep(totpf[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
		        } else {
		            qMage[i,,] <- apply(totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
		            qFage[i,,] <- apply(totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
		            qPropMage[i,,] <- apply(totpm[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		            qPropFage[i,,] <- apply(totpf[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		        }
		    }
		    stotpm <- colSums(totpm, na.rm=TRUE)
		    qM <- apply(stotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		    meansM <- abind(apply(stotpm, 1, mean, na.rm = TRUE), 
		                    apply(stotpm, 1, sd, na.rm = TRUE), along = 0, new.names = c("mean", "sd"))
		    stotpf <- colSums(totpf, na.rm=TRUE)
		    qF <- apply(stotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		    meansF <- abind(apply(stotpf, 1, mean, na.rm = TRUE),
		                    apply(stotpf, 1, sd, na.rm = TRUE), along = 0, new.names = c("mean", "sd"))
		    dimnames(qM)[1:2] <- dimnames(qF)[1:2] <- list(quantiles.to.keep, present.and.proj.years.pop)
		    migr.modified <- migr.modified
		    nr.traj <- nr.traj
		    stotp <- stotpf <- NULL
		})
    } # end of predict.one.country
    
    update.results <- function(cidx, res, bayesPop.prediction) {
        if(length(cidx) == 1 && length(res) > 1) res <- list(res)
        migr.modified <- FALSE
        remove <- c()
        for(i in seq_along(cidx)) {
            idx <- cidx[i]
            rs <- res[[i]]
            country <- country.codes[idx]
            idx.in.pred.overwrite <- which(bayesPop.prediction$countries[,'code'] == country)
            
            if(is.null(rs)) {
                remove <- c(remove, country)
                next
            }

            if(length(idx.in.pred.overwrite)>0) {
                bayesPop.prediction$quantiles[idx.in.pred.overwrite,,] <- rs$PIs
                bayesPop.prediction$traj.mean.sd[idx.in.pred.overwrite,,] <- rs$means
                bayesPop.prediction$traj.mean.sdM[idx.in.pred.overwrite,,] <- rs$meansM
                bayesPop.prediction$traj.mean.sdF[idx.in.pred.overwrite,,] <- rs$meansF
                bayesPop.prediction$quantilesM[idx.in.pred.overwrite,,] <- rs$qM
                bayesPop.prediction$quantilesF[idx.in.pred.overwrite,,] <- rs$qF
                bayesPop.prediction$quantilesMage[idx.in.pred.overwrite,,,] <- rs$qMage
                bayesPop.prediction$quantilesFage[idx.in.pred.overwrite,,,] <- rs$qFage
                bayesPop.prediction$quantilesPropMage[idx.in.pred.overwrite,,,] <- rs$qPropMage
                bayesPop.prediction$quantilesPropFage[idx.in.pred.overwrite,,,] <- rs$qPropFage
            } else { 
                country.row <- UNlocations[countries.idx[idx],c('country_code', 'name')]
                colnames(country.row) <- c('code', 'name')
                #stop("")
                new.names <- as.character(c(bayesPop.prediction$countries$code, country))
                bayesPop.prediction$quantiles <- abind(bayesPop.prediction$quantiles, rs$PIs, 
                                                       along = 1, new.names = list(new.names, NULL, NULL))
                bayesPop.prediction$traj.mean.sd <- abind(bayesPop.prediction$traj.mean.sd, rs$means, 
                                                          along = 1, new.names = list(new.names, NULL, NULL))
                bayesPop.prediction$traj.mean.sdM <- abind(bayesPop.prediction$traj.mean.sdM, rs$meansM, 
                                                           along = 1, new.names = list(new.names, NULL, NULL))
                bayesPop.prediction$traj.mean.sdF <- abind(bayesPop.prediction$traj.mean.sdF, rs$meansF, 
                                                           along = 1, new.names = list(new.names, NULL, NULL))
                bayesPop.prediction$quantilesM <- abind(bayesPop.prediction$quantilesM, rs$qM, 
                                                        along = 1, new.names = list(new.names, NULL, NULL))
                bayesPop.prediction$quantilesF <- abind(bayesPop.prediction$quantilesF, rs$qF, 
                                                        along = 1, new.names = list(new.names, NULL, NULL))
                bayesPop.prediction$quantilesMage <- abind(bayesPop.prediction$quantilesMage, rs$qMage, 
                                                           along = 1, new.names = list(new.names, NULL, NULL, NULL))
                bayesPop.prediction$quantilesFage <- abind(bayesPop.prediction$quantilesFage, rs$qFage, 
                                                           along = 1, new.names = list(new.names, NULL, NULL, NULL))
                bayesPop.prediction$quantilesPropMage <- abind(bayesPop.prediction$quantilesPropMage, rs$qPropMage, 
                                                               along = 1, new.names = list(new.names, NULL, NULL, NULL))
                bayesPop.prediction$quantilesPropFage <- abind(bayesPop.prediction$quantilesPropFage, rs$qPropFage, 
                                                               along = 1, new.names = list(new.names, NULL, NULL, NULL))
                bayesPop.prediction$countries <- rbind(bayesPop.prediction$countries, country.row)
            }
            if(rs$migr.modified) migr.modified <- TRUE
            if(bayesPop.prediction$nr.traj != rs$nr.traj) bayesPop.prediction$nr.traj <- rs$nr.traj
        }
        if(migr.modified) {
            for(what.mig in c('MIGm', 'MIGf')) {
                inp.to.save[[what.mig]] <- inp[[what.mig]]
                inp.to.save$observed[[what.mig]] <- inp$observed[[what.mig]]
            }
            bayesPop.prediction$inputs <- inp.to.save
        }
        if(length(remove) > 0) {
            # remove countries from bayesPop.prediction if present
            idx.in.pred <- which(bayesPop.prediction$countries[,'code'] %in% remove)
            if(length(idx.in.pred) > 0) {
                for(item in c("quantiles", "traj.mean.sd", "traj.mean.sdM", "traj.mean.sdF", 
                              "quantilesM", "quantilesF")) 
                    bayesPop.prediction[[item]] <- bayesPop.prediction[[item]][-idx.in.pred,,, drop = FALSE]
                for(item in c("quantilesMage", "quantilesFage", "quantilesPropMage", "quantilesPropFage"))
                    bayesPop.prediction[[item]] <- bayesPop.prediction[[item]][-idx.in.pred,,,, drop = FALSE]
                bayesPop.prediction$countries <- bayesPop.prediction$countries[-idx.in.pred,, drop = FALSE]
            }
        }
        save(bayesPop.prediction, file=prediction.file)
        return(bayesPop.prediction)
    } # end of updating result
    
    cntries.table <- UNlocations[countries.idx,c('country_code', 'name')]
    colnames(cntries.table)[1] <- "code"
    bayesPop.prediction <- if(!is.null(pred)) .cleanup.pop.before.save(pred, remove.cache= any(country.codes %in% pred$countries[,'code'])) 
        else structure(list(
            nr.traj = nr.traj,	
            # assign empty arrays
            quantiles = PIs_cqp,
            traj.mean.sd = mean_sd,
            quantilesM = quantM, 
            traj.mean.sdM = mean_sdM,
            quantilesF = quantF, 
            traj.mean.sdF = mean_sdF,
            quantilesMage = quantMage, 
            quantilesFage = quantFage, 
            quantilesPropMage = quantPropMage, 
            quantilesPropFage = quantPropFage,
            estim.years=inp$estim.years, 
            proj.years=present.and.proj.years, # includes present period (middle of periods)
            proj.years.pop=present.and.proj.years.pop, # end of periods
            wpp.year = inp$wpp.year,
            inputs = inp.to.save, # save as list because environment takes much more space
            function.inputs=function.inputs,
            countries=cntries.table,
            ages=ages,
            annual = inp$annual), class='bayesPop.prediction')
    
    if(parallel) {
        cl <- create.pop.cluster(nr.nodes, ...)
        clusterExport(cl, exporting.objects, envir=environment())
        cntry.res <- clusterApplyLB(cl, 1:ncountries, predict.one.country, nr.traj = nr.traj, nr_project = nr_project)
        stopCluster(cl)
        bayesPop.prediction <- update.results(1:ncountries, cntry.res, bayesPop.prediction)
    } else {
        for(cidx in 1:ncountries) {
            cntry.res <- predict.one.country(cidx, nr.traj, nr_project)
            bayesPop.prediction <- update.results(cidx, cntry.res, bayesPop.prediction)
        }
	} 
	cat('\nPrediction stored into', outdir, '\n')
	return(bayesPop.prediction)
}

read.pop.file <- function(file, ...) 
	return(read.delim(file=file, comment.char='#', check.names=FALSE, ...))
	
load.wpp.dataset <- function(...)
	bayesTFR:::load.bdem.dataset(...)

load.inputs <- function(inputs, start.year, present.year, end.year, wpp.year, fixed.mx=FALSE, 
                        fixed.pasfr=FALSE, all.countries=TRUE, existing.mig=NULL, 
                        lc.for.hiv = TRUE, lc.for.all = FALSE, annual = FALSE, verbose=FALSE) {
	observed <- list()
	pop.ini.matrix <- pop.ini <- GQ <- list(M=NULL, F=NULL)
	# Get initial population counts
	for(sex in c('M', 'F')) {
		dataset.name <- paste0('pop', sex)
		if(is.null(inputs[[dataset.name]])) 
			POP0 <- load.wpp.dataset(dataset.name, wpp.year)	
		else POP0 <- read.pop.file(inputs[[dataset.name]])
		num.columns <- grep('^[0-9]{4}$', colnames(POP0), value=TRUE) # values of year-columns
		if(!is.element(as.character(present.year), colnames(POP0))) {
			stop('Wrong present.year. ', present.year, ' not available in the ', dataset.name, ' file.\nAvailable years: ',
				paste(num.columns, collapse=', '))
		}
		POP0[is.na(POP0)] <- 0
		if(! 'age' %in% colnames(POP0) || ! 'country_code' %in% colnames(POP0))
		    stop('Columns "age" and "country_code" must be present in the population datasets.')
		num.columns <- num.columns[which(as.integer(num.columns)<= present.year)]
		pop.ini.matrix[[sex]] <- POP0[,num.columns, drop=FALSE]
		dimnames(pop.ini.matrix[[sex]]) <- list(paste(POP0[,'country_code'], POP0[,'age'], sep='_'), 
									as.character(as.integer(num.columns)))
		pop.ini[[sex]] <- POP0[,c('country_code', 'age', present.year)]
		dataset.name <- paste0('GQpop', sex)
		if(!is.null(inputs[[dataset.name]])) {
		    GQ[[sex]] <- read.pop.file(inputs[[dataset.name]])
		    colnames(GQ[[sex]]) <- tolower(colnames(GQ[[sex]]))
		    if(! 'age' %in% colnames(GQ[[sex]]) || ! 'country_code' %in% colnames(GQ[[sex]]) || ! 'gq' %in% colnames(GQ[[sex]]))
		        stop('Columns "age", "country_code" and "gq" must be present in the GQpop datasets.')
		    GQ[[sex]] <- GQ[[sex]][, c("country_code", "age", "gq")]
		}
	}
	POPm0 <- pop.ini[['M']]
	POPf0 <- pop.ini[['F']]
	GQm <- GQ[['M']]
	GQf <- GQ[['F']]

	# Get death rates
	MXm.pred <- MXf.pred <- NULL
	if(is.null(inputs$mxM)) 
		MXm <- load.wpp.dataset('mxM', wpp.year)
	else MXm <- read.pop.file(inputs$mxM)
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
	estim.periods <- names.MXm.data[num.columns[1:present.index]]
	start.year <- cols.starty[start.index]
	
	if(fixed.mx) {
		end.index <- if(annual) which(cols.endy == end.year) else which((cols.endy >= end.year) & (cols.starty < end.year))
		proj.periods <- names.MXm.data[num.columns[(present.index+1):end.index]]
		MXm.pred <- MXm[,c('country_code', 'age', proj.periods)]
	}
	MXm <- MXm[,c('country_code', 'age', estim.periods)]
	if(is.null(inputs$mxF)) 
		MXf <- load.wpp.dataset('mxF', wpp.year)
	else MXf <- read.pop.file(inputs$mxF)
	if(fixed.mx) MXf.pred <- MXf[,c('country_code', 'age', proj.periods)]
	MXf <- MXf[,c('country_code', 'age', estim.periods)]
	
	estim.years <- cols.starty[start.index:present.index]
	if(!annual) estim.years <- estim.years + 3

	# Get sex ratio at birth
	srblist <- .get.srb.data.and.time.periods(inputs$srb, present.year, end.year, wpp.year, annual = annual)
	SRB <- srblist$srb
	observed$SRB <- srblist$obs.srb
	proj.periods <- srblist$proj.periods
	obs.periods <- srblist$obs.periods
	proj.years <- srblist$proj.years

	# Get percentage age-specific fertility rate
	pasfrlist <- .get.pasfr.data(inputs$pasfr, wpp.year, obs.periods, proj.periods, 
	                             include.projection=fixed.pasfr)
	PASFR <- pasfrlist$pasfr
	observed$PASFR <- pasfrlist$obs.pasfr
	
	# Get migration type, migration base year, mx & pasfr patterns
	patterns <- .get.mig.mx.pasfr.patterns(inputs, wpp.year, lc.for.hiv = lc.for.hiv)
	MIGtype <- patterns$mig.type
	MXpattern <- patterns$mx.pattern
	PASFRpattern <- patterns$pasfr.pattern
	
	# Get HIV parameters to be used with hiv.mortmod()
	HIVparams <- .get.hiv.params(inputs)

	# Get age-specific migration
	miginp <- .get.mig.data(inputs, wpp.year, annual, periods = c(estim.periods, proj.periods), 
	                        existing.mig = existing.mig, all.countries = all.countries, pop0 = POPm0,
	                        verbose = verbose)

	MIGm <- miginp[["migM"]]
	MIGf <- miginp[["migF"]]
	
	if(!is.null(obs.periods)) {
		if (!is.null(existing.mig)) { # Migration dataset already exists, e.g. from a priovous simulation for different country
			observed$MIGm <- existing.mig$obsMIGm
			observed$MIGf <- existing.mig$obsMIGf
		} else {
			avail.obs.periods <- is.element(obs.periods, colnames(MIGm))
			observed$MIGm <- MIGm[,c('country_code', 'age', obs.periods[avail.obs.periods])]
			observed$MIGf <- MIGf[,c('country_code', 'age', obs.periods[avail.obs.periods])]
		}
	}
	MIGm <- MIGm[,c('country_code', 'age', proj.periods)]
	MIGf <- MIGf[,c('country_code', 'age', proj.periods)]
	# Get migration trajectories if available
	migpr <- .load.mig.traj(inputs, verbose = verbose)
	migMpred <- migpr$M
	migFpred <- migpr$F
	
	# Get life expectancy
	e0F.wpp.median.loaded <- FALSE
	e0Fpred <- e0Mpred <- NULL
	if(!fixed.mx){
	if(!is.null(inputs$e0F.file)) { # female
		if(inputs$e0F.file == 'median_') {
			e0Fpred <- .load.wpp.traj('e0F', wpp.year, median.only=TRUE)
			e0F.wpp.median.loaded <- TRUE
		} else {
			file.name <-  inputs$e0F.file
			if(!file.exists(file.name))
				stop('File ', file.name, 
					' does not exist.\nSet e0F.sim.dir, e0F.file or change WPP year.')
		 	# comma separated trajectories file
		 	if(verbose) cat('\nLoading ', file.name)
			e0Fpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
			e0Fpred <- e0Fpred[,c('LocID', 'Year', 'Trajectory', 'e0')]
			colnames(e0Fpred) <- c('country_code', 'year', 'trajectory', 'value')
		}
	} else {
		if(!is.null(inputs$e0F.sim.dir)) { 
			if(inputs$e0F.sim.dir == 'median_') {
				e0Fpred <- .load.wpp.traj('e0F', wpp.year, median.only=TRUE)
				e0F.wpp.median.loaded <- TRUE
			} else 
				e0Fpred <- get.e0.prediction(inputs$e0F.sim.dir, mcmc.dir=NA)
		} else e0Fpred <- .load.wpp.traj('e0F', wpp.year)			
	}
	
	if(!is.null(inputs$e0M.file)) { # male
		if(inputs$e0M.file == 'median_')
			e0Mpred <- .load.wpp.traj('e0M', wpp.year, median.only=TRUE)
		else {
			file.name <-  inputs$e0M.file
			if(!file.exists(file.name)) 
				stop('File ', file.name, 
					' does not exist.\nSet e0M.sim.dir, e0M.file or change WPP year.')
			if(verbose) cat('\nLoading ', file.name)
			e0Mpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
			e0Mpred <- e0Mpred[,c('LocID', 'Year', 'Trajectory', 'e0')]
			colnames(e0Mpred) <- c('country_code', 'year', 'trajectory', 'value')
		}
	} else {
		if(!is.null(inputs$e0M.sim.dir)) { 
			if(inputs$e0M.sim.dir == 'joint_') {
				if(e0F.wpp.median.loaded) e0Mpred <- .load.wpp.traj('e0M', wpp.year)
				else {
					if(!has.e0.jmale.prediction(e0Fpred))
						stop('No joint prediction for female and male available. Correct the e0M.sim.dir argument.' )
					e0Mpred <- get.e0.jmale.prediction(e0Fpred)
				}
			} else e0Mpred <- get.e0.prediction(inputs$e0M.sim.dir, mcmc.dir=NA) # independent from female
		} else
			e0Mpred <- .load.wpp.traj('e0M', wpp.year)
	}
	} # end if(!fixed.mx)
	
	# Get TFR
	TFRpred <- .get.tfr.data(inputs, wpp.year, verbose=verbose)
	
	inp <- new.env()
	for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MXm.pred', 'MXf.pred', 'MXpattern', 'SRB',
				'PASFR', 'PASFRpattern', 'MIGtype', 'MIGm', 'MIGf', 'HIVparams', 'GQm', 'GQf',
				'e0Mpred', 'e0Fpred', 'TFRpred', 'migMpred', 'migFpred', 'estim.years', 'proj.years', 'wpp.year', 
				'start.year', 'present.year', 'end.year', 'annual', 'fixed.mx', 'fixed.pasfr', 
				'lc.for.hiv', 'lc.for.all', 'observed'))
		assign(par, get(par), envir=inp)
	inp$pop.matrix <- list(male=pop.ini.matrix[['M']], female=pop.ini.matrix[['F']])
	inp$PASFRnorms <- compute.pasfr.global.norms(inp)
	inp$average.annual <- inputs$average.annual
	return(inp)
}

.get.srb.data.and.time.periods <- function(srb.data, present.year, end.year, wpp.year, annual = FALSE) {
    if(is.null(srb.data)) 
        SRB <- load.wpp.dataset('sexRatio', wpp.year)
    else {
        if(is.character(srb.data)) # file name
            SRB <- read.pop.file(srb.data)
        else SRB <- srb.data
    }
    names.SRB.data <- names(SRB)
    numcol.expr <- if(annual) '^[0-9]{4}$' else '^[0-9]{4}.[0-9]{4}$'
    num.columns <- grep(numcol.expr, names.SRB.data) # index of year-columns
    cols.starty <- as.integer(substr(names.SRB.data[num.columns], 1,4))
    if(!annual) {
        cols.endy <-  as.integer(substr(names.SRB.data[num.columns], 6,9))
        start.index <- which((cols.starty <= present.year) & (cols.endy > present.year))
        end.index <- which((cols.endy >= end.year) & (cols.starty < end.year))
    } else {
        start.index <- which(cols.starty > present.year)[1]
        end.index <- which(cols.starty == end.year)
    }
    
    if(length(end.index) == 0) {
        end.index <- length(num.columns)
        warning('Data for SexRatioAtBirth not available for all projection periods.\nLast projection period set to ', 
                names.SRB.data[num.columns[end.index]])
    }
    proj.periods <- names.SRB.data[num.columns[start.index:end.index]]
    obs.periods <- NULL
    obs.SRB <- NULL
    if(start.index > 1) {
        obs.periods <- names.SRB.data[num.columns[1:(start.index-1)]]
        obs.SRB <- SRB[,c('country_code', obs.periods)]
    }
    SRB <- SRB[,c('country_code', proj.periods)]
    proj.years <- cols.starty[start.index:end.index]
    if(!annual) proj.years <- proj.years + 3
    return(list(srb=SRB, obs.srb=obs.SRB, proj.periods=proj.periods, 
           obs.periods=obs.periods, proj.years=proj.years))
}

.get.pasfr.data <- function(pasfr.data, wpp.year, obs.periods, proj.periods,
                            include.projection=TRUE) {
    if(is.null(pasfr.data)) 
        PASFR <- load.wpp.dataset('percentASFR', wpp.year)
    else {
        if(is.character(pasfr.data)) # file name
            PASFR <- read.pop.file(pasfr.data)
        else PASFR <- pasfr.data
    }
    obs.PASFR <- NULL
    if(!is.null(obs.periods)) {
        avail.obs.periods <- is.element(obs.periods, colnames(PASFR))
        obs.PASFR <- PASFR[,c('country_code', 'age', obs.periods[avail.obs.periods])]
    }
    if(include.projection)
        PASFR <- PASFR[,c('country_code', 'age', proj.periods)]
    else PASFR <- NULL
    return(list(pasfr=PASFR, obs.pasfr=obs.PASFR))
}

.get.mig.template <- function(countries, ages, time.periods, template = NULL, id.col = "country_code"){
    if(!is.null(template)) return(template)
    nages <- length(ages)
    df <- data.table(id = rep(countries, each = nages), age = rep(ages, length(countries)))
    df <- cbind(df, matrix(0, nrow = nrow(df), ncol = length(time.periods),
                           dimnames = list(NULL, time.periods)))
    setnames(df, "id", id.col)
    df
}

migration.totals2age <- function(df, ages = NULL, annual = FALSE, time.periods = NULL, 
                                 schedule = NULL, scale = 1, id.col = "country_code", ...) {
    mig <- totmig <- rc <- NULL
    if(is.null(dim(df))) df <- t(df)
    if(!is.data.table(df)) df <- data.table(df)
    if(is.null(time.periods)) {
        if(is.null(colnames(df))) {
            colnames(df) <- seq_len(ncol(df))
            time.periods <- colnames(df)
        } else {
            numcol.expr <- if(annual) '^[0-9]{4}$' else '^[0-9]{4}.[0-9]{4}$'
            time.periods <- colnames(df)[grep(numcol.expr, colnames(df))]
        }
    } else 
        time.periods <- intersect(colnames(df), as.character(time.periods))
    
    if(is.null(ages))
        ages <- get.age.labels(all.age.index(annual, observed = TRUE), 
                               age.is.index=TRUE, single.year = annual, last.open = TRUE)
        
    if(length(time.periods) == 0) stop("No time periods detected.")
    cntry.missing <- FALSE
    if(is.null(df[[id.col]])) {
        df[, (id.col) := seq_len(nrow(df))]
        cntry.missing <- TRUE
    }
    migtempl <- .get.mig.template(unique(df[[id.col]]), ages, time.periods, id.col = id.col, ...) 
    if(length(ages) != length(unique(migtempl$age))) stop("Argument ages does not correspond to age in template.") 
    
    totmigl <- melt(df[, c(id.col, time.periods), with = FALSE], id.vars = id.col,
                    variable.name = "year", value.name = "totmig")
    migtempll <- melt(migtempl[, c(id.col, "age", time.periods), with = FALSE], 
                      id.vars = c(id.col, "age"), variable.name = "year", value.name = "mig")
    age.idx <- 1:length(ages)
    if(is.null(schedule)) schedule <- rcastro.schedule(annual)
    rcdf <- data.table(age = migtempl[["age"]][age.idx],
                       age.idx = age.idx,
                       rc = scale * schedule[age.idx]/sum(schedule[age.idx]))
    migtmp <- merge(merge(migtempll, totmigl, by = c(id.col, "year"), sort = FALSE), 
                    rcdf, by = "age", sort = FALSE)
    migtmp[, mig := totmig * rc]
    frm <- paste(id.col, "+ age.idx + age ~ year")
    res <- dcast(migtmp, frm, value.var = "mig")
    res[["age.idx"]] <- NULL
    if(cntry.missing) res[[id.col]] <- NULL
    return(res)
}

.get.mig.data <- function(inputs, wpp.year, annual, periods, existing.mig = NULL, 
                          all.countries = TRUE, pop0 = NULL, verbose = FALSE) {
    # Get age-specific migration
    wppds <- data(package=paste0('wpp', wpp.year))
    recon.mig <- NULL
    miginp <- list()
    migtempl <- NULL
    for(sex in c("M", "F")) {
        inpname <- paste0('mig', sex) 
        if(!is.null(inputs[[inpname]])) { # migration given by sex and age
            miginp[[inpname]] <- read.pop.file(inputs[[inpname]])
            next
        }
        # create a template  to be filled with derived migration
        if(is.null(migtempl)) {
            if(!all.countries && !is.null(existing.mig)) { # simulation is run for a subset of countries
                # AND migration dataset already exists, e.g. from a previous simulation for different country
                migtempl <- existing.mig[[paste0('MIG', tolower(sex))]]
            } else {
                # Here create only a dataframe filled with NAs 
                if(!annual)
                    migtempl <- load.wpp.dataset(paste0('migration', sex), 2012) # structure is taken from the wpp2012 dataset
                else {
                    # use countries and ages from the population dataset
                    #migtempl <- cbind(pop0[, c("country_code", "age")], 
                    #                  matrix(0, nrow = nrow(pop0), ncol = length(periods), 
                    #                         dimnames = list(NULL, periods)))
                    migtempl <- .get.mig.template(unique(pop0$country_code), 
                                                  ages = pop0$age[all.age.index(annual, observed = TRUE)],
                                                  time.periods = periods)
                }
                migtempl[,which(!colnames(migtempl) %in% c("country", "country_code", "age"))] <- NA
            }
            migtempl <- data.table(migtempl)
        }
        fname <- paste0(inpname, 't')
        fnametot <- paste0("mig")
        if(!is.null(inputs[[fname]]) || !is.null(inputs[[fnametot]])) { # migration given as totals or totals by sex
            if(!is.null(inputs[[fname]]))
                totmig <- data.table(read.pop.file(inputs[[fname]]))
            else
                totmig <- data.table(read.pop.file(inputs[[fnametot]]))
            migcols <- intersect(colnames(totmig), periods)
            miginp[[inpname]] <- data.frame(migration.totals2age(totmig, ages = migtempl$age[all.age.index(annual, observed = TRUE)],
                                                                 annual = annual, time.periods = migcols, 
                                                                 scale = if(is.null(inputs[[fname]])) 0.5 else 1, # since the totals are sums over sexes
                                                                 template = migtempl), check.names = FALSE)
            next
        }
        # if migration is not given load default datasets
        if(annual) stop("Migration must be given.")
        if(paste0('migration', sex) %in% wppds$results[,'Item']) { # if available in the WPP package
            miginp[[inpname]] <- load.wpp.dataset(paste0('migration', sex), wpp.year)
            next
        }
        if(all.countries) { # reconstruct migration for all countries
            if(is.null(recon.mig)) recon.mig <- age.specific.migration(wpp.year = wpp.year, verbose = verbose)
            miginp[[inpname]] <- recon.mig[[list(M = "male", F = "female")[[sex]]]]
            next
        }
        # Projection for a subset of countries, therefore migration will be recontructed later (in get.country.inputs).
        miginp[[inpname]] <- data.frame(migtempl, check.names = FALSE)
    }
    return(miginp)
}


.get.hiv.params <- function(inputs){
    if(is.null(inputs$hiv.params)) return(NULL)
    params <- read.pop.file(inputs$hiv.params, stringsAsFactors = FALSE)
    # remove country name
    if(any(colnames(params) %in% c("country", "name")))
        params <- params[, -which(colnames(params) %in% c("country", "name"))]
    for(par in c("param"))
        if(! par %in% colnames(params)) stop("Column ", par, " is obligatory in the hiv.params file.")
    if(! "sex" %in% colnames(params)) {
        warning("Column 'sex' is missing in the hiv.params file. The same values will be used for female and male.")
        params <- rbind(cbind(params, sex = "male"), cbind(params, sex = "female"))
    }
    # replace periods by mid-years in column names if needed
    cnames <- colnames(params)
    num.columns <- grep('^[0-9]{4}.[0-9]{4}$', cnames) # index of year-columns
    if(length(num.columns) > 0) {
        cols.start <- as.integer(substr(cnames[num.columns], 1,4))
        colnames(params)[num.columns] <- cols.start + 3
    }
    return(params)
}

.get.mig.mx.pasfr.patterns <- function(inputs, wpp.year, pattern.data = NULL, lc.for.hiv = TRUE) {
    if(is.null(pattern.data)) {
        pattern.file <- if(!is.null(inputs$patterns)) inputs$patterns else inputs$mig.type
        if(is.null(pattern.file)) 
            vwBase <- get(paste0('vwBaseYear', wpp.year))
        else vwBase <- read.pop.file(pattern.file)
    } else vwBase <- pattern.data
    
    if("PasfrNorm" %in% colnames(vwBase) && !is.factor(vwBase$PasfrNorm))
        vwBase$PasfrNorm <- as.factor(vwBase$PasfrNorm)
    
    create.pattern <- function(dataset, columns, char.columns = c()) {
        pattern <- data.frame(dataset[,'country_code'])
        for(col in columns)
            if(col %in% colnames(dataset))
                pattern <- cbind(pattern, dataset[,col])
        for(col in char.columns)
            if(col %in% colnames(dataset))
                pattern <- cbind(pattern, as.character(dataset[,col]), stringsAsFactors = FALSE)
        if(ncol(pattern)==1) pattern <- NULL
        else colnames(pattern) <- c('country_code', c(columns, char.columns)[c(columns, char.columns) %in% colnames(dataset)])
        return(pattern)
    }
    MIGtype <- create.pattern(vwBase, c('ProjFirstYear', 'MigCode'))
    MXpattern <- create.pattern(vwBase, c("AgeMortProjAdjSR", "LatestAgeMortalityPattern", 
                                          "SmoothLatestAgeMortalityPattern", "WPPAIDS", "HIVregion"),
                                char.columns = c("AgeMortalityType", "AgeMortalityPattern", "AgeMortProjMethod1", "AgeMortProjMethod2",
                                                 "AgeMortProjPattern", "AgeMortProjMethodWeights"))
    if(lc.for.hiv) { # replace HIVmortmod with LC
        for(col in c("AgeMortProjMethod1", "AgeMortProjMethod2"))
            if(col %in% colnames(MXpattern) && "HIVmortmod" %in% MXpattern[[col]]) MXpattern[MXpattern[[col]] == "HIVmortmod", col] <- "LC"
    }
    if(! "HIVregion" %in% colnames(MXpattern) && "area_code" %in% colnames(UNlocations)) 
        MXpattern[["HIVregion"]] <- as.integer(UNlocations[match(MXpattern$country_code, UNlocations$country_code), "area_code"] == 903)
    
    PASFRpattern <- create.pattern(vwBase, c("PasfrNorm", paste0("Pasfr", .remove.all.spaces(levels(vwBase$PasfrNorm)))))
    return(list(mig.type=MIGtype, mx.pattern=MXpattern, pasfr.pattern=PASFRpattern))
}

.get.tfr.data <- function(inputs,  wpp.year, verbose=FALSE) {
    trajectory <- NULL
  if(!is.null(inputs$tfr.file)) {
    if(inputs$tfr.file == 'median_')
      TFRpred <- .load.wpp.traj('tfr', wpp.year, median.only=TRUE)
    else {
      file.name <- inputs$tfr.file
      if(!file.exists(file.name))
        stop('File ', file.name, 
             ' does not exist.\nSet tfr.sim.dir, tfr.file or change WPP year.')
      if(verbose) cat('\nLoading ', file.name, '\n')
      if(!is.null(inputs$tfr.file.type) && inputs$tfr.file.type == "w") { # file in tab-separated wide format
          TFRpred <- data.table(read.csv(file=file.name, comment.char='#', check.names=FALSE, sep = "\t"))
          if('LocID' %in% colnames(TFRpred))
              setnames(TFRpred, 'LocID', 'country_code')
          TFRpred <- melt(TFRpred, id.vars = "country_code", 
                          measure.vars = setdiff(colnames(TFRpred), c("country_code", "country_name", "name", "last.observed")),
                          variable.name = "year")
          TFRpred[, trajectory := 1]
          TFRpred <- as.data.frame(TFRpred)
      } else { # file in comma-separated long format 
        TFRpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
        TFRpred <- TFRpred[,c('LocID', 'Year', 'Trajectory', 'TF')]
        colnames(TFRpred) <- c('country_code', 'year', 'trajectory', 'value')
      } 
    } 
  } else {
    if(!is.null(inputs$tfr.sim.dir)) 
      TFRpred <- get.tfr.prediction(inputs$tfr.sim.dir, mcmc.dir=NA)
    else TFRpred <- .load.wpp.traj('tfr', wpp.year)
  }
    return(TFRpred)
}

.set.inp.migration.if.needed <- function(inputs, inpc, country) {
	# If the age-specific migration values in "inputs" are NA 
	# (e.g. because it was reconstructed later in the code), 
	# replace those with values from the country-specific inputs.
	migr.modified <- FALSE
	for(what.mig in c('MIGm', 'MIGf')) {
		cidx <- which(inputs[[what.mig]]$country_code==country)
		cols <- intersect(colnames(inputs[[what.mig]]), colnames(inpc[[what.mig]]))
		if(any(!is.na(inputs[[what.mig]][cidx,cols]))) next
		inputs[[what.mig]][cidx,cols] <- inpc[[what.mig]][,cols]
		cidx <- which(inputs$observed[[what.mig]]$country_code==country)
		cols <- intersect(colnames(inputs$observed[[what.mig]]), colnames(inpc$observed[[what.mig]]))
		inputs$observed[[what.mig]][cidx,cols] <- inpc$observed[[what.mig]][,cols]
		migr.modified <- TRUE
	}
	return(migr.modified)
}

.load.wpp.traj <- function(type, wpp.year, median.only=FALSE) {
	dataset.obs <- dataset.low <- dataset.high <- NA
	if(type %in% c('e0F', 'e0M')){
		if(wpp.year < 2010) stop('e0 projections not available for wpp 2008.')
		#if(wpp.year == 2010) 
		dataset <- paste(type, 'proj', sep='')
		#else { # wpp 2012
		#	dataset <- paste(type, 'projMed', sep='')
		#}
		dataset.obs <- type
	} else { # tfr
		if(wpp.year < 2010) dataset <- type
		else {
			dataset <- paste(type, 'projMed', sep='')
			dataset.obs <- type
			if(!median.only) {
				dataset.low <- paste(type, 'projLow', sep='')
				dataset.high <- paste(type, 'projHigh', sep='')
			}
		}
	}
	if(!is.na(dataset.obs)) 
		pred.obs <- bayesTFR:::load.bdem.dataset(dataset.obs, wpp.year)
		
	pred.all <- NULL
	itraj <- 1
	for(dataset.name in c(dataset, dataset.low, dataset.high)) {
		if(is.na(dataset.name)) next
		pred <- bayesTFR:::load.bdem.dataset(dataset.name, wpp.year)
		remove.cols <- which(colnames(pred) %in% c('name', 'country', 'last.observed'))
		pred <- pred[,-remove.cols]
		if(!is.na(dataset.obs)) {
			remove.cols <- which(colnames(pred.obs) %in% c('name', 'country', 'last.observed', 
								colnames(pred)[-which(colnames(pred)=='country_code')]))
			pred <- merge(pred.obs[,-remove.cols], pred, by='country_code')
		}
		ncols <- ncol(pred)
		nonnum.idx <- which(colnames(pred)=='country_code')
		cnames <- colnames(pred)[-nonnum.idx]
		colnames(pred)[-nonnum.idx] <- paste(type, cnames, sep='')
		pred.long <- reshape(pred, direction='long', varying=(1:ncols)[-nonnum.idx], v.names=type, times=cnames)
		pred.long <- cbind(pred.long, year=as.integer(substr(pred.long$time,1,4))+3, trajectory=itraj)
		pred.long <- pred.long[,c('country_code', 'year', 'trajectory', type)]
		pred.all <- rbind(pred.all, pred.long)
		itraj <- itraj + 1
	}
	colnames(pred.all) <- c('country_code', 'year', 'trajectory', 'value')
	return(pred.all)
}

.load.mig.traj <- function(inputs, verbose = FALSE) {
    migMpred <- migFpred <- NULL
    migtrajcols <- list(LocID = "country_code", Year = "year", Trajectory = "trajectory", Age = "age", Migration = "value")
    for(sex in c('M', 'F', '')) {
        if(is.null(inputs[[paste0('mig', sex, 'traj')]])) next
        file.name <- inputs[[paste0('mig', sex, 'traj')]]
        if(!file.exists(file.name))
            stop('File ', file.name, ' does not exist.')
        # comma separated trajectories file
        if(verbose) cat('\nLoading ', file.name)
        
        migpred.raw <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
        cols.to.keep <- intersect(names(migtrajcols), colnames(migpred.raw))
        if(length((miss <- setdiff(setdiff(names(migtrajcols), "Age"), cols.to.keep)))>0)
            stop("Columns ", paste(miss, collapse = ", "), "are missing from ", file.name)
        migpred <- migpred.raw[, cols.to.keep]
        colnames(migpred) <- unlist(migtrajcols[cols.to.keep])
        if(sex == '') { # total for both sexes in one file; split it into two objects
            for(sext in c('M', 'F')){
                var.name <- paste0('mig',sext, 'pred')
                if(!is.null(get(var.name))) next
                migpreds <- migpred
                migpreds$value <- migpreds$value/2
                assign(var.name, migpreds)
            }
        } else {
            var.name <- paste0('mig',sex, 'pred')
            assign(var.name, migpred)
        }
    }
    return(list(M = migMpred, F = migFpred))
}

.get.migration.traj <- function(pred, par, country) {
		cidx <- pred$inputs[[par]][,'country_code'] == country 
		idx <- cidx & is.element(pred$inputs[[par]][,'year'], pred$inputs$proj.years)
		if(sum(idx) == 0) return(NULL)
		migdf <- pred$inputs[[par]][idx,-1]
		utrajs <- sort(unique(migdf$trajectory))
		ntrajs <- length(utrajs)
		if(! "age" %in% colnames(migdf)){ # need to disaggregate into age-specific trajectories
		    dfw <- dcast(data.table(migdf), trajectory ~ year)
		    adf <- migration.totals2age(dfw, annual = pred$inputs$annual, time.periods = colnames(dfw)[-1],
		                                id.col = "trajectory")
		    migdf <- melt(adf, value.name = "value", variable.name = "year", id.vars = c("trajectory", "age"))
		}
		migdf$age <- gsub("^\\s+|\\s+$", "", migdf$age) # trim leading and trailing whitespace
		lyears <- length(pred$inputs$proj.years)
		lage <- all.age.length(pred$inputs$annual, observed = TRUE)
		sorted.df <- data.frame(year=rep(pred$inputs$proj.years, each=ntrajs*lage), trajectory=rep(rep(utrajs, each=lage), times=lyears),
									age = get.age.labels(all.ages(pred$inputs$annual, observed = TRUE), last.open=TRUE, single.year = pred$inputs$annual))
		# this is to get rows of the data frame in a particular order
		migdf <- merge(sorted.df, migdf, sort=FALSE)
		res <- array(migdf$value, dim=c(lage, ntrajs, lyears))
		dimnames(res) <- list(1:lage,  NULL, pred$inputs$proj.years)
		return(res)	
}

.pasfr.norm.name <- function(norms)
	return(paste0('Pasfr', .remove.all.spaces(norms)))
	
compute.pasfr.global.norms <- function(inputs) {
	pattern <- inputs$PASFRpattern
	if(is.null(pattern) || !('PasfrNorm' %in% colnames(pattern))) return(NULL)
	norms <- .pasfr.norm.name(levels(pattern$PasfrNorm))
	result <- list()
	for(norm in norms) {
		tpasfr <- NULL
		countries <- pattern$country_code[which(pattern[[norm]]==1)]
		if(length(countries) == 0) next
		ccounter <- rep(0, )
		for(country in countries) {
			pasfr <- .get.par.from.inputs('PASFR', inputs$observed, country)
			pasfr <- .fill.pasfr.ages(pasfr, fert.ages(inputs$annual), check.length.only = !inputs$annual)
			if(is.null(ccounter)) ccounter <- rep(0, ncol(pasfr)) # deals with missing years for some countries
            is.not.observed <- apply(pasfr, 2, function(x) any(is.na(x)))
            if(any(is.not.observed))  # fill NA with 0 
                pasfr[,is.not.observed] <- 0
            ccounter <- ccounter + !is.not.observed
			tpasfr <- if(is.null(tpasfr)) pasfr else tpasfr + pasfr
		}
		tpasfr <- tpasfr/matrix(ccounter, nrow = nrow(tpasfr), ncol = ncol(tpasfr), byrow = TRUE)
		result[[norm]] <- scale(tpasfr, center=FALSE, scale=colSums(tpasfr))*100
	}
	return(result)
}

kantorova.pasfr <- function(tfr, inputs, norms, proj.years, tfr.med, annual = FALSE, nr.est.points = if(annual) 15 else 3) {
	logit <- function(x) log(x/(1-x))
	inv.logit <- function(x) exp(x)/(1+exp(x))
	fac.mac.start <- fert.ages(annual)[1]
	if(annual) {
	    fac.mac.start <- fac.mac.start + 0.5
	    by <- 1
	} else {
	    fac.mac.start <- fac.mac.start + 2.5
	    by <- 5
	}
	compute.mac <- function(x) {
		factors <- seq(fac.mac.start, by=by, length=dim(x)[1])
		mac <- rep(0, dim(x)[2])
        for(iage in 1:dim(x)[1]) 
        	mac <- mac + x[iage,]*factors[iage]/100.
        return(mac)
	}
	update.by.mac <- function(x, sp3i) {
		if(sp3i >= dim(x)[2] || all(is.na(x))) return(x)
		phase3i <- seq(sp3i, dim(x)[2])
		mac <- compute.mac(x[,phase3i]*100)		
		mac.norm <- compute.mac(matrix(gnorm, ncol=1))
		maxi <- which.max(mac)
		if(mac[maxi] <= mac.norm)
			return(x)
		x[,phase3i[maxi:length(phase3i)]] <- x[,phase3i[maxi]]
		return(x)
	}
	pattern <- inputs$PASFRpattern
	min.value <- 1e-6
	pasfr.obs <- inputs$observed$PASFR
	
	years <- as.integer(names(tfr))
	if(length(years)==0)
		years <- sort(seq(proj.years[length(proj.years)], length=length(tfr), by=-by))
	lyears <- length(years)
	years.long <- c(years, seq(years[lyears]+by, by=by, length=75/by)) # up to 2175
	tobs <- lyears - length(proj.years)
	end.year <- years[lyears]
	end.phase2 <- bayesTFR:::find.lambda.for.one.country(tfr, lyears, annual = annual)
	start.phase3 <- end.phase2 + 1
	if(start.phase3 > lyears) { # Phase 3 not observed until the end of projection (Case 2)
		if(tfr[lyears] > 1.8) { # regress the last 20 years to approximate start of Phase 3
		    nrpoints <- 20/by - 1 # one point because of the way it is used
            tbeyond <- 50/by
			df <- data.frame(tfr=tfr[(lyears-nrpoints):lyears], time=(lyears-nrpoints):lyears)
			reg <- lm(tfr~time, df)
			if(reg$coefficients[2] < -1e-3) {# use only if it has negative slope 
				start.phase3 <- min(round((1.8-reg$coefficients[1])/reg$coefficients[2],0)+1, lyears + tbeyond)
			} else {
				start.phase3 <- lyears + tbeyond
			}			
		}
		#endT <- years.long[min(start.phase3+5, length(years.long))]
		endT <- years.long[start.phase3 + 25/by]
	} else { # Case 1
		smaller.than.median <- tfr[start.phase3:lyears] < tfr.med
		if(all(smaller.than.median)) { # t_u does not exist
			#endT <- years.long[max(start.phase3+10, tobs+10)]
			endT <- years.long[max(lyears, start.phase3 + 25/by)]
		} else { # t_u exists
			first.larger <- which(!smaller.than.median)[1] + start.phase3 - 1
			endT <- years.long[max(first.larger, tobs + 20/by)]
		}
	}
	#endT <- years.long[max(start.phase3+5, tobs+5)] # no upper bound
	
	# Trend towards global model pattern
	startTi <- which(years == proj.years[1])
	gnorm <- norms[[.pasfr.norm.name(
	    if(is.null(pattern)) "Global Norm" else pattern[,'PasfrNorm'])]]
	gnorm <- gnorm[, ncol(gnorm)] # global norm from the last time period 
	asfr1 <- asfr2 <- res.asfr <- matrix(0, nrow=length(gnorm), ncol=length(proj.years))
	
	t.r <- if(startTi == 1) years[1] - by else years[startTi-1]
	tau.denominator <- endT - t.r
	p.r <- pasfr.obs[,ncol(pasfr.obs)]/100. # last observed pasfr
	if(any(is.na(p.r))) stop("Observed PASFR necessary to estimate future trends.")
	p.r <- pmax(p.r, min.value)
	p.r <- p.r/sum(p.r)
	logit.pr <- logit(p.r)
	logit.dif <- logit(gnorm/100.) - logit.pr
	#stop("")
	for(t in 1:ncol(asfr1)){
		asfr1[,t] <- logit.pr + min((years[t+tobs] - t.r)/tau.denominator, 1)*logit.dif
	}
	asfr1 <- inv.logit(asfr1)
	asfr1 <-  scale(asfr1, center=FALSE, scale=colSums(asfr1))
	
	# Continuing of observed trend
	if(startTi < nr.est.points+1) { # the years vector does not include all the observed data
	    yd <- years[1] - by * (nr.est.points - startTi)
	} else yd <- years[startTi - nr.est.points]
	p.e <- pasfr.obs[,ncol(pasfr.obs)-nr.est.points+1]/100.
	if(any(is.na(p.e))) stop("Not enough data on PASFR available to estimate future trends.")
	p.e <- pmax(p.e, min.value)
	p.e <- p.e/sum(p.e)

	tau.denominator2 <- t.r - yd
	logit.dif <- logit.pr - logit(p.e)
	for(t in 1:ncol(asfr2)){
		asfr2[,t] <- logit.pr + ((years[t+tobs] - t.r)/tau.denominator2) *logit.dif
	}
	asfr2 <- inv.logit(asfr2)
	asfr2 <-  scale(asfr2, center=FALSE, scale=colSums(asfr2))
	
	# combining the two trends
	logit.asfr1 <- logit(asfr1)
	logit.asfr2 <- logit(asfr2)
	for(t in 1:ncol(res.asfr)){
		tau <- min((years[t+tobs] - t.r)/tau.denominator, 1)
		res.asfr[,t] <- tau*logit.asfr1[,t] + (1-tau)*logit.asfr2[,t]
	}
	res.asfr <- inv.logit(res.asfr)
	res.asfr <- scale(res.asfr, center=FALSE, scale=colSums(res.asfr))
	#stop("")
	if(start.phase3 <= lyears) res.asfr <- update.by.mac(res.asfr, max(1, start.phase3-tobs))
	return(res.asfr)
}

.get.par.from.inputs <- function(par, inputs, country, convert.to.matrix = TRUE) {
	if(is.null(inputs[[par]])) return (NULL)
	idx <- inputs[[par]][,'country_code'] == country
	if(sum(idx)==0) return (NULL)
	res <- inputs[[par]][idx,,drop=FALSE]
	res2 <- res[, !is.element(colnames(res), c('country_code', 'age')),drop=FALSE]
	if(convert.to.matrix) {
	    res2 <- as.matrix(res2)
	    if('age' %in% colnames(res)) rownames(res2) <- res[, 'age']
	}
    return (res2)
}

.fill.pasfr.ages <- function(dat, fages, check.length.only = FALSE){
    if(is.null(dat)) return(dat)
    if(check.length.only){
        if(nrow(dat) != length(fages))
            stop('PASFR dataset contains ages that are not allowed.\nAllowed: ', paste(fages, collapse = ", "),
                 '\nUsed: ', paste(rownames(dat), collapse = ", "))
        return(dat)
    }
    # fill-in ages if not complete
    if(all(fages %in% rownames(dat))) return(dat)
    datf <- matrix(0, ncol = ncol(dat), nrow = length(fages), dimnames = list(fages, colnames(dat)))
    datf[rownames(dat), ] <- dat
    return(datf)
}

get.country.inputs <- function(country, inputs, nr.traj, country.name) {
	inpc <- list()
	obs <- list()
	for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MXpattern', 'SRB',
				'PASFR', 'PASFRpattern', 'MIGtype', 'MIGm', 'MIGf', 'MXm.pred', 'MXf.pred', 
				'GQm', 'GQf')) {
		inpc[[par]] <- .get.par.from.inputs(par, inputs, country)
		obs[[par]] <- .get.par.from.inputs(par, inputs$observed, country)
	}
	repr.ages <- fert.ages(inputs$annual)
	inpc[['PASFR']] <- .fill.pasfr.ages(inpc[['PASFR']], repr.ages, check.length.only = !inputs$annual)
	obs[['PASFR']] <- .fill.pasfr.ages(obs[['PASFR']], repr.ages, check.length.only = !inputs$annual)
    if(inputs$fixed.pasfr && is.null(inpc[['PASFR']])) {
        warning('No PASFR projection for ', country.name, '. No population projection generated.')
	    return(NULL)
    }
	if(!inputs$fixed.pasfr && is.null(obs[['PASFR']])) {
	    warning('No observed PASFR for ', country.name, '. No population projection generated.')
	    return(NULL)
	}
	if(is.null(inpc[['MXm']]) || is.null(inpc[['MXf']])) {
	    warning('No mortality data for ', country.name, '. No population projection generated.')
	    return(NULL)
	}

	for(par in c('MXpattern', 'PASFRpattern', 'HIVparams')) { # keep these datasets in data.frame format
	    inpc[[par]] <- .get.par.from.inputs(par, inputs, country, convert.to.matrix = FALSE)
	}
	inpc[['MIGBaseYear']] <- inpc[['MIGtype']][,'ProjFirstYear']
	inpc[['MIGtype']] <- inpc[['MIGtype']][,'MigCode']
	# generate sex and age-specific migration if needed
	if((!is.null(inpc[['MIGm']]) && all(is.na(inpc[['MIGm']]))) || (!is.null(inpc[['MIGf']]) && all(is.na(inpc[['MIGf']])))) {
		mig.recon <- age.specific.migration(wpp.year=inputs$wpp.year, countries=country, verbose=FALSE)
		mig.pair <- list(MIGm="male", MIGf="female")
		for(what.mig in names(mig.pair)) {
			if(!is.null(inpc[[what.mig]]) && all(is.na(inpc[[what.mig]]))) {
				# extact predicted migration
				cols <- intersect(colnames(mig.recon[[mig.pair[[what.mig]]]]), colnames(inpc[[what.mig]]))
				inpc[[what.mig]][,cols] <- as.matrix(mig.recon[[mig.pair[[what.mig]]]][,cols])
				rownames(inpc[[what.mig]]) <- rownames(mig.recon[[mig.pair[[what.mig]]]])
				# extact observed migration
				if(!is.null(obs[[what.mig]])) {
					cols <- intersect(colnames(mig.recon[[mig.pair[[what.mig]]]]), colnames(obs[[what.mig]]))
					obs[[what.mig]][,cols] <- as.matrix(mig.recon[[mig.pair[[what.mig]]]][,cols])
					rownames(obs[[what.mig]]) <- rownames(mig.recon[[mig.pair[[what.mig]]]])
				}
			}
		}
	}

	what.traj <- list(TFRpred='TFR', e0Mpred='male e0', e0Fpred='female e0')
	medians <- list()
	lyears <- length(inputs$proj.years)
	for(par in names(what.traj)) {
		if (!is.data.frame(inputs[[par]])) next
		cidx <- inputs[[par]][,'country_code'] == country 
		idx <- cidx & is.element(inputs[[par]][,'year'], inputs$proj.years)
		if(sum(idx) == 0) {
			warning('No ', what.traj[[par]], ' trajectories for ', country.name, 
				'. No population projection generated.')
			return(NULL)
		}
		df <- inputs[[par]][idx,-1]
		utrajs <- sort(unique(df$trajectory))
		ntrajs <- length(utrajs)
		if(ntrajs > 1){
			sorted.df <- data.frame(year=rep(inputs$proj.years, times=ntrajs), trajectory=rep(utrajs, each=lyears))
			# this is to get rows of the data frame in a particular order (align years and trajectories)
			df <- merge(sorted.df, df, sort=FALSE)
		}
		inpc[[par]] <- matrix(df[, 'value'], nrow=lyears)
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
	e <- new.env()
	e$inputs <- inputs
	for(sex in c('M', 'F')) {
		par <- paste0('mig', sex, 'pred')
		if(is.null(inputs[[par]])) next
		inpc[[par]] <- .get.migration.traj(e, par, country)
		if(is.null(inpc[[par]])) next
		medians[[par]] <- apply(inpc[[par]], c(1,3), quantile, 0.5, na.rm = TRUE)
	}
	inpc$migMmedian <- medians$migMpred
	inpc$migFmedian <- medians$migFpred
	
	if(is.null(inpc$TFRpred)) {
		inpc$TFRpred <- get.tfr.trajectories(inputs$TFRpred, country)
		if(is.null(inpc$TFRpred)) {
			warning('No TFR trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)	
		}
		country.obj <- get.country.object(country, inputs$TFRpred$mcmc.set$meta)
		medians$TFRpred <- bayesTFR::get.median.from.prediction(inputs$TFRpred, country.obj$index, country.obj$code)[-1]
		obs.tfr <- bayesTFR:::get.tfr.reconstructed(inputs$TFRpred$tfr_matrix_reconstructed, inputs$TFRpred$mcmc.set$meta)
		obs$TFRpred <- obs.tfr[1:if(!is.null(inputs$TFRpred$present.year.index)) inputs$TFRpred$present.year.index else nrow(obs.tfr),country.obj$index]
	} 
	inpc$TFRhalfchild <- bayesTFR:::get.half.child.variant(median=medians$TFRpred, increment=c(0.25, 0.4, 0.5))
	if(!all(is.element(inputs$proj.years, colnames(inpc$TFRhalfchild))))
		stop('Mismatch in projection periods of TFR and the target projection years.')
	inpc$TFRhalfchild <- inpc$TFRhalfchild[,as.character(inputs$proj.years)]
	
	if(!inputs$fixed.mx){
	if(is.null(inpc$e0Mpred)) {
		inpc$e0Mpred <- get.e0.trajectories(inputs$e0Mpred, country)
		if(is.null(inpc$e0Mpred)) {
			warning('No male e0 trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)	
		}
		country.obj <- get.country.object(country, inputs$e0Mpred$mcmc.set$meta)
		medians$e0Mpred <- bayesTFR::get.median.from.prediction(inputs$e0Mpred, country.obj$index, country.obj$code)
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
		medians$e0Fpred <- bayesTFR::get.median.from.prediction(inputs$e0Fpred, country.obj$index, country.obj$code)
		obs.e0F <- bayesLife:::get.e0.reconstructed(inputs$e0Fpred$e0.matrix.reconstructed, inputs$e0Fpred$mcmc.set$meta)
		obs$e0Fpred <- obs.e0F[1:if(!is.null(inputs$e0Fpred$present.year.index)) inputs$e0Fpred$present.year.index else nrow(obs.e0F),country.obj$index]

	}
	if(!all(is.element(inputs$proj.years, names(medians$e0Fpred))))
		stop('Mismatch in projection periods of female e0 and the target projection years.')
	inpc$e0Fmedian <- medians$e0Fpred[as.character(inputs$proj.years)]
	}
	# select trajectories
	indices <- list()
	nr.traj <- min(nr.traj, max(ncol(inpc$e0Mpred), ncol(inpc$e0Fpred), ncol(inpc$TFRpred), 
							    ncol(inpc$migMpred), ncol(inpc$migFpred)))
	for (par in c('TFRpred', 'e0Fpred', 'migMpred')) {
		if(is.null(inpc[[par]])) next
		traj.available <- ncol(inpc[[par]])
		if(traj.available == nr.traj) {
			indices[[par]] = 1:nr.traj
			next
		}
		if(traj.available > nr.traj) # select equidistantly
			indices[[par]] <- get.traj.index(nr.traj, inpc[[par]]) 
		else # re-sample
			indices[[par]] <- sample(1:traj.available, nr.traj, replace=TRUE)
	}
	# dependent trajectories
	dependencies <- list(e0Fpred='e0Mpred', migMpred='migFpred')
	for(par in names(dependencies)) {
		if(is.null(inpc[[dependencies[[par]]]])) {
			if(par == 'migMpred') {
				if((is.null(inpc[[par]]) || ncol(inpc[[par]])==1)) next
				stop('Female migration trajectories must match male migration.')
			}
			next
		}		
		if(par == 'migMpred' && (is.null(inpc[[par]]) || (!is.null(inpc[[par]]) && ncol(inpc[[par]])==1))) {
			indices[[dependencies[[par]]]] <- rep(1, nr.traj)
			next # using one trajectory, male is default
		}
		traj.available <- ncol(inpc[[par]])
		if (traj.available == ncol(inpc[[dependencies[[par]]]])) # same number of trajectories, therefore same indices
			indices[[dependencies[[par]]]] <- indices[[par]]
		else
			stop('Trajectories for female ', list(e0Fpred='life expectancy', migFpred='migration')[[par]], ' do not match the male ones: ',
				traj.available, ' <> ', ncol(inpc[[dependencies[[par]]]]), ' for ', country.name, '. No population projection generated.')
	}
	for (par in c('TFRpred', 'e0Mpred', 'e0Fpred')) { # these are matrices
		if(is.null(inpc[[par]])) next
		inpc[[par]] <- inpc[[par]][,indices[[par]], drop=FALSE]
		if(tolower(substr(par, 1,3)) %in% tolower(inputs$average.annual) && !inputs$annual) { # average annual data to 5-years data
		    years <- as.integer(rownames(inpc[[par]]))
		    year.ranges <- range(years[years %% 5 == 0])
		    mid.points <- c(0, seq(year.ranges[1]-2, year.ranges[2]+3, by = 5))
		    brks <- seq(year.ranges[1]-5, year.ranges[2] + 5, by = 5)
		    year.bin <- findInterval(years, brks, left.open = TRUE)
		    vals <- apply(inpc[[par]], 2, function(v) aggregate(v, by = list(year.bin), FUN = mean, na.rm = TRUE)[,"x"])
		    rownames(vals) <- mid.points[-c(1, length(mid.points))]
		    inpc[[par]] <- vals
		}
		inpc[[par]] <- inpc[[par]][as.character(inputs$proj.years),, drop=FALSE]
	}
	inpc$mig.nr.traj <- 1
	for(par in c('migMpred', 'migFpred')) { # age-specific, thus 3-d arrays
		if(is.null(inpc[[par]])) next
		inpc[[par]] <- inpc[[par]][,indices[[par]], , drop=FALSE]
		inpc$mig.nr.traj <- length(indices[[par]])
	}
	for(par in c("GQm", "GQf")) {
	    if(is.null(inpc[[par]])) next
	    # match ages
	    age.labels <- get.age.labels(all.ages(inputs$annual, observed = TRUE))
	    if(!all(rownames(inpc[[par]]) %in% age.labels))
	        stop("Mismatch in age labels for ", par, "\nAllowed labels: ", paste(age.labels, collapse = ", "))
	    gq <- rep(0, length(age.labels))
	    names(gq) <- age.labels
	    gq[rownames(inpc[[par]])] <- inpc[[par]]
	    # expand from 100+ to 130+
	    gq <- c(gq, rep(0, all.age.length(inputs$annual, observed = FALSE) - length(gq)))
	    inpc[[par]] <- gq
	}
	    
	inpc$observed <- obs
	inpc$trajectory.indices <- indices
	return(inpc)
}

get.traj.index <- function(nr.traj, traj, traj.dim=2, sample=FALSE) {
	traj.tot <- dim(traj)[traj.dim]
	if(nr.traj >= traj.tot) 
		return(if(sample) sample(1:traj.tot, nr.traj, replace=TRUE) else 1:traj.tot)
	return(if(sample) sample(1:traj.tot, nr.traj) else as.integer(seq(1, traj.tot, length=nr.traj)))
}

rotateLC <- function(e0, bx, bux, axM, axF, e0u=102, p=0.5) {
	# Rotation from Li, Lee, Gerland, 2013
	lmin <- -12
	lmax <- 0
	machine.max <- log(.Machine$double.xmax)
	machine.min <- log(.Machine$double.xmin)
	npred <- length(e0)
	kranges <- ax <- list()
	ku <- rep(NA, npred)
	for(sex in 1:2)
		kranges[[sex]] <- list(ku=ku, kl=ku)
	ax <- list(axM, axF)
	#if(!is.null(dim(bx)) && length(dim(bx))==2) { # aids country; bx is already a matrix
	#	Bxt <- bx
	#} else {
		wt <- (e0 - 80)/(e0u-80)
		wst <- (0.5*(1+(sin(pi/2*(2*wt-1)))))^p
		Bxt <- matrix(NA, nrow=length(bx), ncol=npred)
		for(t in 1:npred) {
			Bxt[,t] <- switch(cut(e0[t], c(0, 80, 102, 9999), labels=FALSE, right=FALSE),
					bx, 
					(1-wst[t])*bx + wst[t]*bux,
					bux)
		}
	#}
	bx.lim=rbind(apply(Bxt, 2, function(x) min(x[x>0])), apply(Bxt, 2, function(x) max(x[x>0])))
	for(sex in 1:2) {
		kranges[[sex]]$kl <- pmin((lmax - apply(ax[[sex]], 2, min))/bx.lim[2,], machine.max)
		kranges[[sex]]$ku <- pmax((lmin - apply(ax[[sex]], 2, max))/bx.lim[1,], machine.min)
	}
	return(list(bx=Bxt, kranges=kranges))
}

.pattern.value <- function(name, pattern, default = NULL, na.means.missing = FALSE) {
    if(is.null(pattern)) return(default)
    val <- if(name %in% colnames(pattern)) pattern[, name] else default
    if(na.means.missing && is.na(val)) val <- default
    return(val)
}

.hiv.mortality <- function(e0m, e0f, country, region, params = NULL) {
    npred <- length(e0m)
    male.mx <- female.mx <- matrix(NA, nrow = 22, ncol = npred, 
                                   dimnames = list(NULL, names(e0m)))
    prevF <- prevM <- rep(3, npred)
    names(prevF) <- names(prevM) <- names(e0m)
    if(!is.null(params)) {
        prev.cols <- names(e0m)[names(e0m) %in% names(params)]
        prevF[prev.cols] <- as.numeric(params[params$param == "prev" & params$sex == "female", 
                               prev.cols])
        prevM[prev.cols] <- as.numeric(params[params$param == "prev" & params$sex == "male", 
                                   prev.cols])
    }
    for(i in 1:ncol(male.mx)) {
        fct <- "hiv.mortmod"
        male.mx[,i] <- do.call(fct, list(e0m[i], prev = prevM[i], sex = 0, region = region))
        female.mx[,i] <- do.call(fct, list(e0f[i], prev = prevF[i], sex = 1, region = region))
    }
    return(list(male = list(mx = male.mx), female = list(mx = female.mx)))
}

.prepare.for.mortality.projection <- function(pattern, mxKan, hiv.params = NULL, lc.for.all = FALSE, annual = FALSE) {
    if(lc.for.all) {
        meth1 <- "LC"
        meth2 <- ""
    } else {
        meth1 <- .pattern.value("AgeMortProjMethod1", pattern, "LC")
        meth2 <- .pattern.value("AgeMortProjMethod2", pattern, "")
    }
    if(is.na(meth2)) meth2 <- ""
    args <- list()
    if("MLT" %in% c(meth1, meth2)) {
        mlttype <- .pattern.value("AgeMortProjPattern", pattern, NULL)
        if(is.null(mlttype)) {
            warning("Column for MLT type (AgeMortProjPattern) is missing. CD_West used.")
            mlttype <- "CD_West"
        }
        args[["MLT"]] <- list(type = mlttype, nx = if(annual) 1 else 5)
    }
    if("PMD" %in% c(meth1, meth2)) {
        adj.code <- .pattern.value("AgeMortProjAdjSR", pattern, 0)
        args[["PMD"]] <- list(
            mxm0 = mxKan$male$mx.orig[,ncol(mxKan$male$mx)],
            mxf0 = mxKan$female$mx.orig[,ncol(mxKan$female$mx)],
            interp.rho = TRUE, keep.lt = TRUE,
            sexratio.adjust = adj.code == 1,
            adjust.sr.if.needed = adj.code == 3,
            nx = if(annual) 1 else 5
        )
    }
    if("LC" %in% c(meth1, meth2)) {
        args[["LC"]] <- list(lc.pars = mxKan, keep.lt = TRUE, constrain.all.ages = TRUE)
    }
    if("HIVmortmod" %in% c(meth1, meth2)) {
        callstr <- 'requireNamespace("HIV.LifeTables")'
        pkg.available <-  eval(parse(text = callstr))
        if(!pkg.available)
            stop("Method HIVmortmod requires the HIV.LifeTables package to be installed. It is available on GitHub in the PPgP/HIV.LifeTables repository.")
        callstr <- "data('HIVModelLifeTables', package = 'HIV.LifeTables')"  # need to do this because requireNamespace does not attach lazy data
        eval(parse(text = callstr))
        args[["HIVmortmod"]] <- list(region = .pattern.value("HIVregion", pattern, 1),
                                     params = hiv.params
        )
        if(! meth2 %in% c("HIVmortmod", "")) {
            warning("HIVmortmod cannot be combined with other methods.")
        }
        meth1 <- "HIVmortmod"
        meth2 <- ""
    }
    if("LogQuad" %in% c(meth1, meth2)) args[["LQ"]] <- list() # no arguments for logquad
    if(meth2 != "") args$weights <- eval(parse(text = .pattern.value("AgeMortProjMethodWeights", pattern, c(1, 0.5))))
    args$meth1 <- meth1
    args$meth2 <- meth2
    return(args)
}

project.mortality <- function (eopm, eopf, npred, ..., mortcast.args = NULL, annual = FALSE, verbose=FALSE, debug=FALSE) {
    args <- if(is.null(mortcast.args)) .prepare.for.mortality.projection(..., annual = annual) else mortcast.args
    nage <- all.age.length(annual, observed = FALSE)
    min.age.groups <- lt.age.length(annual, observed = FALSE)
    kann.proj.ages <- if(annual) 100:130 else seq(100, 130, by = 5)
    kann.est.ages <- if(annual) 90:99 else seq(80, 95, by = 5)
    if(args$meth2 == "") { # apply a single method 
        res <- switch(args$meth1, 
            LC = do.call("mortcast", c(list(eopm, eopf), args[["LC"]])),
            PMD = do.call("copmd", c(list(eopm, eopf), args[["PMD"]])),
            MLT = do.call("mltj", c(list(eopm, eopf), args[["MLT"]])),
            HIVmortmod = do.call(".hiv.mortality", c(list(eopm, eopf), args[["HIVmortmod"]])),
            LogQuad = do.call("logquadj", c(list(eopm, eopf), args[["LQ"]]))
            )
        res <- MortCast:::.apply.kannisto.if.needed(res, min.age.groups = min.age.groups,
                                                    proj.ages = kann.proj.ages, est.ages = kann.est.ages)
    } else { # combination of two methods
        res <- mortcast.blend(eopm, eopf, meth1 = tolower(args$meth1),
                              meth2 = tolower(args$meth2), 
                              weights = args$weights,
                              min.age.groups = min.age.groups,
                              meth1.args = args[[args$meth1]], meth2.args = args[[args$meth2]],
                              kannisto.args = list(proj.ages = kann.proj.ages, est.ages = kann.est.ages))
    }
    # consolidate results which can be in different formats from the different methods
    if(!"mx" %in% names(res))
        res <- list(mx = list(res$male$mx, res$female$mx), sr = list(res$male$sr, res$female$sr),
                    LLm = list(res$male$Lx, res$female$Lx), lx = list(res$male$lx, res$female$lx))
    res$male$sex <- 1
    res$female$sex <- 2
    if(is.null(res$sr[[1]]) || nrow(res$sr[[1]]) != nage) {# compute survival
        srinput <- list(male = list(sex = 1, mx = res$mx[[1]]),
                        female = list(sex = 2, mx = res$mx[[2]]))
        res <- survival.fromLT(npred, srinput, annual = annual, verbose=verbose, debug=debug)
    }
    return(res)
}


survival.fromLT <- function (npred, mxKan, annual = FALSE, observed = FALSE, 
                             include01 = annual, verbose=FALSE, debug=FALSE) {
    nage <- all.age.length(annual, observed = observed)
    nagemx <- lt.age.length(annual, observed = observed)
    sr <- LLm <- lx <- list(matrix(0, nrow=nage, ncol=npred), matrix(0, nrow=nage, ncol=npred))
    Mx <- list(matrix(0, nrow=nagemx, ncol=npred), matrix(0, nrow=nagemx, ncol=npred))
    sx <- rep(0, nage)
    for (mxYKan in list(mxKan$female, mxKan$male)) { # iterate over male and female
    	Mx[[mxYKan$sex]] <- mxYKan$mx
    	sex <- c('Male', 'Female')[mxYKan$sex]
    	for(time in 1:npred) {
			res <- LifeTableMx(mxYKan$mx[,time], sex=sex, abridged = !annual)
			if(!include01) {
			    # collapse first two age groups
			    LLm[[mxYKan$sex]][,time] <- .collapse.Lx(res)  
			    sr[[mxYKan$sex]][,time] <- .collapse.sx(res)			
			    lx[[mxYKan$sex]][,time] <- .collapse.lx(res)
			} else {
			    LLm[[mxYKan$sex]][,time] <- res$Lx
			    sr[[mxYKan$sex]][,time] <- res$sx	
			    lx[[mxYKan$sex]][,time] <- res$lx
			}
		}
    }
	return(list(sr=sr, LLm=LLm, mx=Mx, lx=lx))    
}

runKannisto <- function(inputs, start.year, lc.for.all = FALSE, ...) {
	# extend mx, get LC ax,bx,k1
	KannistoAxBx.joint(inputs$MXm, inputs$MXf, start.year=start.year, mx.pattern=inputs$MXpattern, 
	                   compute.AxBx = lc.for.all || any(c(.pattern.value("AgeMortProjMethod1", inputs$MXpattern, "LC"),
	                                        .pattern.value("AgeMortProjMethod2", inputs$MXpattern, "", na.means.missing = TRUE)) == "LC"), 
	                   ...)
}

runKannisto.noLC <- function(inputs, ...) {
	# extend mx
	KannistoAxBx.joint(inputs$MXm.pred, inputs$MXf.pred, compute.AxBx=FALSE, ...)
}


KannistoAxBx.joint <- function(male.mx, female.mx, start.year=1950, mx.pattern=NULL, ax.latest.periods=99, npred=19, 
								joint=TRUE, compute.AxBx=TRUE, annual = FALSE)  {
	# Extending mx to age 130 using Kannisto model and mx 80-99, OLS
	rownames(male.mx) <- rownames(female.mx) <- if(!annual) c(0,1, seq(5, by=5, length=nrow(male.mx)-2)) else 0:(nrow(male.mx)-1)
	proj.ages <- if(annual) 100:130 else seq(100, 130, by = 5)
	est.ages <- if(annual) 90:99 else seq(80, 95, by = 5)
	if(joint) {
	    kann <- cokannisto(male.mx, female.mx, est.ages = est.ages, proj.ages = proj.ages)
	    Mxe.m <- kann$male
	    Mxe.f <- kann$female
	} else {
	    Mxe.m <- kannisto(male.mx, est.ages = est.ages, proj.ages = proj.ages)
	    Mxe.f <- kannisto(female.mx, est.ages = est.ages, proj.ages = proj.ages)
	}
	result <- list(male=list(mx=Mxe.m, mx.orig = male.mx, sex = 1), 
	               female=list(mx=Mxe.f, mx.orig = female.mx, sex = 2))
	if(!compute.AxBx) return(result)
	#Get Lee-Cater Ax and Bx
	ne <- ncol(Mxe.m)
	old.age <- all.age.length(annual, observed = TRUE)
	nmx <- lt.age.length(annual, observed = FALSE)
	years <- as.integer(substr(colnames(male.mx),1,4))
	year.step <- if(annual) 1 else 5
	first.year <- years[1]
	has.nas.in.old.ages <- FALSE
	if(any(is.na(male.mx[old.age, ]))) {# remove columns that have NAs for old ages
	    first.year <- years[which(!is.na(male.mx[old.age,]))[1]]
	    start.year <- max(start.year, first.year)
	    has.nas.in.old.ages <- TRUE
	}
	ns <- which(years == start.year)
	if(length(ns)==0) stop('start.year must be between ', first.year, ' and ', years[ne])
    model.bx <- .pattern.value("AgeMortalityType", mx.pattern, "") == "Model life tables"
    avg.ax <- .pattern.value("LatestAgeMortalityPattern", mx.pattern, 1) == 0
    smooth.ax <-  !avg.ax && .pattern.value("SmoothLatestAgeMortalityPattern", mx.pattern, 0) == 1
    is.aids.country <- .pattern.value("WPPAIDS", mx.pattern, 0) == 1
    if(is.aids.country) {
    	avg.ax <- FALSE
    	smooth.ax <- TRUE
    	aids.idx <- if(!has.nas.in.old.ages) which(years < 1985) else 1:length(years)
    	aids.npred <- min((2100-(as.integer(years[ne])+year.step))/year.step, npred)
    }
    if(!avg.ax && !is.null(lpat <- .pattern.value("LatestAgeMortalityPattern", mx.pattern, NULL))) {
        # lpat should not be zero because of the !avg.ax condition, but it can be negative for removing time periods
        ax.latest.periods <- lpat 
    }
    mlt.bx <- NULL
    if(model.bx) {
    	bx.env <- new.env()
    	data(MLTbx, envir = bx.env)
    	bx.pattern <- .pattern.value("AgeMortalityPattern", mx.pattern, "UN_General")
    	if (!bx.pattern %in% rownames(bx.env$MLTbx)) bx.pattern <- "UN_General"
    	mlt.bx <- as.numeric(bx.env$MLTbx[bx.pattern,])
    }
    #this.ns <- if(any(is.na(result$male$mx[,ns:ne]))) 
    #    ns + sum(apply(result$male$mx[,ns:ne], 2, function(z) all(is.na(z))))
    #else ns
    length.mx <- length(ns:ne)
    if(ax.latest.periods < 0) { # remove ax.latest.periods from the end
        ax.index <- (1:length.mx)[-(max(length.mx+ax.latest.periods+1, 2):length.mx)] # at least one period should stay in
    } else { # take the ax.latest.periods latest time periods
        if(avg.ax || ax.latest.periods == 0) ax.index <- 1:length.mx
        else {
            ax.ns <- max(length.mx - ax.latest.periods+1, 1)
            ax.index <- ax.ns:length.mx
        }
    }
    lc.est <- lileecarter.estimate(result$male$mx[,ns:ne], result$female$mx[,ns:ne],
                                   ax.index = ax.index, ax.smooth = smooth.ax, nx = year.step)

    if(is.aids.country) { # modify ax and bx
        for(sex in c('male', 'female')) {
    	    lMxe <- log(result[[sex]]$mx)
		    axt <- matrix(lc.est[[sex]]$ax, nrow=nmx, ncol=npred)
			ax.end <- apply(lMxe[,aids.idx, drop=FALSE], 1, sum, na.rm=TRUE)/length(aids.idx)
			ax.end.sm <- smooth.spline(ax.end[1:old.age], df=11)$y
    		ax.end[2:old.age] <- ax.end.sm[2:old.age] # keep value of the first age group
			for (i in 1:nmx) { # linear interpolation to the average ax ending in 2050; after that the avg ax is used
				axt[i,1:aids.npred] <- approx(c(1,aids.npred), c(lc.est[[sex]]$ax[i], ax.end[i]), xout=1:aids.npred)$y
				if(aids.npred < npred)
					axt[i,(aids.npred+1):npred] <- ax.end[i]	
			}
    		lc.est[[sex]]$axt <- axt
        }
    }
    if(model.bx) {
        names(mlt.bx) <- names(lc.est$male$bx)
        lc.est$male$bx <- lc.est$female$bx <- mlt.bx
        lc.est$bx <- (lc.est$male$bx + lc.est$female$bx)/2
        lc.est$ultimate.bx <- ultimate.bx(lc.est$bx)
    }
    # merge results
    #sex.code <- list(male = 1, female = 2)
    for(sex in c('male', 'female')) {
        lc.est[[sex]] <- c(lc.est[[sex]], result[[sex]])
        #lc.est[[sex]]$sex <- sex.code[[sex]]
    }
	return(lc.est)
}


StoPopProj <- function(npred, inputs, LT, asfr, mig.pred=NULL, mig.type=NULL, country.name=NULL, 
                       keep.vital.events=FALSE, annual = FALSE) {
    change.by.gq <- function(gq, pop, factor = -1){
        pop <- pop + factor * gq
    }
    nagecat <- all.age.length(annual, observed = FALSE)
    nbagecat <- fert.age.length(annual)
	popm <- popf <- matrix(0, nrow=nagecat, ncol=npred+1)
	popm[,1] <- c(inputs$POPm0, rep(0, nagecat - nrow(inputs$POPm0)))
	popf[,1] <- c(inputs$POPf0, rep(0, nagecat - nrow(inputs$POPf0)))
	use.gq <- FALSE
	if(!is.null(inputs$GQm)) {
	    popm[,1] <- change.by.gq(inputs$GQm, popm[,1])
	    use.gq <- TRUE
	}
	if(!is.null(inputs$GQf)) {
	    popf[,1] <- change.by.gq(inputs$GQf, popf[,1])
	    use.gq <- TRUE
	}
	totp <- c(sum(popm[,1]+popf[,1]), rep(0, npred))
	btageM <- btageF <- matrix(0, nrow=nbagecat, ncol=npred) # births by age of mother and sex of child
	deathsM <- deathsF <- matrix(0, nrow=nagecat, ncol=npred)
	nproj <- npred
	migM <- as.matrix(if(!is.null(mig.pred[['M']])) mig.pred[['M']] else inputs[['MIGm']])
	migF <- as.matrix(if(!is.null(mig.pred[['F']])) mig.pred[['F']] else inputs[['MIGf']])
	observed <- 0

	res <- .C("CCM", as.integer(observed), as.integer(!annual), as.integer(nproj), 
	            as.numeric(migM), as.numeric(migF), nrow(migM), ncol(migM), as.integer(mig.type),
		        srm=LT$sr[[1]], srf=LT$sr[[2]], asfr=as.numeric(as.matrix(asfr)), 
		        srb=as.numeric(as.matrix(inputs$SRB)), 
		        Lm=LT$LLm[[1]], Lf=LT$LLm[[2]], lxm=LT$lx[[1]], lxf=LT$lx[[2]],
		        nages = as.integer(nagecat), nfages = as.integer(nbagecat), 
		        fstart = as.integer(fert.age.index(annual)[1]),
		        popm=popm, popf=popf, totp=totp,
		        btagem=as.numeric(btageM), btagef=as.numeric(btageF), 
		        deathsm=as.numeric(deathsM), deathsf=as.numeric(deathsF)
		        )

	if(use.gq) {
	    if(!is.null(inputs$GQm)) res$popm <- change.by.gq(inputs$GQm, res$popm, factor = 1)
	    if(!is.null(inputs$GQf)) res$popf <- change.by.gq(inputs$GQf, res$popf, factor = 1)
	    res$totp <- colSums(res$popm + res$popf)
	}
	
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

compute.observedVE <- function(inputs, pop.matrix, mig.type, mxKan, country.code, estim.years, annual = FALSE) {
	obs <- inputs$observed
	if(is.null(obs$PASFR)) return(NULL)
	npasfr <- nrow(obs$PASFR)
	nest <- max(min(length(obs$TFRpred), ncol(obs$PASFR), sum(!is.na(obs$MIGm[1,])), length(estim.years),
	            ncol(pop.matrix[[1]])-1), 1)
	estim.years <- estim.years[(length(estim.years)-nest+1):length(estim.years)]
	pasfr <- obs$PASFR[,(ncol(obs$PASFR)-nest+1):ncol(obs$PASFR), drop=FALSE]
	tfr <- obs$TFRpred[(length(obs$TFRpred)-nest+1):length(obs$TFRpred)]
	mig.data <- list(as.matrix(obs$MIGm[,(ncol(obs$MIGm)-nest+1):ncol(obs$MIGm)]), 
					as.matrix(obs$MIGf[,(ncol(obs$MIGf)-nest+1):ncol(obs$MIGf)]))
	asfr <- pasfr/100.
	for(i in 1:npasfr) asfr[i,] <- tfr * asfr[i,]
	
	maxage <- all.age.length(annual = annual, observed = TRUE)
	reprod.age <- fert.age.index(annual)
	nfertages <- fert.age.length(annual)
	pop <- D10 <- list()
	nmx <- ncol(inputs$MXm)
	mx <-  list(inputs$MXm[,(nmx-nest+1):nmx, drop=FALSE], inputs$MXf[,(nmx-nest+1):nmx, drop=FALSE])
	srb <- obs$SRB[(length(obs$SRB)-nest+1):length(obs$SRB)]
	srb.ratio <- srb / (1 + srb)
	sr <- deaths <- list(matrix(0, nrow=maxage, ncol=nest), matrix(0, nrow=maxage, ncol=nest))
	births <- list(matrix(0, nrow=nfertages, ncol=nest), matrix(0, nrow=nfertages, ncol=nest))
	for(sex in 1:2) {
	    pop[[sex]] <- matrix(0, nrow=maxage, ncol=nest+1)
	    rownames(pop[[sex]]) <- rownames(mig.data[[sex]])
		popage <- get.pop.observed.with.age(NULL, country.code, sex=c('male', 'female')[sex], 
						data=pop.matrix, annual = annual)
		popage <- popage$data[popage$age.idx,(ncol(popage$data)-nest):ncol(popage$data), drop = FALSE]
		pop[[sex]][1:nrow(popage), 1:ncol(popage)] <- as.matrix(popage)
		pop[[sex]][is.na(pop[[sex]])] <- 0 # set pop for ages with NA to 0
	}
	nobs <- ncol(popage)
    totp <- rep(0, ncol(pop[[1]])) # not used for observed data
	LTinputs <- list(male = list(sex = 1, mx = mx[[1]]), female = list(sex = 2, mx = mx[[2]]))
	LT <- survival.fromLT(nest, LTinputs, annual = annual, observed = TRUE)

	res <- .C("CCM", as.integer(nobs), as.integer(!annual), as.integer(nest), 
	              as.numeric(mig.data[[1]]), as.numeric(mig.data[[2]]), 
	              nrow(mig.data[[1]]), ncol(mig.data[[1]]), as.integer(mig.type),
	              srm=LT$sr[[1]], srf=LT$sr[[2]], asfr=as.numeric(as.matrix(asfr)), 
	              srb=as.numeric(as.matrix(srb)), 
	              Lm=LT$LLm[[1]], Lf=LT$LLm[[2]], lxm=LT$lx[[1]], lxf=LT$lx[[2]],
	              nages = as.integer(maxage), nfages = as.integer(nfertages), fstart = as.integer(reprod.age[1]),
	              popm=as.numeric(pop[[1]]), popf=as.numeric(pop[[2]]), totp=totp,
	              btagem=as.numeric(births[[1]]), btagef=as.numeric(births[[2]]), 
	              deathsm=as.numeric(deaths[[1]]), deathsf=as.numeric(deaths[[2]])
	            ) 
    deaths[[1]] <- matrix(res$deathsm, ncol = nest)
    deaths[[2]] <- matrix(res$deathsf, ncol = nest)
    births[[1]] <- matrix(res$btagem, ncol = nest)
    births[[2]] <- matrix(res$btagef, ncol = nest)
    colnames(deaths[[1]]) <- colnames(deaths[[2]]) <- colnames(births[[1]]) <- colnames(births[[2]]) <- colnames(asfr) <- estim.years
    rownames(deaths[[1]]) <- rownames(deaths[[2]]) <- all.ages(annual = annual, observed = TRUE)
    rownames(births[[1]]) <- rownames(births[[2]]) <- rownames(asfr) <- fert.ages(annual = annual)

	res <- list(btm=births[[1]], btf=births[[2]], 
				deathsm=deaths[[1]], deathsf=deaths[[2]], asfert=asfr, pasfert=pasfr, 
				mxm=mx[[1]], mxf=mx[[2]])
	for(par in names(res))
		res[[par]] <- abind(res[[par]], NULL, along=3) # add trajectory dimension
	return(res)
}

write.pop.projection.summary <- function(pop.pred, what=NULL, expression=NULL, output.dir=NULL, ...) {
	pred <- as.environment(pop.pred) # to avoid lots of copying and allowing attaching cache
	if (is.null(output.dir)) output.dir <- pop.output.directory(pred)
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.what <- c('pop', 'popsex', 'popsexage', 'popage', 'births', 'birthssex', 'birthsage', 'birthssexage', 
			'deaths', 'deathssex', 'deathsage', 'deathssexage', 'srsexage', 'fertility', 'fertilityage', 'pfertilityage')
	what <- if(is.null(what)) all.what else match.arg(what, all.what, several.ok=TRUE)
	params <- list()
	if(!is.null(expression)) {
		what <- 'expression'
		params <- list(expression=expression)
	}
	for(summary.type in what) 
		do.call(paste0('write.', summary.type), c(list(pred, output.dir=output.dir), params, ...))
}

write.pop <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, ...)
	
write.popsex <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, what.log='population', ...)
	
write.popsexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, what.log='population', ...)
	
write.popage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, what.log='population', ...)
	
write.births <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='births', 
			file.suffix='births', what.log='total births', ...)
	
write.birthssex <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, vital.event='births', 
			file.suffix='births', what.log='births', ...)

write.birthsage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='births', 
			file.suffix='births', what.log='births', ...)

write.birthssexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='births', 
			file.suffix='births', what.log='births', ...)

write.deaths <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='deaths', 
			file.suffix='deaths', what.log='total deaths', ...)
	
write.deathssex <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths', ...)

write.deathsage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths', ...)
	
write.deathssexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths', ...)
			
write.srsexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='survival', 
			file.suffix='sr', what.log='survival ratio', digits=litem('digits', list(...), 4))
			
write.fertility <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='fertility', 
			file.suffix='asfr', what.log='fertility rate', digits=litem('digits', list(...), 4))
			
write.fertilityage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='fertility', 
			file.suffix='asfr', what.log='fertility rate', digits=litem('digits', list(...), 4))
			
write.pfertilityage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='pasfr', 
			file.suffix='pasfr', what.log='percent fertility rate', digits=litem('digits', list(...), 4))
	
write.expression <- function(pop.pred, expression, output.dir, file.suffix='expression', 
							expression.label=expression, include.observed=FALSE, digits=NULL, 
							adjust=FALSE, adj.to.file=NULL, end.time.only=FALSE) {
	cat('Creating summary file for expression ', expression, ' ...\n')
	header <- list(country.name='country_name',  country.code='country_code', variant='variant')
	variant.names <- c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95')
	nr.var <- length(variant.names)
	pred.period <- get.pop.prediction.periods(pop.pred, end.time.only=end.time.only)
	nr.proj <- length(pred.period)
	nr.obs <- 0
	if(include.observed) {
		obs.period <- get.pop.observed.periods(pop.pred, end.time.only=end.time.only)
		nr.obs <- length(obs.period)-1
		obs.period <- obs.period[-(nr.obs+1)] # remove the last one because the same as the first projection period
		for(iyear in 1:nr.obs) header[[paste0('year', iyear)]] <- obs.period[iyear]
	}
	for(iyear in 1:nr.proj) header[[paste0('year', iyear+nr.obs)]] <- pred.period[iyear]
	col.names <- grep('year', names(header), value=TRUE)
	result <- NULL
	if(include.observed) {
		#res <- get.pop.observed.from.expression.all.countries(expression, pop.pred, time.index=1:nr.obs)
	    res <- c()
	    for(iyear in 1:nr.obs)
		    res <- cbind(res, get.pop.from.expression.all.countries(expression, pop.pred, 
		                                                            time.index=iyear, observed = TRUE))
		#copy the same data into the variant rows 
		result <- matrix(NA, nrow=nrow(res)*5, ncol=ncol(res))
		for(i in 1:5) result[seq(i,by=5, length=nrow(res)),] <- res
	}
	if(adjust && is.null(pop.pred$adjust.env)) pop.pred$adjust.env <- new.env()
	for(iyear in 1:nr.proj) {	
		result <- cbind(result, as.vector(t(get.pop.from.expression.all.countries(expression, pop.pred, 
						quantiles=c(0.5, 0.1, 0.9, 0.025, 0.975), time.index=iyear, adjust=adjust, adj.to.file=adj.to.file))))
	}
	if(!is.null(digits)) result <- round(result, digits)
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
							what.log='total population', include.observed=FALSE, digits=0, adjust=FALSE, 
							end.time.only=FALSE) {
	cat('Creating summary file of ', what.log, ' ')
	if(bysex) cat('by sex ')
	if(byage) cat('by age ')
	cat ('...\n')
	header <- list(country.name='country_name',  country.code='country_code', variant='variant')
	if(bysex) header[['sex']] <- 'sex'
	if(byage) header[['age']] <- 'age'
	variant.names <- c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95')
	nr.var <- length(variant.names)
	if(missing(end.time.only)) end.time.only <- is.null(vital.event)
	pred.period <- get.pop.prediction.periods(pop.pred, end.time.only=end.time.only)
	#if(!is.null(vital.event)) pred.period <- pred.period[2:length(pred.period)]
	nr.obs <- 0
	if(include.observed) {
		obs.period <- get.pop.observed.periods(pop.pred, end.time.only=end.time.only)
		nr.obs <- length(obs.period)-1
		obs.period <- obs.period[-(nr.obs+1)] # remove the last one because the same as the first projection period
		for(iyear in 1:nr.obs) header[[paste0('year', iyear)]] <- obs.period[iyear]
	}
	nr.proj <- length(pred.period)
	for (i in 1:nr.proj) 
		header[[paste0('year', i+nr.obs)]] <- pred.period[i]
	col.names <- grep('year', names(header), value=TRUE)
	result <- NULL
	sex.index <- c(TRUE, FALSE, FALSE)
	if(bysex) sex.index <- !sex.index
	age.index <- c(TRUE, rep(FALSE, length(pop.pred$ages)))
	if(byage) age.index <- !age.index
	ages <- 1:length(pop.pred$ages)
	if(adjust && is.null(pop.pred$adjust.env)) pop.pred$adjust.env <- new.env()
	observed.data <- NULL
	for (country in 1:nrow(pop.pred$countries)) {
		country.obj <- get.country.object(country, country.table=pop.pred$countries, index=TRUE)
		for(sex in c('both', 'male', 'female')[sex.index]) {
			if(!is.null(vital.event)) {
			 	sum.over.ages <- age.index[1]
			 	if(include.observed) 
			 		observed <- get.popVE.trajectories.and.quantiles(pop.pred, country.obj$code, event=vital.event, 
										sex=sex, age='all', sum.over.ages=sum.over.ages, is.observed=TRUE)
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
					this.result <- cbind(this.result, age=rep(get.age.labels(pop.pred$ages, single.year = pop.pred$annual)[age], nr.var))
				}
				if(is.null(vital.event)) {
					if(include.observed) 
						observed.data <- get.pop.observed(pop.pred, country.obj$code, sex=sex, age=age)
					quant <- get.pop.trajectories(pop.pred, country.obj$code, nr.traj=0, sex=sex, age=age, adjust=adjust)$quantiles
					traj <- NULL
					reload <- TRUE
				} else { # vital event
					quant <- traj.and.quantiles$quantiles
					traj <- traj.and.quantiles$trajectories
					if(include.observed)
						observed.data <- observed$trajectories[,,1]
					if(!sum.over.ages) {
						quant <- quant[age-subtract.from.age,,]
						traj <- traj[age-subtract.from.age,,]
						if(include.observed) {
							if (age-subtract.from.age > nrow(observed.data)) # because observed goes only up to 100+
								observed.data <- rep(0, ncol(observed.data))
							else
								observed.data <- observed.data[age-subtract.from.age,]
						}
					}
					reload <- FALSE
					#stop('')
				}
				proj.result <- round(rbind(
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, q=0.5, 
											trajectories=traj, reload=reload, sex=sex, age=age), 
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, pi=80, 
											trajectories=traj, reload=reload, sex=sex, age=age),
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, pi=95, 
											trajectories=traj, reload=reload, sex=sex, age=age)),
					digits)
				if(!is.null(observed.data)) {
					# put it into the same shape as proj.result minus the last observed
					observed.data <- round(rbind(observed.data, NULL), digits)
					observed.data <- observed.data[rep(1, nrow(proj.result)), -ncol(observed.data)]
					proj.result <- cbind(observed.data, proj.result)
				}
				colnames(proj.result) <- col.names
				#stop('')
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

LifeTableMxCol <- function(mx, colname=c('Lx', 'lx', 'qx', 'mx', 'dx', 'Tx', 'sx', 'ex', 'ax'), ...){
	colname <- match.arg(colname)
	if(is.null(dim(mx))) return(.doLifeTableMxCol(mx, colname, ...))
	return(apply(mx, 2, .doLifeTableMxCol, colname=colname, ...))
}

.collapse.sx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	# sx does not need to be collapsed, as all elements refer to a five-year age group
	# for collapsed format remove last age group so to assure same length as all other columns after collapsing	
	#sx.start <- c(LT$sx[1:2], (LT$Lx[1] + LT$Lx[2])/5)[age05]
	#return(c(sx.start, LT$sx[-(1:2)]))
	return(switch(sum(age05)+1,
				LT$sx[c(-1, -length(LT$sx))], # unuseful case of all values FALSE
				LT$sx[-length(LT$sx)], # one value TRUE (collapsed life table)
				LT$sx, # two values TRUE (abridged life table)
				c(LT$sx[1], LT$sx) # unuseful case of all values TRUE
				))
}

.collapse.dx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	dx.start <- c(LT$dx[1:2], LT$dx[1] + LT$dx[2])[age05]
	return(c(dx.start, LT$dx[-(1:2)]))
}

.collapse.ax <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	ax.start <- c(LT$ax[1:2], ((LT$Lx[1]+LT$Lx[2]-5*LT$lx[3])/(LT$dx[1]+LT$dx[2])))[age05]
	return(c(ax.start, LT$ax[-(1:2)]))
}

.collapse.Tx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	Tx.start <- c(LT$Tx[c(1,2,1)])[age05]
	return(c(Tx.start, LT$Tx[-(1:2)]))
}

.collapse.ex <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	ex.start <- c(LT$ex[c(1,2,1)])[age05]
	return(c(ex.start, LT$ex[-(1:2)]))
}

.collapse.Lx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	Lx.start <- c(LT$Lx[1:2], LT$Lx[1] + LT$Lx[2])[age05]
	return(c(Lx.start, LT$Lx[-(1:2)]))
}

.collapse.lx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	lx.start <- LT$lx[c(1,2,1)][age05]
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

.doLifeTableMxCol <- function(mx, colname, age05=c(FALSE, FALSE, TRUE), abridged = TRUE, ...) {
	# age05 determines the inclusion of ages 0-1, 1-4, 0-4
	LT <- LifeTableMx(mx, abridged = abridged, ...)
	if(!abridged || all(age05==c(TRUE, TRUE, FALSE))) { # no collapsing
		result <- LT[,colname]
		names(result) <- rownames(LT)
	} else {
		result <- do.call(paste('.collapse', colname, sep='.'), list(LT, age05=age05))
		names(result) <- c(c('0-1', '1-4', '0-4')[age05], rownames(LT)[-(1:2)])
	}
	return(result)
}


LifeTableMx <- function(mx, sex=c('Male', 'Female', 'Total'), include01=TRUE, 
                        abridged = TRUE, radix = 1, open.age = 130){
	# The first two elements of mx must correspond to 0-1 and 1-4. 
	# If include01 is FALSE, the first two age groups of the results are collapsed to 0-5
    sex <- tolower(match.arg(sex))
    LT <- MortCast::life.table(mx, sex = sex, abridged = abridged, radix = radix, open.age = open.age)
    if(!include01 && abridged) {
        if(all(is.na(LT$ax))) return(LT[-2,])
        age05 <- c(FALSE, FALSE, TRUE)
        LTres <- data.frame(age=LT$age[-2])
        for(colname in setdiff(colnames(LT), "age"))
            LTres[[colname]] <- do.call(paste('.collapse', colname, sep='.'), list(LT, age05=age05))
        rownames(LTres) <- rownames(LT[-2,])
        rownames(LTres)[1] <- "0-4"
        LT <- LTres
    }
	return(LT)
}

unblock.gtk.if.needed <- function(status, gui.opts=list()) {
	# This is to unblock the GUI, if the run is invoked from bayesDem
	# and pass info about its status
	# In such a case the gtk libraries are already loaded
	if(getOption('bDem.Poppred', default=FALSE)) {
		gui.opts$bDem.Poppred.status <- status
		unblock.gtk('bDem.Poppred', gui.opts)
	}
}

unblock.gtk <- function(...) bayesTFR:::unblock.gtk(...)

create.pop.cluster <- function(nr.nodes, ...) {
	cl <- makeCluster(nr.nodes, ...)
	lib.paths <- .libPaths()
	clusterExport(cl, c("lib.paths"), envir=environment())
	clusterEvalQ(cl, {.libPaths(lib.paths); library(bayesPop)})
	return(cl)
}

age.specific.migration <- function(wpp.year=2019, years=seq(1955, 2100, by=5), countries=NULL, smooth=TRUE, 
									rescale=TRUE, ages.to.zero=18:21,
									write.to.disk=FALSE, directory=getwd(), file.prefix="migration", 
									depratio=wpp.year == 2015, verbose=TRUE) {
	# Reconstruct sex- and age-specific net migration using a residual method using wpp data on population
	# and other available indicators. It is scaled to the total net migration for each country. 
	# It is not balanced over the world. Due to rounding issues, often it results in zig-zags over ages,
	# therefore it is smoothed (in a double pass through the smoother).
	# If raw residuals are desired, set smooth=FALSE, rescale=FALSE and ages.to.zero=c().
	if(verbose) {
		status.text <- paste('Reconstructing sex- and age-specific migration from', paste0('wpp', wpp.year, ' '))
		cat('\n', status.text)
	}
	popm0 <- load.wpp.dataset("popM", wpp.year)
	popm0.num.cols <- grep('^[0-9]{4}$', colnames(popm0), value=TRUE) # values of year-columns
	popf0 <- load.wpp.dataset("popF", wpp.year)
	popf0.num.cols <- grep('^[0-9]{4}$', colnames(popf0), value=TRUE)
	popmproj <- load.wpp.dataset("popMprojMed", wpp.year)
	popmproj.num.cols <- grep('^[0-9]{4}$', colnames(popmproj), value=TRUE)
	popfproj <- load.wpp.dataset("popFprojMed", wpp.year)
	popfproj.num.cols <- grep('^[0-9]{4}$', colnames(popfproj), value=TRUE)
	sexrat <- load.wpp.dataset("sexRatio", wpp.year)
	sexrat.num.cols <- grep('^[0-9]{4}', colnames(sexrat), value=TRUE)
	mxm <- load.wpp.dataset("mxM", wpp.year)
	mxm.num.cols <- grep('^[0-9]{4}', colnames(mxm), value=TRUE)
	mxf <- load.wpp.dataset("mxF", wpp.year)
	mxf.num.cols <- grep('^[0-9]{4}', colnames(mxf), value=TRUE)
	mig <- load.wpp.dataset("migration", wpp.year)
	mig.num.cols <- grep('^[0-9]{4}', colnames(mig), value=TRUE)
	tfrproj <- .load.wpp.traj('tfr', wpp.year, median.only=TRUE)
	pasfr <- load.wpp.dataset("percentASFR", wpp.year)
	pasfr.num.cols <- grep('^[0-9]{4}', colnames(pasfr), value=TRUE)
	vwBase <- get(paste0('vwBaseYear', wpp.year))[,c('country_code', 'MigCode')]
	pop.first.country <- popm0[popm0$country_code == mig$country_code[1],]
	max.ages <- nrow(pop.first.country)
	ages <- 1:max.ages
	age.labels <- get.age.labels(ages, age.is.index=TRUE, last.open=TRUE)
	years.periods <- paste(years-5, years, sep="-")
	lyears <- length(years)
	if(is.null(countries)) {
		countries <- mig$country_code
		# filter out non-countries
		if(!exists("UNlocations")) 
            bayesTFR:::load.bdem.dataset('UNlocations', wpp.year, envir=globalenv())
		locs <- UNcountries()
		#locs <- bayesTFR:::load.bdem.dataset('UNlocations', wpp.year, envir=globalenv())
		countries <- countries[countries %in% locs]
	} else mig <- mig[which(mig$country_code %in% countries),]
	depratio.correction <- FALSE
	if (depratio == TRUE || is.character(depratio)) {
		# if it's character it is a name of an rda file; if it's TRUE, take the default file.
		# must have objects depratioM and depratioF
		# which are data frames with columns country_code, period and three dependency ratio 
		# columns (for age groups 0-4, 5-9, 10-14).
		# They represent ratios of that age group to age group 20-25.
		edr <- new.env()
		if(is.character(depratio))
			load(depratio, envir=edr)
		else do.call("data", list(paste0("migdepratio_", wpp.year), envir=edr))
		if(!exists("depratioM", envir=edr) || !exists("depratioF", envir=edr))
			stop("The depratio object must contain objects called depratioM and depratioF\nContains: ", paste(ls(edr), collapse=", "))
		if(ncol(edr$depratioM) < 5 || ncol(edr$depratioF) < 5 || !all(c('country_code', "period") %in% colnames(edr$depratioM))
				|| !all(c('country_code', "period") %in% colnames(edr$depratioF)))
			stop("Objects depratioM and depratioF must contain at least 5 columns (country_code, period and three dependency ratio columns).")
		ratio.colsM <- (1:ncol(edr$depratioM))[-which(colnames(edr$depratioM) %in% c('country_code', "period"))]
		ratio.colsF <- (1:ncol(edr$depratioF))[-which(colnames(edr$depratioM) %in% c('country_code', "period"))]
		depratio.correction <- TRUE
	}
	all.migM <- all.migF <- NULL
	lcountries <- length(countries)
	for(icountry in 1:lcountries) {
		if(verbose && interactive()) cat('\r', status.text, round(icountry/lcountries*100), '%')
		country <- countries[icountry]
		country.name <- as.character(mig[mig$country_code==country, 'name'])
		# filter country data
		popm.obs <- popm0[popm0$country_code==country, popm0.num.cols]
		popf.obs <- popf0[popf0$country_code==country, popf0.num.cols]
		pop1m <- cbind(popm.obs, popmproj[popmproj$country_code==country, popmproj.num.cols])
		pop1f <- cbind(popf.obs, popfproj[popfproj$country_code==country, popfproj.num.cols])
		tfra <- tfrproj[tfrproj$country_code==country,]
		asfr <- pasfr[pasfr$country_code==country, pasfr.num.cols]
		sr <- sexrat[sexrat$country_code==country, sexrat.num.cols] 
		mortM <- mxm[mxm$country_code==country, mxm.num.cols]
		mortF <- mxf[mxf$country_code==country, mxf.num.cols]
		totmig <- mig[mig$country_code==country, mig.num.cols]
		mtype <- vwBase[vwBase$country_code==country,'MigCode']
		if(length(mtype)==0) mtype <- 9
		this.all.migM <- this.all.migF <- data.frame(
				country_code=rep(country, max.ages), name=rep(country.name, max.ages), age=age.labels)
		for(iyear in 1:lyears) {
			year <- years[iyear]
			year.col <- years.periods[iyear]
			year.char <- as.character(year)
			pop0m <- pop1m[,as.character(year-5)]
			pop0f <- pop1f[,as.character(year-5)]
			mortMy <- mortM[,year.col]
			mortFy <- mortF[,year.col]
			sxm <- get.survival(matrix(mortMy, ncol=1), sex="Male")[,1,1]
      		sxf <- get.survival(matrix(mortFy, ncol=1), sex="Female")[,1,1]
			totmigy <- round(totmig[,year.col],3)
			if(totmigy == 0) netmigM <- netmigF <- rep(0, max.ages)
			else {			
				B2 <- sum((pop1f[4:10,year.char] + pop0f[4:10])/2 * tfra[tfra$year==year-2,'value'] * asfr[,year.col]/100)
				netmigM <- c(NA, pop1m[2:max.ages,year.char] - (pop0m[1:(max.ages-1)] * sxm[2:max.ages]))
				netmigF <- c(NA, pop1f[2:max.ages,year.char] - (pop0f[1:(max.ages-1)] * sxf[2:max.ages]))
				B2m <- B2 * sr[,year.col]/(1+sr[,year.col])
				netmigM[1] <- pop1m[1,year.char] - B2m * sxm[1]
				netmigF[1] <- pop1f[1,year.char] - (B2 - B2m) * sxf[1]
				migdata <- list(M=netmigM, F=netmigF)
				sxdata <- list(M=sxm, F=sxf)
				for(sex in c('M', 'F')) {
				    # In wpp2017, for some past years population is reported only up to 85+. 
				    # Set migration of the open age group to 0. 
				    if(any(is.na(pop1m[1:max.ages]))) {
				        # find the index of the first NA 
				        ina <- which(is.na(migdata[[sex]])==TRUE)[1]
				        migdata[[sex]][ina-1] = 0
				    }
					if(mtype == 0) { 
					  # Migration distributed across the time interval.
					  # In projections in this case, the migration is derived as 
					  # M'_a = (M_a + M_{a-1}*sx_a)/2, M'_0 = M_0/2
					  # Thus, here is the reverse of that. 
					  # However, it can yield zig-zags, which are removed in the smoothing step.
					  migdata[[sex]][1] <- 2*migdata[[sex]][1]
						for(i in 2:max.ages) {
					       migdata[[sex]][i] <- 2*migdata[[sex]][i] - migdata[[sex]][i-1]*sxdata[[sex]][i]
						}
					  #stop('')
					}
					migdata[[sex]][ages.to.zero] <- 0
					if(smooth) { #smoothing					
						for(izig in 1:2) { # two passes of smoothing
							# are there significant zig-zags?						
							tops <- ((migdata[[sex]][2:(max.ages-1)] > migdata[[sex]][1:(max.ages-2)] & 
											migdata[[sex]][2:(max.ages-1)] > migdata[[sex]][3:max.ages]) | 
								(migdata[[sex]][2:(max.ages-1)] < migdata[[sex]][1:(max.ages-2)] & migdata[[sex]][2:(max.ages-1)] < migdata[[sex]][3:max.ages])) & 
								abs(diff(migdata[[sex]])[1:(max.ages-2)]/(totmigy/100)) > 0.05
							cs <- cumsum(c(TRUE, tops)) # consider first point as top
							if(any(cs[3:(max.ages-1)] > 2 & cs[3:(max.ages-1)] - cs[1:(max.ages-3)] > 1))  { # at least 3 neighboring tops
								migdata[[sex]] <- smooth.spline(migdata[[sex]], df=10)$y # smooth
								migdata[[sex]][ages.to.zero] <- 0
							} else break
						}
					}
				}
				netmigM <- migdata[['M']]
				netmigF <- migdata[['F']]
				if(depratio.correction) {
					# correct dependency ratio
					cntry <- country
					rowM <- edr$depratioM[edr$depratioM$country_code==cntry & edr$depratioM$period==year.col, ratio.colsM]
					if(nrow(rowM) > 0 && !any(is.na(rowM))) netmigM[1:3] <- as.double(netmigM[5]*rowM)
					rowF <- edr$depratioF[edr$depratioF$country_code==cntry & edr$depratioF$period==year.col, ratio.colsF]
					if(nrow(rowF) > 0 && !any(is.na(rowF))) netmigF[1:3] <- as.double(netmigF[5]*rowF)
				}
				if(rescale) {
					s <- sum(netmigM + netmigF)
					netmigM <- netmigM/s * totmigy
					netmigF <- netmigF/s * totmigy
				}
			}
			this.all.migM <- cbind(this.all.migM, netmigM)
			this.all.migF <- cbind(this.all.migF, netmigF)
		}
		colnames(this.all.migM) <- colnames(this.all.migF) <- c('country_code', 'name', 'age', years.periods)
		all.migM <- rbind(all.migM, this.all.migM)
		all.migF <- rbind(all.migF, this.all.migF)
	}
	if(write.to.disk) {
		write.table(all.migM, file=file.path(directory, paste0(file.prefix, "M.txt")), sep='\t', row.names=FALSE)
		write.table(all.migF, file=file.path(directory, paste0(file.prefix, "F.txt")), sep='\t', row.names=FALSE)
		if(verbose) cat('\nMigration files written into ', file.path(directory, paste0(file.prefix, "X.txt")))
	}
	if(verbose) cat('\n')
	return(invisible(list(male=all.migM, female=all.migF)))
}

rcastro.schedule <- function(annual = FALSE) {
    if(annual) return(rcastro1.schedule())
    return(rcastro5.schedule())
}

rcastro5.schedule <- function()
    c(0.06133, 0.02667, 0.02067, 0.10467, 0.188, 
      0.18067, 0.13733, 0.09533, 0.064, 0.04267, 
      0.028, 0.01867, 0.012, 0.008, 0.00533, 
      0.00333, 0.00333, 0, 0, 0, 
      0)

rcastro1.schedule <- function()
    c(0.00734, 0.00734, 0.00734, 0.00694, 0.00654, 0.00614, 0.00573, 0.00533, 0.0051, 0.00486, 
      0.00461, 0.00437, 0.00413, 0.0075, 0.01085, 0.01421, 0.01758, 0.02093, 0.02456, 0.02819, 
      0.03182, 0.03545, 0.03908, 0.03888, 0.03869, 0.03849, 0.0383, 0.0381, 0.03627, 0.03444, 
      0.03261, 0.03078, 0.02894, 0.02697, 0.02499, 0.02302, 0.02104, 0.01907, 0.01782, 0.01656, 
      0.0153, 0.01405, 0.0128, 0.01195, 0.01109, 0.01024, 0.00939, 0.00854, 0.00795, 0.00736, 
      0.00677, 0.00619, 0.0056, 0.00522, 0.00486, 0.00448, 0.00411, 0.00373, 0.00347, 0.0032, 
      0.00293, 0.00267, 0.0024, 0.00224, 0.00208, 0.00192, 0.00176, 0.0016, 0.00149, 0.00139, 
      0.00128, 0.00117, 0.00107, 0.00098, 0.00091, 0.00083, 0.00074, 0.00067, 0.00067, 0.00067, 
      0.00067, 0.00067, 0.00067, 0.00053, 4e-04, 0.00027, 0.00013, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0)