has.pop.prediction <- function(sim.dir) {
	if(file.exists(file.path(sim.dir, 'predictions', 'prediction.rda'))) return(TRUE)
	return(FALSE)
}

get.pop.prediction <- function(sim.dir) {
	############
	# Returns an object of class bayesPop.prediction
	############
	output.dir <- file.path(sim.dir, 'predictions')
	pred.file <- file.path(output.dir, 'prediction.rda')
	if(!file.exists(pred.file)) {
		warning('File ', pred.file, ' does not exist.')
		return(NULL)
	}
	load(file=pred.file)
	bayesPop.prediction$output.directory <- output.dir
	return(bayesPop.prediction)
}

summary.bayesPop.prediction <- function(object, country=NULL, sex=c('both', 'male', 'female'), compact=TRUE, ...) {
	if(!is.null(country)) country <- get.country.object(country, country.table=object$countries)
	sex <- match.arg(sex)
	if (sex == 'male') {
		object$quantiles <- object$quantilesM
		object$traj.mean.sd <- object$traj.mean.sdM
	} else {
		if (sex == 'female') {
			object$quantiles <- object$quantilesF
			object$traj.mean.sd <- object$traj.mean.sdF
		} else sex <- 'both'
	}
	res <- bayesTFR:::get.prediction.summary.data(object, 
				unchanged.pars=c('nr.traj', 'estim.years'), 
				country=country, compact=compact)
	res$nr.countries <- nrow(object$countries)
	res$sex <- sex
	class(res) <- 'summary.bayesPop.prediction'
	return(res)
}

print.summary.bayesPop.prediction <- function(x, digits = 5, ...) {
	cat('\nProjections:', length(x$projection.years), '(', x$projection.years[1], '-', 
					x$projection.years[length(x$projection.years)], ')')
	cat('\nInitial time point:', x$estim.years[length(x$estim.years)])
	cat('\nObserved time points:', length(x$estim.years), '(', x$estim.years[1], '-', 
					x$estim.years[length(x$estim.years)], ')')
	cat('\nTrajectories:', x$nr.traj)
	cat('\nNumber of countries:', x$nr.countries)
	if(!is.null(x$country.name)) {
		cat('\n\nCountry:', x$country.name, '\n')
		cat('\nProjected Population')
		if (x$sex != 'all') cat(' for', x$sex)
		cat(':\n')
		print(x$projections, digits=digits, ...)
	}
}


get.countries.table.bayesPop.prediction <- function(object, ...) 
	return(object$countries)