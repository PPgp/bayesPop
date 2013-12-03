\name{get.pop}
\alias{get.pop}
\title{
Evaluation of Basic Components of Expressions
}
\description{
Evaluates a basic component of an expression and results in a four-dimensional array.
}
\usage{
get.pop(object, pop.pred, aggregation = NULL, observed = FALSE, ...)
}
\arguments{
  \item{object}{Character string giving the basic component of an expression (see \link{pop.expressions}).}
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{aggregation}{If the basic component is to be evaluated on an aggregated prediction object, this argument gives the name of the aggregation (corresponds argument \code{name} in \code{\link{pop.aggregate}}). By default, the function searches for available aggregations and gives priority to the one called \dQuote{independence}.}
  \item{observed}{Logical. Determines if the evaluation uses observed data (\code{TRUE}) or predictions (\code{FALSE}).}
  \item{\dots}{Additional argument passed to the underlying functions. Currently, this is only used for a \dQuote{P}opulation component and \code{observed=FALSE}. It can be either \code{nr.traj} giving the number of trajectories or logical \code{typical.trajectory}.}
}

\details{
This function is rather an internal function used for evaluation after an expression is decomposed into basic components. However, it can be also used for debugging purposes, to obtain results from parts of an expression.
}
\note{If \code{object} contains the wildcard \dQuote{XXX}, the function only works on observed data, i.e. \code{observed} must be \code{TRUE}.}
\value{
Array of four dimensions (country x age x time x trajectory). See \link{pop.expressions} for more details.
}

\author{
Hana Sevcikova
}

\seealso{
\link{pop.expressions}
}
\examples{
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir, write.to.cache=FALSE)

# observed female of Netherlands by age; 1x21x13x1 array
popFNL <- get.pop("PNL_F{}", pred, observed=TRUE)

# observed migration for all countries in the prediction object,
# here 2 countries; 2x1x13x1 array
migAll <- get.pop("GXXX", pred, observed=TRUE)

# projection population for Ecuador with 3 trajectories; 
# 1x1x19x3 array
popEcu <- get.pop("P218", pred, observed=FALSE)

# Expression "PNL_F{} / PNL_M{}" evaluated on projections
# is internally replaced by
FtoM <- get.pop("PNL_F{}", pred) / get.pop("PNL_M{}", pred)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ manip }

