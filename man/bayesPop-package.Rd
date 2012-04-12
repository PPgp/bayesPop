\name{bayesPop-package}
\alias{bayesPop-package}
\alias{bayesPop}
\docType{package}
\title{
Probabilistic Population Projection
}
\description{
The package allows to generate population projections for all countries of the world using several probabilistic components, such as total fertility rate (TFR) and life expectancy.
}
\details{
\tabular{ll}{
Package: \tab bayesPop\cr
Type: \tab Package\cr
Version: \tab 0.2-1\cr
Date: \tab 2012-04-04\cr
License: \tab GPL (>=2)\cr
URL: \tab \url{http://bayespop.csss.washington.edu}\cr 
}
The main function is called \code{\link{pop.predict}}. It uses trajectories of TFR from the \pkg{bayesTFR} package and life expectancy from the \pkg{bayesLife} package and for each trajectory it computes a population projection using the Cohort component method. It results in probabilistic age and sex specific projections. Various plotting functions are available for results visualization (\code{\link{pop.trajectories.plot}}, \code{\link{pop.pyramid}}, \code{\link{pop.trajectories.pyramid}}), as well as a summary function (\code{\link{summary.bayesPop.prediction}}).
}
\author{
Hana Sevcikova, Adrian Raftery

Maintainer: Hana Sevcikova <hanas@uw.edu>
}
%\references{
%}

\keyword{ package }
\seealso{
\code{\link[bayesTFR:bayesTFR-package]{bayesTFR}}, \code{\link[bayesLife:bayesLife-package]{bayesLife}}
}
\examples{
sim.dir <- tempfile()
# Generates population projection for one country
country <- 'Netherlands'
pred <- pop.predict(countries=country, nr.traj=10, output.dir=sim.dir)
summary(pred, country)
pop.trajectories.plot(pred, country)
pop.pyramid(pred, country)
pop.pyramid(pred, country, year=2100, age=1:26)
unlink(sim.dir, recursive=TRUE)
}
