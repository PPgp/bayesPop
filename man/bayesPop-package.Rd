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
Version: \tab 6.0-0\cr
Date: \tab 2016-1-15\cr
License: \tab GPL (>=2)\cr
URL: \tab \url{http://bayespop.csss.washington.edu}\cr 
}
The main function is called \code{\link{pop.predict}}. It uses trajectories of TFR from the \pkg{bayesTFR} package and life expectancy from the \pkg{bayesLife} package and for each trajectory it computes a population projection using the Cohort component method. It results in probabilistic age and sex specific projections. Various plotting functions are available for results visualization (\code{\link{pop.trajectories.plot}}, \code{\link{pop.pyramid}}, \code{\link{pop.trajectories.pyramid}}), as well as a summary function (\code{\link{summary.bayesPop.prediction}}).
}


\author{
Hana Sevcikova, Adrian Raftery, Thomas Buettner

Maintainer: Hana Sevcikova <hanas@uw.edu>
}
\references{
A. E. Raftery,  N. Li, H. Sevcikova, P. Gerland, G. K. Heilig (2012). Bayesian probabilistic population projections for all countries. Proceedings of the National Academy of Sciences 109:13915-13921.

P. Gerland, A. E. Raftery, H. Sevcikova, N. Li, D. Gu, T. Spoorenberg, L. Alkema, B. K. Fosdick, J. L. Chunn, N. Lalic, G. Bay, T. Buettner, G. K. Heilig,  J. Wilmoth (2014). World Population Stabilization Unlikely This Century. Science 346:234-237.

H. Sevcikova, N. Li, V. Kantorova, P. Gerland and A. E. Raftery (2015). Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. arXiv:1503.05215. \url{http://arxiv.org/abs/1503.05215}
%Working Paper 150. Center for Statistics and the Social Sciences, University of Washington. \url{http://www.csss.washington.edu/Papers/2015/wp150.pdf} 
}


\keyword{ package }
\seealso{
\code{\link[bayesTFR:bayesTFR-package]{bayesTFR}}, \code{\link[bayesLife:bayesLife-package]{bayesLife}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
# Generates population projection for one country
country <- 'Netherlands'
pred <- pop.predict(countries=country, output.dir=sim.dir)
summary(pred, country)
pop.trajectories.plot(pred, country)
pop.pyramid(pred, country)
pop.pyramid(pred, country, year=2100, age=1:26)
unlink(sim.dir, recursive=TRUE)
}

# Here are commands needed to run probabilistic projections
# from scratch, i.e. including TFR and life expectancy.
# Note that running the first four commands 
# (i.e. predicting TFR and life expectancy) can take 
# LONG time (up to several days; see below for possible speed-up). 
# For a toy simulation, set the number of iterations (iter) 
# to a small number.
\dontrun{
sim.dir.tfr <- 'directory/for/TFR'
sim.dir.e0 <-  'directory/for/e0'
sim.dir.pop <- 'directory/for/pop'

# Estimate TFR parameters (speed-up by including parallel=TRUE)
run.tfr.mcmc(iter='auto', output.dir=sim.dir.tfr, seed=1)

# Predict TFR (if iter above < 4000, reduce burnin and nr.traj accordingly)
tfr.predict(sim.dir=sim.dir.tfr, nr.traj=2000, burnin=2000)

# Estimate e0 parameters (females) (speed-up by including parallel=TRUE)
# Can be run independently of the two commands above
run.e0.mcmc(sex='F', iter='auto', output.dir=sim.dir.e0, seed=1)

# Predict female and male e0	
# (if iter above < 22000, reduce burnin and nr.traj accordingly)
e0.predict(sim.dir=sim.dir.e0, nr.traj=2000, burnin=20000)

# Population prediction
pred <- pop.predict(output.dir=sim.dir.pop, verbose=TRUE, 
    inputs = list(tfr.sim.dir=sim.dir.tfr, 
                  e0F.sim.dir=sim.dir.e0, e0M.sim.dir='joint_'))
pop.trajectories.plot(pred, 'Madagascar', nr.traj=50, sum.over.ages=TRUE)
pop.trajectories.table(pred, 'Madagascar')
}
}
% reproducing example data in the package:
% pred <- pop.predict(countries=c(528,218), nr.traj=5, output.dir=sim.dir, replace.output=TRUE)
