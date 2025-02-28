\name{peak.probability}
\alias{peak.probability}

\title{
Probability of Peaks in Population Indicators
}
\description{
For a given indicator and a country, the function computes the probability of a peak happening before a given year, as well as a range of years between which a peak happens 
with given probability.
}
\usage{
peak.probability(pop.pred, country = NULL, expression = NULL, year = NULL, 
    pi = 95, verbose = TRUE, \dots)
}
\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Name or numerical code or ISO-2 or ISO-3 character code of a country. If given, population is used as an indicator and the \code{expression} argument is ignored.}
  \item{expression}{Expression defining an indicator. For syntax see \code{\link{pop.expressions}}. It must be defined by time (i.e. either without or with square brackets, and no curly braces). Only used if \code{country} is not speicified.}
  \item{year}{Used for computing the probability of a peak happenning before \code{year}.}
  \item{pi}{Probability between 0 and 100. Used for selecting a range of years between which a peak happens with probability given by this argument.}
  \item{verbose}{Logical. If \code{TRUE}, results are printed.}
  \item{\dots}{Additional arguments passed to the underlying functions. If \code{country} is given, these are arguments passed to \code{\link{pop.trajectories}}, e.g. \code{sex}, \code{age} or \code{adjust}. If the indicator is given via \code{expression}, it can be e.g. \code{adj.to.file}.}
}
\details{
Given an indicator, the function computes two quantities:
\itemize{
\item{probability that the indicator reaches its peak before given \code{year};}
\item{range of years between which a peak happens with the given probability \code{pi}.}
}
The indicator can be either population (if \code{country} is given), or it can be any expression defined as a function of time (see \code{\link{pop.expressions}}).
}

\value{
List with elements:
\item{prob.peak.less.given.year}{Probability that the indicator reaches its peak before \code{year}.}
\item{given.year}{The value of \code{year}.}
\item{peak.quantiles}{The lower bound, the upper bound and the median of years defining a time interval in which a peak happens with the given probability \code{pi}}.
\item{all.prob.peak.by.time}{Data frame containing the probability of peak happening in each projected year, as well as the corresponding cummulative probability. 
Years in which no peak is projected are not included.}
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{pop.expressions}}
}
\examples{
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir, write.to.cache=FALSE)

# probability that population of Netherlands peaks before 2040 
# and between which years it will peak with probablity 80%
peak.probability(pred, "NL", year = 2040, pi = 80)

# check visually with  
# pop.trajectories.plot(pred, "NL")

# the same for female of age 45-49
peak.probability(pred, "NL", year = 2040, pi = 80, sex = "female", age = 10)

# probability of a peak for the potential support ratio in Ecuador
peak.probability(pred, expression = "PEC[5:13]/PEC[14:27]")

# check visually that it already peaked
# pop.trajectories.plot(pred, expression = "PEC[5:13]/PEC[14:27]")
}

% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ univar }

