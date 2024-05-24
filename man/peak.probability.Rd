\name{peak.probability}
\alias{peak.probability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability of Peaks in Projected Indicators
}
\description{
For a given indicator in a country, the function computes the probability of a peak happening before given year and the earliest year a peak happens with given probability.
}
\usage{
peak.probability(pop.pred, country = NULL, expression = NULL, year = NULL, 
    prob = 0.5, verbose = TRUE, \dots)
}
\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Name or numerical code or ISO-2 or ISO-3 character code of a country. If given, population is used as an indicator and the \code{expression} argument is ignored.}
  \item{expression}{Expression defining the indicator. For syntax see \code{\link{pop.expressions}}. It must be define by time (i.e. either without or with brackets, and no curly braces). Only used if \code{country} is not speicified.}
  \item{year}{Integer used for computing the probability of the peak happenning before \code{year}.}
  \item{prob}{Probability between 0 and 1 to be used for selecting the earliest year that the peak happens.}
  \item{verbose}{Logical. If \code{TRUE}, results are printed.}
  \item{\dots}{Additional arguments passed to the underlying functions. If \code{country} is given, these are arguments passed to \code{\link{pop.trajectories}}, e.g. \code{sex}, \code{age} or \code{adjust}. If the indicator is given via \code{expression}, it can be e.g. \code{adj.to.file}.}
}
\details{
Given a projected indicator, the function computes two quantities:
\itemize{
\item{probability that the indicator reaches its peak before given \code{year};}
\item{the earliest year by which a peak happens with the given probability \code{prob}.}
}
The indicator can be either population (if \code{country} is given), or it can be any expression defined as a function of time (see \code{\link{pop.expressions}}).
}

\value{
List with elements:
\item{prob.peak.less.given.year}{Probability that the indicator reaches its peak before \code{year}.}
\item{given.year}{The value of \code{year}.}
\item{min.peak.year.given.prob}{The earliest year by which a peak happens with the given probability given by \code{prob}}.
\item{given.prob}{The value of \code{prob}.}
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
# and in which year it will peak the earliest with probablity 0.6
peak.probability(pred, "NL", year = 2040, prob = 0.6)

# check visually with  
# pop.trajectories.plot(pred, "NL")

# probability of a peak for the potential support ratio in Ecuador
peak.probability(pred, expression = "PEC[5:13]/PEC[14:27]")

# check visually (it already peaked!) 
# pop.trajectories.plot(pred, expression = "PEC[5:13]/PEC[14:27]")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
