\name{summary.bayesPop.prediction}
\alias{summary.bayesPop.prediction}
\alias{print.summary.bayesPop.prediction}

\title{
Summary of Probabilistic Population Projection
}
\description{
Summary of an object \code{\link{bayesPop.prediction}} created using the \code{\link{pop.predict}} function. The summary contains the mean, standard deviation and several commonly used quantiles of the simulated trajectories.
}
\usage{
\method{summary}{bayesPop.prediction}(object, country = NULL, 
    sex = c("both", "male", "female"), compact = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Country name or code. It can also be given as ISO-2 or ISO-3 characters. If it is \code{NULL}, only meta information included.}
  \item{sex}{One of \dQuote{both} (default), \dQuote{male}, or \dQuote{female}. If it is not \dQuote{both}, the summary is given for sex-specific trajectories.}
  \item{compact}{Logical switching between a smaller and larger number of displayed quantiles.}
  \item{\dots}{A list of further arguments.}
}

\author{
Hana Sevcikova
}


\seealso{
\code{\link{bayesPop.prediction}}
}
\examples{
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir)
summary(pred, "Netherlands")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ print }
