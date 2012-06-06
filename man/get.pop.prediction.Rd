\name{get.pop.prediction}
\alias{get.pop.prediction}
\alias{has.pop.prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessing Prediction Object
}
\description{
Function \code{get.pop.prediction} retrieves results of a prediction from disk and creates an object of class \code{\link{bayesPop.prediction}}. Function \code{has.pop.prediction} checks an existence of such results.
}
\usage{
get.pop.prediction(sim.dir, aggregation = NULL)

has.pop.prediction(sim.dir)
}

\arguments{
  \item{sim.dir}{Directory where the prediction is stored. It should correspond to the value of the \code{output.dir} argument used in the \code{\link{pop.predict}} function.}
  \item{aggregation}{If given, the prediction object is considered to be an aggregation and both arguments are passed to \code{\link{get.pop.aggregation}}.}
}
\details{
The \code{\link{pop.predict}} function stores resulting trajectories into a directory called \code{output.dir}/prediction. Here the argument \code{sim.dir} should correspond to \code{output.dir} (i.e. without the \dQuote{prediction} part).
}
\value{
Function \code{has.pop.prediction} returns a logical indicating if a prediction exists. 
        
Function \code{get.pop.prediction} returns an
object of class \code{\link{bayesPop.prediction}}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesPop.prediction}}, \code{\link{get.pop.aggregation}}
}
\examples{
sim.dir <- file.path(.find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir)
summary(pred)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

