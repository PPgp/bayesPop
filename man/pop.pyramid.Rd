\name{pop.pyramid}
\alias{pop.pyramid}
\alias{pop.pyramid.all}
\alias{pop.trajectories.pyramid}
\alias{pop.trajectories.pyramid.all}
\docType{data}
\title{
Probabilistic Population Pyramid
}
\description{
Functions for plotting probabilistic population pyramid for any country and projection period.
}
\usage{
pop.pyramid(pop.pred, country, year = NULL, pi = c(80, 95), 
    proportion = FALSE, main = NULL, age = 1:21)
    
pop.pyramid.all(pop.pred, year = NULL,
    output.dir = file.path(getwd(), 'pop.pyramid'),
	output.type = "png", verbose = FALSE, \dots)
	
pop.trajectories.pyramid(pop.pred, country, year = NULL, pi = c(80, 95), 
    nr.traj = NULL, proportion = FALSE, main = NULL, age = 1:21)
    
pop.trajectories.pyramid.all(pop.pred, year = NULL,
    output.dir = file.path(getwd(), 'pop.traj.pyramid'),
	output.type = "png", verbose = FALSE, \dots)
}
\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Name or numerical code of a country.}
  \item{year}{Year within the projection period to be plotted. Default is the start year of the prediction. In the functions \code{pop.pyramid.all} and \code{pop.trajectories.pyramid.all}, the \code{year} argument can be a vector of years.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{proportion}{Logical. If \code{TRUE} the pyramid contains the distribution of rates of age-specific counts and population totals.}
  \item{main}{Titel of the plot. By default it is the country name and projection year.}
  \item{age}{Integer vector of age indices. Value 1 corresponds to age 0-4, value 2 corresponds to age 5-9 etc. Last available age goup is 130+ which corresponds to index 27. The purpose of this argument here is mainly to controls the top of the pyramid.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.}
  \item{\dots}{Arguments passed to the underlying functions.}
}
\details{
For a given country and projection year the \code{pop.pyramid} function generates a population pyramid that includes given confidence intervals. \code{pop.pyramid.all} creates such pyramids for all countries for which a projection is available and for all years given by the \code{year} argument.

Functions \code{pop.trajectories.pyramid} and \code{pop.trajectories.pyramid.all} work in a similar way. The only difference is that the pyramid is created using the population trajectories.   
}

\author{
Hana Sevcikova, Adrian Raftery, using feedback from Sam Clark and the bayesPop group.
}

\seealso{
\code{\link{pop.trajectories.plot}}, \code{\link{bayesPop.prediction}}, \code{\link{summary.bayesPop.prediction}}
}
\examples{
sim.dir <- file.path(.find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir)
pop.pyramid(pred, "Netherlands", 2045)
pop.trajectories.pyramid(pred, "Netherlands", 2045, age=1:25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

