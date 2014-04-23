\name{pop.map}
\alias{pop.map}
\alias{get.pop.map.parameters}
\alias{pop.map.gvis}

\title{
World Map of Population Measures
}
\description{
Generates a world map of various population measures for a given quantile and a projection or observed period.
}
\usage{
pop.map(pred, sex = c("both", "male", "female"), age = "all", expression = NULL, \dots)

get.pop.map.parameters(pred, expression = NULL, sex = c('both', 'male', 'female'), 
    age = 'all', range = NULL, nr.cats = 50, same.scale = TRUE, quantile = 0.5, \dots)
    
pop.map.gvis(pred, \dots)
}
\arguments{
  \item{pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{sex}{One of \dQuote{both}  (default), \dQuote{male} or \dQuote{female}. By default the male and female counts are summed up. This argument is only used if \code{expression} is \code{NULL}.}
  \item{age}{Either a character string \dQuote{all} (default) or an integer vector of age indices. Value 1 corresponds to age 0-4, value 2 corresponds to age 5-9 etc. Last age goup \eqn{130+} corresponds to index 27. This argument is only used if \code{expression} is \code{NULL}.}
  \item{expression}{Expression defining the population measure to be plotted. For syntax see \code{\link{pop.expressions}}. The country components of the expression should be given as \dQuote{XXX}.}
  \item{range}{Range of the population measure to be displayed. It is of the form \code{c(}\var{min}, \var{max}\code{)}.}
  \item{nr.cats}{Number of color categories.}
  \item{same.scale}{Logical controlling if maps for all years of this prediction object should be on the same color scale.}
   \item{quantile}{Quantile for which the map should be generated. It must be equal to one of the values in \code{dimnames(pred$quantiles[[2]])}, i.e. 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 1. Value 0.5 corresponds to the median.}
  \item{\dots}{Additional arguments passed to the underlying functions. In \code{pop.map}, these are \code{quantile}, \code{year}, \code{projection.index}, \code{device}, \code{main}, and \code{device.args} (see \code{\link[bayesTFR]{tfr.map}}). In \code{pop.map.gvis}, these are all arguments that can be passed to \code{\link[bayesTFR]{tfr.map.gvis}}. In addition, the first two functions accept
arguments passed to the \code{\link[rworldmap]{mapCountryData}} function of the \pkg{rworldmap} package.
}
}
\details{
\code{pop.map} creates a single map for the given time period and quantile. If the package \pkg{fields} is installed, a color bar legend at the botom of the map is created.

Function \code{get.pop.map.parameters} can be used in combination with \code{pop.map}. It sets breakpoints for the color scheme.

Function \code{pop.map.gvis} creates an interactive map using the \pkg{googleVis} package and opens it in an internet browser. It also generates a table of the mapped values that can be sorted by columns interactively in the browser. 
}

\section{Performance and Caching}{
	If the \code{expression} argument or a non-standard combination of sex and age is used, quantiles are computed on the fly. In such a case, trajectory files for all countries have to be loaded from disk, which can be quite time expensive. Therefore a simple caching mechanism was added to the prediction object which allows re-using data from previously used expressions. The prediction object points to an environment called \code{cache} which is a collection of data arrays that are results of evaluating expressions. The space-trimmed expressions are the names of the \code{cache} entries. Every time a map function is called, it is checked if the corresponding expression is contained in the \code{cache}. If it is not the case, the quantiles are computed on the fly, otherwise the existing values are taken. 
	
When computing on the fly, the function tries to process it in parallel if possible, using the package \pkg{parallel}. In such a case, the computation is split into \eqn{n} nodes where \eqn{n} is either the number of cores detected automatically (default), or the value of \code{getOption('cl.cores')}. Use \code{options(cl.cores=n)} to modify the default. If a sequential processing is desired, set \code{cl.cores} to 1. 

The cache data are also stored on disk, namely in the simulation directory of the prediction object. By default, every update of the cache in memory is also updated on the disk. Thus, data expression results can be re-used in multiple R sessions. Function \code{\link{pop.cleanup.cache}} deletes the content of the cache. This behaviour can be turned off by setting the argument \code{write.to.cache=FALSE} in the \code{\link{get.pop.prediction}} function. We use this settings in the examples throughout this manual whenever the example data from the installation directory is used, in order to prevent writing into the installation directory.
}


\value{
\code{get.pop.map.parameters} returns a list with elements:
  \item{pred}{The object of class \code{\link{bayesPop.prediction}} used in the function.}
  \item{quantile}{Value of the argument \code{quantile}.}
  \item{catMethod}{If the argument \code{same.scale} is \code{TRUE}, this element contains breakpoints for categorization. Otherwise, it is \code{NULL}.}
  \item{numCats}{Number of categories.}
  \item{coulourPalette}{Subset of the rainbow palette, starting from dark blue and ending at red.}
  \item{\dots}{Additional arguments passed to the function.}
}
\author{
Hana Sevcikova
}

\seealso{
\code{\link[bayesTFR]{tfr.map}}
}
\examples{
\dontrun{
##########################
# This example only makes sense if there is a simulation 
# for all countries. Below, only two countries are included,
# so the map is useless.
##########################
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir=sim.dir, write.to.cache=FALSE)
# Uses heat colors with seven categories by default
pop.map(pred, sex='female', age=4:10)
# Female population in child-bearing age as a proportion of totals
pop.map(pred, expression='PXXX_F[4:10] / PXXX')
# The same with more colors
params <- get.pop.map.parameters(pred, expression='PXXX_F[4:10] / PXXX')
do.call('pop.map', params)
# Another projection year on the same scale
do.call('pop.map', c(list(year=2043), params))
}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

