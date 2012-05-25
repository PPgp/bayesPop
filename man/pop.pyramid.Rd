\name{pop.pyramid}
\alias{pop.pyramid}
\alias{pop.pyramidAll}
\alias{pop.pyramid.bayesPop.prediction}
\alias{pop.pyramid.bayesPop.pyramid}
\alias{pop.trajectories.pyramid}
\alias{pop.trajectories.pyramidAll}
\alias{pop.trajectories.pyramid.bayesPop.prediction}
\alias{pop.trajectories.pyramid.bayesPop.pyramid}
\alias{plot.bayesPop.pyramid}
\alias{bayesPop.pyramid}
\alias{get.bPop.pyramid}
\alias{get.bPop.pyramid.bayesPop.prediction}
\alias{get.bPop.pyramid.data.frame}
\alias{get.bPop.pyramid.matrix}
\alias{get.bPop.pyramid.list}
\docType{data}
\title{
Probabilistic Population Pyramid
}
\description{
Functions for plotting probabilistic population pyramid. \code{pop.pyramid} creates a classic pyramid using rectangles; \code{pop.trajectories.pyramid} creates a pyramid using vertical lines (possibly derived from population trajectories). They can be used for a prediction object created with this package, or for any user-defined sex- and age-specific dataset. In such a case, function \code{get.bPop.pyramid} should be used to translate user-defined data into a \code{bayesPop.pyramid} object.
}
\usage{
\method{pop.pyramid}{bayesPop.prediction}(pop.pred, country, year = NULL, pi = c(80, 95), 
    proportion = FALSE, age = 1:21, draw.past.year = FALSE, plot = TRUE, \dots)
    
\method{pop.pyramid}{bayesPop.pyramid}(pyramid, main = NULL, show.legend = TRUE, ann = TRUE, 
    pyr1.par = list(border='black', col=NA, density=NULL),
	pyr2.par = list(border='violet', col='violet', density=20), col.pi = NULL, \dots)
    
pop.pyramidAll(pop.pred, year = NULL,
    output.dir = file.path(getwd(), 'pop.pyramid'),
	output.type = "png", verbose = FALSE, \dots)
	
\method{pop.trajectories.pyramid}{bayesPop.prediction}(pop.pred, country, year = NULL, pi = c(80, 95), 
    nr.traj = NULL, proportion = FALSE, age = 1:21, draw.past.year = FALSE, plot = TRUE, \dots)
    
\method{pop.trajectories.pyramid}{bayesPop.pyramid}(pyramid, main = NULL, show.legend = TRUE, 
    ann = TRUE, col = c('red', 'red', 'black', 'gray'), lwd = c(2, 2, 2, 1), \dots)
    
pop.trajectories.pyramidAll(pop.pred, year = NULL,
    output.dir = file.path(getwd(), 'pop.traj.pyramid'),
	output.type = "png", verbose = FALSE, \dots)
	
\method{plot}{bayesPop.pyramid}(pyramid, \dots)

\method{get.bPop.pyramid}{bayesPop.prediction}(pop.pred, country, year = NULL, pi = c(80, 95), 
    proportion = FALSE, age = 1:21, nr.traj = 0, draw.past.year = FALSE, sort.pi=TRUE, \dots)
    
\method{get.bPop.pyramid}{data.frame}(data, main.label = NULL, legend = 'observed', 
    is.proportion = FALSE, ages = NULL, pop.max = NULL, \dots)
    
\method{get.bPop.pyramid}{matrix}(data, \dots)

\method{get.bPop.pyramid}{list}(data.list, main.label = NULL, 
    legend1 = 'observed', legend2 = 'comparison', \dots)
}
\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Name or numerical code of a country.}
  \item{year}{Year within the projection or estimation period to be plotted. Default is the start year of the prediction. In the functions \code{pop.pyramid.all} and \code{pop.trajectories.pyramid.all}, the \code{year} argument can be a vector of years.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{proportion}{Logical. If \code{TRUE} the pyramid contains the distribution of rates of age-specific counts and population totals.}
  \item{age}{Integer vector of age indices. Value 1 corresponds to age 0-4, value 2 corresponds to age 5-9 etc. Last available age goup is 130+ which corresponds to index 27. The purpose of this argument here is mainly to control the height of the pyramid.}
  \item{draw.past.year}{If \code{TRUE}, present year is shown in the pyramid. In addition to a logical value it can be a numeric value giving a past year to be drawn in the pyramid for comparison purposes.}
  \item{plot}{If \code{FALSE}, nothing is plotted. It can be used to retrieve the pyramid object without drawing it.}
  \item{pyramid}{Object of class \code{bayesPop.pyramid}. See below.}
  \item{main}{Titel of the plot. By default it is the country name and projection year if known.}
  \item{show.legend}{Logical controlling if the plot legend is drawn.}
  \item{ann}{Logical controlling if any annotation (main and legend) is plotted.}
  \item{pyr1.par, pyr2.par}{List of graphical parameters (color, border and density) for drawing the pyramid rectangles, for the first and second pyramid, respectively (see Details).}
  \item{col.pi}{Vector of colors for drawing the probability boxes. If it is given, it must be of the same length as \code{pi}.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.}
  \item{col, lwd}{Vector of size four giving the color/line width for 1. median, 2. quantiles, 3. comparison data, 4. trajectories.}
  \item{sort.pi}{Logical controlling if the probability intervals are sorted in decreasing order. This has an effect on the order in which they are plotted and thus on overlapping of pyramid boxes. By default the largest intervals are plotted first.}
  \item{data}{Data frame or matrix. It must have columns called \dQuote{male} and \dQuote{female}. The row names will determine the age labels, unless the \code{ages} argument is given.}
  \item{main.label, legend, legend1, legend2}{Optional arguments for the main title, legend of the main pyramimd and the comparison pyramid.}
  \item{is.proportion}{Either logical, indicating if the values in \code{data} are proportions, or \code{NA} in which case the proportions are computed.}
  \item{ages}{Vector of age labels. It must be of the same length as the number of rows of \code{data}. If it is not given, the age labels are considered to be the row names of \code{data}.}
  \item{pop.max}{Maximum value to be drawn in the pyramid. If it is not given, \code{max(data)} is taken.}
  \item{data.list}{List with one to three element. The first two are data frames or matrices. The first is used for the main pyramid, the second for a comparison pyramid. The third optional element is a list of quantiles. It must be in the same format as bayesPop$quantiles (see Section Value below).}
  \item{\dots}{Arguments passed to the underlying functions.}
}
\details{
For a given country and projection year the \code{pop.pyramid} and \code{pop.trajectories.pyramid} functions generate a population pyramid that includes given confidence intervals. \code{pop.pyramid.all} and \code{pop.trajectories.pyramid.all} create such pyramids for all countries for which a projection is available and for all years given by the \code{year} argument.

The functions can generate one or two pyramids in one plot. The first (main) one is usually the median of a future year prediction, but it can also be the current year or any population estimates. The second one serves the purpose of comparing two pyramids, and it must be data from the past (i.e. without uncertainty). It is drawn on top of the main pyramid. For example, one can use it to compare a future prediction with the present, or two different time points in the past, or two different geographies.

The core of these functions operates on a \code{bayesPop.pyramid} object which is automatically created when called with a \code{bayesPop.prediction} object. If used with a user-defined data set, one has to convert such data into \code{bayesPop.pyramid} using the function \code{get.bPop.pyramid} (see an example below). In such a case, one can simply use the \code{plot} function which then calls \code{pop.pyramid}.
}

\value{
	\code{pop.pyramid}, \code{pop.trajectories.pyramid} and \code{get.bPop.pyramid} return an object of class \code{bayesPop.pyramid} which is a list with the following components:
	\item{label}{Label used for the main titel.}
	\item{pyr1}{List with elements \code{legend} (text used in the legend) and \code{data} which is a data frame with columns \dQuote{female} and \dQuote{male} and one row per age group. Row names are used for labeling the y-axis. This component is used for drawing the main pyramid.}
	\item{pyr2}{The same structure as \code{pyr1}. It is used for drawing a comparison pyramid.}
	\item{quantiles}{List with one element per probability interval whose names are the values of the intervals. Each element is again a list with components \code{low} and \code{high} which have the same structure as \code{pyr1$data} and contain the lower and upper bounds of the corresponding interval.}
	\item{trajectories}{List with elements \dQuote{female} and \dQuote{male} each of which is a matrix of size number of age categories x number of trajectories. This is only used by the \code{pop.trajectories.pyramid} function.}
	\item{is.proportion}{Logical indicating if values in the various data frames in this object are proportions or raw values.}
	\item{pop.max}{Maximum value for the x-axis.}
}

\author{
Hana Sevcikova, Adrian Raftery, using feedback from Sam Clark and the bayesPop group.
}

\seealso{
\code{\link{pop.trajectories.plot}}, \code{\link{bayesPop.prediction}}, \code{\link{summary.bayesPop.prediction}}
}
\examples{
# pyramids for bayesPop prediction objects
##########################################
sim.dir <- file.path(.find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir)
pop.pyramid(pred, "Netherlands", 2045, draw.past.year=TRUE)
dev.new()
pop.trajectories.pyramid(pred, "Netherlands", 2045, age=1:25, draw.past.year=TRUE)
# using manual manipulation of the data: e.g. show only the prob. intervals 
pred.pyr <- get.bPop.pyramid(pred, country="Madagascar", year=2090, age=1:27)
pred.pyr$pyr1 <- NULL
plot(pred.pyr)

# pyramids for user-defined data
################################
# this example dataset conains population estimates for the Washington state and King county (Seattle area) in 2011
data <- read.table(file.path(.find.package("bayesPop"), "ex-data", "popestimates_WAKing.txt"), header=TRUE, row.names=1)
# extract data for two pyramids and put it into the right format
WA <- data[,c('WA.male', 'WA.female')]
colnames(WA) <- c('male', 'female')
King <- data[,c('King.male', 'King.female')]
colnames(King) <- c('male', 'female')
# create and plot a bayesPop.pyramid object
pyramid <- get.bPop.pyramid(list(WA, King), legend1='Washington', legend2='King')
plot(pyramid, main='Population in 2011')
# show data as proportions
pyramid.prop <- get.bPop.pyramid(list(WA, King), is.proportion=NA, legend1='Washington', legend2='King')
pop.pyramid(pyramid.prop, pyr1.par=list(col='lightgreen', border='lightgreen', density=30), 
    pyr2.par=list(col='darkred', border='darkred', density=50), main='Population in 2011 (proportions)')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

