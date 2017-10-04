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
Functions for plotting probabilistic population pyramid. \code{pop.pyramid} creates a classic pyramid using rectangles; \code{pop.trajectories.pyramid} creates one or more pyramids using vertical lines (possibly derived from population trajectories). They can be used to view a prediction object created with this package, or any user-defined sex- and age-specific dataset. For the latter, function \code{get.bPop.pyramid} should be used to translate user-defined data into a \code{bayesPop.pyramid} object.
}
\usage{
\method{pop.pyramid}{bayesPop.prediction}(pop.object, country, year = NULL, 
    indicator = c("P", "B", "D"), pi = c(80, 95), 
    proportion = FALSE, age = 1:21, plot = TRUE, pop.max = NULL, \dots)
    
\method{pop.pyramid}{bayesPop.pyramid}(pop.object, main = NULL, show.legend = TRUE, 
    pyr1.par = list(border="black", col=NA, density=NULL, height=0.9),
    pyr2.par = list(density = -1, height = 0.3), 
    col.pi = NULL, ann = par("ann"), axes = TRUE, grid = TRUE, 
    cex.main = 0.9, cex.sub = 1, cex = 1, cex.axis = 1, \dots)
    
pop.pyramidAll(pop.pred, year = NULL,
    output.dir = file.path(getwd(), "pop.pyramid"),
    output.type = "png", one.file = FALSE, verbose = FALSE, \dots)
	
\method{pop.trajectories.pyramid}{bayesPop.prediction}(pop.object, country, year = NULL, 
    indicator = c("P", "B", "D"), pi = c(80, 95), nr.traj = NULL, 
    proportion = FALSE, age = 1:21, plot = TRUE, pop.max = NULL, \dots)
    
\method{pop.trajectories.pyramid}{bayesPop.pyramid}(pop.object, main = NULL, show.legend = TRUE, 
    col = rainbow, col.traj = "#00000020", lwd = 2, ann = par("ann"), axes = TRUE, 
    grid = TRUE, cex.main = 0.9, cex.sub = 1, cex = 1, cex.axis = 1, \dots)
    
pop.trajectories.pyramidAll(pop.pred, year = NULL,
    output.dir = file.path(getwd(), "pop.traj.pyramid"),
    output.type = "png", one.file = FALSE, verbose = FALSE, \dots)
	
\method{plot}{bayesPop.pyramid}(x, \dots)

\method{get.bPop.pyramid}{bayesPop.prediction}(data, country, year = NULL, 
    indicator = c("P", "B", "D"), pi = c(80, 95), 
    proportion = FALSE, age = 1:21, nr.traj = 0, sort.pi=TRUE, pop.max = NULL, \dots)
    
\method{get.bPop.pyramid}{data.frame}(data, main.label = NULL, legend = "observed", 
    is.proportion = FALSE, ages = NULL, pop.max = NULL, 
    LRmain = c("Male", "Female"), LRcolnames = c("male", "female"), CI = NULL, \dots)
    
\method{get.bPop.pyramid}{matrix}(data, \dots)

\method{get.bPop.pyramid}{list}(data, main.label = NULL, legend = NULL, CI = NULL, \dots)
}
\arguments{
  \item{pop.object}{Object of class \code{\link{bayesPop.prediction}} or \code{bayesPop.pyramid} (see Value section).}
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{x}{Object of class \code{bayesPop.pyramid}.}
  \item{data}{Data frame, matrix, list or object of class \code{\link{bayesPop.prediction}}. For data frame and matrix, it must have columns defined by \code{LRcolnames} (\dQuote{male} and \dQuote{female} by default). The row names will determine the age labels. For lists, it can be a collection of such data frames. The names of the list elements are used for legend, unless \code{legend} is given.}
  \item{country}{Name or numerical code of a country.}
  \item{year}{Year within the projection or estimation period to be plotted. Default is the start year of the prediction. It can also be a vector of years. \code{pop.pyramid} draws the first two, \code{pop.trajectories.pyramid} draws all of them.
  	In the functions \code{pop.pyramidAll} and \code{pop.trajectories.pyramidAll}, the \code{year} argument can be a list of years, in which case the pyramids are created for all elements in the list.}
  	\item{indicator}{One of the characters \dQuote{P} (population), \dQuote{B} (births), \dQuote{D} (deaths) determining the pyramid indicator.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{proportion}{Logical. If \code{TRUE} the pyramid contains the distribution of rates of age-specific counts and population totals.}
  \item{age}{Integer vector of age indices. Value 1 corresponds to age 0-4, value 2 corresponds to age 5-9 etc. Last available age goup is 130+ which corresponds to index 27. The purpose of this argument here is mainly to control the height of the pyramid.}
  \item{plot}{If \code{FALSE}, nothing is plotted. It can be used to retrieve the pyramid object without drawing it.}
  \item{main}{Titel of the plot. By default it is the country name and projection year if known.}
  \item{show.legend}{Logical controlling if the plot legend is drawn.}
  \item{pyr1.par, pyr2.par}{List of graphical parameters (color, border, density and height) for drawing the pyramid rectangles, for the first and second pyramid, respectively (see Details). The \code{height} component should be a number between 0 (corresponds to a line) and 1 (for non-overelapping rectangles). If \code{density} is \code{NULL}, the rectangles are transparent, see the argument \code{density} in \code{\link{rect}}.}
  \item{col.pi}{Vector of colors for drawing the probability boxes. If it is given, it must be of the same length as \code{pi}.}
   \item{ann}{Logical controlling if any annotation (main and legend) is plotted.}
   \item{axes}{Logical controlling if axes are plotted.}
   \item{grid}{Logical controlling if grid lines are plotted.}
   \item{cex.main, cex.sub, cex, cex.axis}{Magnification to be used for the title, secondary titles on the right and left panels, legend and axes, respectively.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{one.file}{Logical. If \code{TRUE} the output is put into one single file, by default a PDF.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.}
  \item{col}{Colors generating function. It is called with an argument giving the number of pyramids to be plotted. Each color is then used for one pyramid, including its confidence intervals.}
  \item{col.traj}{Color used for trajectories. If more than one pyramid is drawn with its trajectories, this can be a vector of the size of number of pyramids.}
  \item{lwd}{Line width for the pyramids.}
  \item{sort.pi}{Logical controlling if the probability intervals are sorted in decreasing order. This has an effect on the order in which they are plotted and thus on overlapping of pyramid boxes. By default the largest intervals are plotted first.}
  \item{main.label}{Optional argument for the main title.}
  \item{legend}{Legend to be used. In case of multiple pyramids, this can be a vector for each of them. If not given and \code{data} is a list, names of the list elements are taken as legend.}
  \item{is.proportion}{Either logical, indicating if the values in \code{data} are proportions, or \code{NA} in which case the proportions are computed.}
  \item{ages}{Vector of age labels. It must be of the same length as the number of rows of \code{data}. If it is not given, the age labels are considered to be the row names of \code{data}.}
  \item{pop.max}{Maximum value to be drawn in the pyramid. If it is not given, \code{max(data)} is taken.}
  \item{LRmain}{Vector of character strings giving the secondary titles for the left and right panel, respectively.}
  \item{LRcolnames}{Vector of character strings giving the column names of data to be used for the left and right panel of the pyramid, respectively.}
  \item{CI}{Confidence intervals. It should be of the same format as the \code{bayesPop.pyramid$CI} object, see below.}
  \item{\dots}{Arguments passed to the underlying functions.}
}
\details{
The \code{pop.pyramid} function generates one or two population pyramids in one plot. The first (main) one is usually the median of a future year prediction, but it can also be the current year or any population estimates. The second one serves the purpose of comparing two pyramids with one another and is drawn on top of the main pyramid. For example, one can use it to compare a future prediction with the present, or two different time points in the past, or two different geographies. The main pyramid can have confidence intervals associated with it, which are also plotted. If \code{pop.pyramid} is called on a \code{bayesPop.prediction} object, the main and secondary pyramid, respectively, is generated from data of a time period given by the first and second element, respectively, of the \code{year} argument. In such a case, confidence intervals only of the first year are shown. Thus, it makes sense to set the first year to be a prediction year and the second year to an observed time period. If \code{pop.pyramid} is called on a \code{bayesPop.pyramid} object, data in the first and second element, respectively, of the \code{bayesPop.pyramid$pyramid} list are used, and only the first element of \code{bayesPop.pyramid$CI} is used.

Pyramids generated via the \code{pop.trajectories.pyramid} function have different appearance and therefore more than two pyramids can be put into one figure. Furthermore, confidence intervals of more than one pyramid can be shown. Thus, all elements of \code{bayesPop.pyramid$pyramid} and \code{bayesPop.pyramid$CI} are plotted. In addition, single trajectories given in \code{bayesPop.pyramid$trajectories} can be shown by setting the argument \code{nr.traj} larger than 0. 

Both, \code{pop.pyramid} and \code{pop.trajectories.pyramid} (if called with a \code{bayesPop.prediction} object) use data from one country. 
Functions \code{pop.pyramidAll} and \code{pop.trajectories.pyramidAll} create such pyramids for all countries for which a projection is available and for all years given by the \code{year} argument which should be a list. In this case, one pyramid figure (possibly containing multiple pyramids) is created for each country and each element of the \code{year} list.

The core of these functions operates on a \code{bayesPop.pyramid} object which is automatically created when called with a \code{bayesPop.prediction} object. If used with a user-defined data set, one has to convert such data into \code{bayesPop.pyramid} using the function \code{get.bPop.pyramid} (see an example below). In such a case, one can simply use the \code{plot} function which then calls \code{pop.pyramid}.
}

\value{
	\code{pop.pyramid}, \code{pop.trajectories.pyramid} and \code{get.bPop.pyramid} return an object of class \code{bayesPop.pyramid} which is a list with the following components:
	\item{label}{Label used for the main titel.}
	\item{pyramid}{List of pyramid data, one element per pyramid. Each component is a data frame with at least two columns, containing data for the left and right panels of the pyramid. Their names must correspond to \code{LRcolnames} (see below). There is one row per age group and the row names are used for labeling the y-axis. Names of the list elements are used in the legend.}
	\item{CI}{List of lists of confidence intervals with one element per pyramid. The order corresponds to the order in the \code{pyramid} component and it is \code{NULL} if the corresponding pyramid does not have confidence intervals. Each element is a list with one element per probability interval whose names are the values of the intervals. Each element is again a list with components \code{low} and \code{high} which have the same structure as \code{pyramid} and contain the lower and upper bounds of the corresponding interval.}
	\item{trajectories}{List of lists of trajectories with one element per pyramid. As in the case of \code{CI}, it is ordered the same way as the \code{pyramid} component and is \code{NULL} if the corresponding pyramid does not have any trajectories to be shown. Each element is again a list with two components, one for the left part and one for the right part of the pyramid. Their names correspond to \code{LRcolnames} and each of them is a matrix of size number of age categories x number of trajectories. This is only used by the \code{pop.trajectories.pyramid} function.}
	\item{is.proportion}{Logical indicating if values in the various data frames in this object are proportions or raw values.}
	\item{pop.max}{Maximum value for the x-axis.}
	\item{LRmain}{Vector of character strings determining the titles for the left and right panels, respectively.}
	\item{LRcolnames}{Vector of character strings determining the column names in \code{pyramid}, \code{CI} and \code{trajectories} used to plot data into the left and right panel, respectively.}
}

\author{
Hana Sevcikova, Adrian Raftery, using feedback from Sam Clark and the bayesPop group at the University of Washington.
}

\seealso{
\code{\link{pop.trajectories.plot}}, \code{\link{bayesPop.prediction}}, \code{\link{summary.bayesPop.prediction}}
}
\examples{
# pyramids for bayesPop prediction objects
##########################################
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir)
pop.pyramid(pred, "Netherlands", c(2045, 2010))
dev.new()
pop.trajectories.pyramid(pred, "Netherlands", c(2045, 2010, 1960), age=1:25, proportion=TRUE)
# using manual manipulation of the data: e.g. show only the prob. intervals 
pred.pyr <- get.bPop.pyramid(pred, country="Ecuador", year=2090, age=1:27)
pred.pyr$pyramid <- NULL
plot(pred.pyr)

# pyramids for user-defined data
################################
# this example dataset contains population estimates for the Washington state and King county 
# (Seattle area) in 2011
data <- read.table(file.path(find.package("bayesPop"), "ex-data", "popestimates_WAKing.txt"), 
    header=TRUE, row.names=1)
# extract data for two pyramids and put it into the right format
head(data)
WA <- data[,c("WA.male", "WA.female")]; colnames(WA) <- c("male", "female")
King <- data[,c("King.male", "King.female")]; colnames(King) <- c("male", "female")
# create and plot a bayesPop.pyramid object
pyramid <- get.bPop.pyramid(list(WA, King), legend=c("Washington", "King"))
plot(pyramid, main="Population in 2011", pyr2.par=list(height=0.7, col="violet", border="violet"))
# show data as proportions
pyramid.prop <- get.bPop.pyramid(list(WA, King), is.proportion=NA, legend=c("Washington", "King"))
pop.pyramid(pyramid.prop, main="Population in 2011 (proportions)",
    pyr1.par=list(col="lightgreen", border="lightgreen", density=30), 
    pyr2.par=list(col="darkred", border="darkred", density=50))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

