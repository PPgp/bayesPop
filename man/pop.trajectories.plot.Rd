\name{pop.trajectories.plot}
\alias{pop.trajectories.plot}
\alias{pop.trajectories.plotAll}
\alias{pop.trajectories.table}
\alias{pop.byage.plot}
\alias{pop.byage.plotAll}
\alias{pop.byage.table}

\title{
Output of Probabilistic Population Projection
}
\description{
The functions plot and tabulate the distribution of population projection for a given country, or for all countries, including the median and given probability intervals.
}
\usage{
pop.trajectories.plot(pop.pred, country = NULL, expression = NULL, pi = c(80, 95), 
    sex = c("both", "male", "female"), age = "all", sum.over.ages = TRUE, 
    half.child.variant = FALSE, nr.traj = NULL, typical.trajectory = FALSE,
    main = NULL, dev.ncol = 5, lwd = c(2, 2, 2, 2, 1), 
    col = c("black", "red", "red", "blue", "#00000020"), show.legend = TRUE, 
    ann = par("ann"), xshift = 0, \dots)
    
pop.trajectories.plotAll(pop.pred, 
    output.dir=file.path(getwd(), "pop.trajectories"),
    output.type="png", expression = NULL, verbose=FALSE, \dots)
    
pop.trajectories.table(pop.pred, country = NULL, expression = NULL, pi = c(80, 95), 
    sex = c("both", "male", "female"), age = "all", half.child.variant = FALSE,  
    xshift = 0, \dots)
    
pop.byage.plot(pop.pred, country = NULL, year = NULL, expression = NULL, 
    pi = c(80, 95), sex = c("both", "male", "female"), 
    half.child.variant = FALSE, nr.traj = NULL, typical.trajectory=FALSE,
    xlim = NULL, ylim = NULL, xlab = "", ylab = "Population projection", 
    main = NULL, lwd = c(2,2,2,1), col = c("red", "red", "blue", "#00000020"),
    show.legend = TRUE, add = FALSE, ann = par("ann"), type = "l", pch = NA, 
    pt.cex = 1, \dots)
    
pop.byage.plotAll(pop.pred, 
    output.dir=file.path(getwd(), "pop.byage"),
    output.type="png", expression = NULL, verbose=FALSE, \dots)

pop.byage.table(pop.pred, country = NULL, year = NULL, expression = NULL, 
    pi = c(80, 95), sex = c("both", "male", "female"), 
    half.child.variant = FALSE)
}
\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Name or numerical code of a country. It can also be given as ISO-2 or ISO-3 characters.}
  \item{expression}{Expression defining the population measure to be plotted. For syntax see \code{\link{pop.expressions}}. For \code{pop.trajectories.plot}, \code{pop.trajectories.table}, \code{pop.byage.plot} and \code{pop.byage.table} the basic components of the expression must be country-specific. For \code{pop.trajectories.plotAll} and \code{pop.byage.plotAll} the country part should be given as \dQuote{XXX}. In addition, expressions passed into \code{pop.byage.plot} and \code{pop.byage.table} must contain curly braces (i.e. be age specific).}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{sex}{One of \dQuote{both}  (default), \dQuote{male} or \dQuote{female}. By default the male and female projections are summed up.}
  \item{age}{Either a character string \dQuote{all} (default) or an integer vector of age indices. In a five year simulation, value 1 corresponds to age 0-4, value 2 corresponds to age 5-9 etc. Last age goup \eqn{130+} corresponds to index 27. In an annual simulation, the age indices 1, 2, 3, \dots, 131 corrrespond to ages 0, 1, 2, \dots, \eqn{130+}.}
  \item{sum.over.ages}{Logical.  If \code{TRUE}, the values are summed up over given age groups. Otherwise there is a separate plot for each age group.}
  \item{half.child.variant}{Logical. If TRUE the United Nations \dQuote{+/-0.5 child} variant computed with fertility \eqn{+/- 0.5*} TFR median and the median of life expectancy is shown.
}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.}
  \item{typical.trajectory}{Logical. If \code{TRUE} one trajectory is shown that has the smallest distance to the median.}
  \item{xlim, ylim, xlab, ylab, main, ann, pt.cex}{Graphical parameters passed to the \code{plot} function.}
  \item{xshift}{Constant added to the x-axis (year).}
  \item{dev.ncol}{Number of column for the graphics device if \code{sum.over.ages} is \code{FALSE}. If the number of age groups is smaller than \code{dev.ncol}, the number of columns is automatically decreased.}
    \item{lwd, col}{For the first three functions it is a vector of five elements giving the line width and color for: 1. observed data, 2. median, 3. quantiles, 4. half-child variant, 5. trajectories. For functions that show results by age it is a vector of four elements - as above without the first item (observed data).}
    \item{type, pch}{Currently works for plotting by age only. It is a vector of four elements giving the plot type and point type for: 1. median, 2. quantiles, 3. half-child variant, 4. trajectories. The last element of the array is recycled.}
    \item{show.legend}{Logical controlling whether the legend should be drawn.}
  \item{\dots}{Additional graphical arguments. Functions \code{pop.trajectories.plotAll} and \code{pop.byage.plotAll} accept also any arguments of \code{pop.trajectories.plot} and \code{pop.byage.plot}, respectively, except \code{country}.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{year}{Any year within the time period to be outputted.}
  \item{add}{Logical specifying if the plot should be added to an existing graphics.}
}
\details{
\code{pop.trajectories.plot} plots trajectories of population projection by time for a given country. \cr \code{pop.trajectories.table} gives the same output as a table. \code{pop.trajectories.plotAll} creates a set of graphs (one per country) that are stored in \code{output.dir}. The projections can be visualized separately for each sex and age groups, or summed up over both sexes and/or given age groups. This is controlled by the arguments \code{sex}, \code{age} and \code{sum.over.ages}.

\code{pop.byage.plot} and \code{pop.byage.table} plots/tabulate the posterior distribution by age for a given country and time period. \code{pop.byage.plotAll} creates such plots for all countries. 

The median and given probability intervals are computed using all available trajectories. Thus, \code{nr.traj} does not influence those values - it is used only to control the number of trajectories plotted.

If plotting results of an expression and the function fails, to debug obtain values of that expression using the functions \code{\link{get.pop.ex}} (for \code{pop.trajectories.plot}) and \code{\link{get.pop.exba}} (for \code{pop.byage.plot}).
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesPop.prediction}}, \code{\link{summary.bayesPop.prediction}}, \code{\link{pop.pyramid}}, \code{\link{pop.expressions}}, \code{\link{get.pop}}
}
\examples{
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir)
pop.trajectories.plot(pred, country="Ecuador", pi=c(80, 95))
pop.trajectories.table(pred, country="ECU", pi=c(80, 95))
# female population of Ecuador in child bearing ages (by time)
pop.trajectories.plot(pred, expression="PEC_F[4:10]") 
# Population by age in Netherands for two different years 
pop.byage.plot(pred, country="Netherlands", year=2050)
pop.byage.plot(pred, expression="PNL{}", year=2000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
