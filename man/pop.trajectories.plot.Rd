\name{pop.trajectories.plot}
\alias{pop.trajectories.plot}
\alias{pop.trajectories.plot.all}
\alias{pop.trajectories.table}

\title{
Output of Probabilistic Population Projection
}
\description{
The functions plot and tabulate the distribution of population projection for a given country, or for all countries, including the median and given probability intervals.
}
\usage{
pop.trajectories.plot(pop.pred, country, pi = c(80, 95), 
    sex = c("both", "male", "female"), age = "all", sum.over.ages = FALSE, 
    half.child.variant = FALSE, nr.traj = NULL, typical.trajectory = FALSE,
    xlim = NULL, ylim = NULL, xlab = "Year", ylab = "Population projection", 
    main = NULL, dev.ncol = 5, lwd = c(2, 2, 2, 2, 1), 
    col = c('black', 'red', 'red', 'blue', 'gray'), show.legend = TRUE, \dots)
    
pop.trajectories.plot.all(pop.pred, 
    output.dir=file.path(getwd(), "pop.trajectories"),
    output.type="png", verbose=FALSE, \dots)
    
pop.trajectories.table(pop.pred, country, pi = c(80, 95), 
    sex = c("both", "male", "female"), age = "all",
    half.child.variant = FALSE)
}
\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Name or numerical code of a country.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{sex}{One of \dQuote{both}  (default), \dQuote{male} or \dQuote{female}. By default the male and female projections are summed up.}
  \item{age}{Either a character string \dQuote{all} (default) or an integer vector of age indices. Value 1 corresponds to age 0-4, value 2 corresponds to age 5-9 etc. Last age goup is \eqn{130+} which corresponds to index 27. The default value of \dQuote{all} corresponds to indices \code{1:20}, i.e. ages from 0 to 100.}
  \item{sum.over.ages}{Logical.  If \code{TRUE}, the values are summed up over given age groups. Otherwise there is a separate plot for each age group.}
  \item{half.child.variant}{Logical. If TRUE the United Nations \dQuote{+/-0.5 child} variant computed with fertility \eqn{+/- 0.5*} TFR median and the median of life expectancy is shown.
}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.}
  \item{typical.trajectory}{Logical. If \code{TRUE} one trajectory is shown that has the smallest distance to the median.}
  \item{xlim, ylim, xlab, ylab, main}{Graphical parameters passed to the \code{plot} function.}
  \item{dev.ncol}{Number of column for the graphics device if \code{sum.over.ages} is \code{FALSE}. If the number of age groups is smaller than \code{dev.ncol}, the number of columns is automatically decreased.}
    \item{lwd, col}{Vector of five elements giving the line width and color for: 1. observed data, 2. median, 3. quantiles, 4. half-child variant, 5. trajectories.}
    \item{show.legend}{Logical controlling whether the legend should be drawn.}
  \item{\dots}{Additional graphical arguments. Function \code{pop.trajectories.plot.all} accepts also any arguments of \code{pop.trajectories.plot} (except \code{country}).}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
}
\details{
\code{pop.trajectories.plot} plots trajectories of population projection for a given country. \cr \code{pop.trajectories.table} gives the same output as a table. \code{pop.trajectories.plot.all} creates a set of graphs (one per country) that are stored in \code{output.dir}. The projections can be visualized separately for each sex and age groups, or summed up over both sexes and/or given age groups. This is controlled by the arguments \code{sex}, \code{age} and \code{sum.over.ages}.
The median and given probability intervals are computed using all available trajectories. Thus, \code{nr.traj} does not influence those values - it is used only to control the number of trajectories plotted.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesPop.prediction}}, \code{\link{summary.bayesPop.prediction}}, \code{\link{pop.pyramid}}
}
\examples{
sim.dir <- file.path(.find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir)
pop.trajectories.plot(pred, country="Malaysia", pi=c(80, 95))
pop.trajectories.table(pred, country="Malaysia", pi=c(80, 95))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
