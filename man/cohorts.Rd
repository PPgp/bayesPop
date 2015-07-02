\name{pop.cohorts.plot}
\alias{pop.cohorts.plot}
\alias{cohorts}

\title{
Extracting and Plotting Cohort Data 
}
\description{
Extracts and plots population counts or results of expressions by cohorts. 
}
\usage{
cohorts(pop.pred, country = NULL, expression = NULL, pi = c(80, 95))
	
pop.cohorts.plot(pop.pred, country = NULL, expression = NULL, cohorts = NULL, 
    cohort.data = NULL, pi = c(80, 95), dev.ncol = 5, show.legend = TRUE, 
    legend.pos = "bottomleft", ann = par("ann"), add = FALSE, xlab = "", ylab = "",  
    main = NULL, xlim = NULL, ylim = NULL, col = "red", \dots)
}

\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{country}{Name or numerical code of a country. If it is not given, \code{expression} must be specified.}
  \item{expression}{Expression defining the population measure to be plotted. For syntax see \code{\link{pop.expressions}}. It must be country-specific, i.e. \dQuote{XXX} is not allowed, and it must contain curly braces, i.e. be age specific.}
   \item{pi}{Probability interval. It can be a single number or an array.}
  \item{cohorts}{Years of the cohorts to be plotted. By default, 10 future cohorts (starting from the last observed one) are used. It can be a single number or an array.}
  \item{cohort.data}{List with the cohort data obtained via the \code{cohorts} function. If it is not given, function \code{cohorts} is called internally, but by passing this argument the processing is faster.}
  \item{dev.ncol}{Number of column for the graphics device.}
  \item{show.legend}{Logical controlling whether the legend should be drawn.}
  \item{legend.pos}{Position of the legend passed to the \code{\link{legend}} function.}
  \item{ann, xlab, ylab, main, xlim, ylim, col, \dots}{Graphical parameters passed to the \code{plot} function.}
  \item{add}{Logical specifying if the plot should be added to an existing graphics.}
}
\details{
\code{pop.cohorts.plot} plots all cohorts passed in the \code{cohorts} argument on the same scale of the \eqn{y}-axis. 
}
\value{
Function \code{cohorts} returns a list where each element corresponds to one cohort. Each cohort element is a matrix with columns corresponding to years and rows corresponding to the median (first row) and quantiles of the given probability intervals.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pop.trajectories.plot}}, \code{\link{pop.byage.plot}}, \code{\link{pop.expressions}}
}
\examples{
    sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
    pred <- get.pop.prediction(sim.dir)
    # Population cohorts
    pop.cohorts.plot(pred, "Netherlands")
    # plot specific cohorts using expression (must contain {})
    pop.cohorts.plot(pred, expression="P528{}", cohorts=c(1960, 1980, 2000, 2020))
    # the same as
    cohort.data <- cohorts(pred, expression="P528{}")
    pop.cohorts.plot(pred, cohort.data=cohort.data, cohorts=c(1960, 1980, 2000, 2020))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ manip }
