\name{write.pop.projection.summary}
\alias{write.pop.projection.summary}
\alias{write.pop.trajectories}

\title{
Writing Projection Summary and Trajectory Files
}
\description{
Functions for creating ASCII files containing projection summaries, such as the median, the lower and upper bound of the 80 and 95\% probability intervals, respectively, as well as containing individual trajectories.
}
\usage{
write.pop.projection.summary(pop.pred, what = NULL, expression = NULL, 
    output.dir = NULL, \dots)
    
write.pop.trajectories(pop.pred, expression = "PXXX", 
    output.file = "pop_trajectories.csv", byage = FALSE, 
    observed = FALSE,  wide = FALSE, digits = NULL,
    include.name = FALSE, sep = ",", na.rm = TRUE, \dots)
}

\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{what}{A character vector specifying what kind of projection to write. Total population is specified by \dQuote{pop}. Vital events are specified by \dQuote{births}, \dQuote{deaths}, \dQuote{sr} (survival rate), \dQuote{fertility} and \dQuote{pfertility} (percent fertility). Each of these strings can (some must) have a suffix \dQuote{sex} and/or \dQuote{age} if sex- and/or age-specific measure is desired. For example, \dQuote{popage}, \dQuote{birthssexage}, \dQuote{deaths}, \dQuote{deathssex}, are all valid values. Note that for survival, only \dQuote{srsexage} is allowed. For percent fertility, only \dQuote{pfertilityage} is allowed. Suffix \dQuote{sex} cannot be used in combination with \dQuote{fertility}. Moreover, \dQuote{fertility} (without age) corresponds to the total fertility rate. If the argument is \code{NULL}, all valid combinations are used.  The argument is not used if \code{expression} is given. Note that vital events can be only used if the prediction object contains vital events, i.e. if it was generated with the \code{keep.vital.events} argument being \code{TRUE} (see \code{\link{pop.predict}}).}
  \item{expression}{Expression defining the measure to be written. If it is not \code{NULL}, argument \code{what} is ignored. For expression syntax see \code{\link{pop.expressions}}. The country components of the expression should be given as \dQuote{XXX}.}
  \item{output.dir}{Directory in which the resulting files will be stored. If \code{NULL} \code{pop.pred$output.directory} is used.}
  \item{output.file}{File name to write the trajectories into.}
  \item{byage}{Logical indicating if the expression is defined by age, i.e. if it includes curly braces (\code{TRUE}), of if it is defined by time (\code{FALSE}), see \code{\link{pop.expressions}} for more detail on the expression syntax.}
  \item{observed}{Logical indicating if observed data should be written (\code{TRUE}) or projected trajectories (\code{FALSE}).}
  \item{wide}{Logical indicating if the data format should be wide. By default, trajectories are written in long format.}
  \item{digits}{To how many decimal digits should the indicator be rounded. By default no rounding takes place.}
  \item{include.name}{Logical indicating if country names should be included in the dataset.}
  \item{sep}{The field separator string.}
  \item{na.rm}{Logical indicating if records with \code{NA} values should be included in the dataset.}
  \item{\dots}{For \code{write.pop.projection.summary}, these are:
  \itemize{
  \item if \code{expression} is given, then one can use here \code{file.suffix} (defines the file suffix) and/or \code{expression.label} which defaults to the actual expression and is put as the first line in the resulting file;
  \item logical \code{include.observed} determines if observed data should be included;
  \item integer \code{digits} defines the number of decimal places in the resulting file;
  \item for 5-year projections, logical \code{end.time.only} determines if the time columns should be in form of time periods (as XXXX-YYYY) or just the end years (YYYY);
  \item logical \code{adjust} determines if the numbers should be adjusted; in such a case, \code{adj.to.file} and \code{allow.negative.adj} give the file name to which to adjust and a switch if negatives are allowed for the adjustments, respectively.
 }
  For \code{write.pop.trajectories}, these are arguments passed to \code{\link{get.pop.ex}} (if \code{byage} is \code{FALSE}) or \code{\link{get.pop.exba}} (if \code{byage} is \code{TRUE}). }
}
\details{
The \code{write.pop.projection.summary} function creates one file per value of \code{what}, or \code{expression}, called \file{projection_summary_}\emph{suffix}\file{.csv}, where \emph{suffix} is either \code{what} or, if an expression is given, the value of \code{file.suffix}. It is a comma-separated table with the following columns:
\itemize{
	\item \dQuote{country_name}: country name 
	\item \dQuote{country_code}: country code 
	\item \dQuote{variant}: name of the variant, such as \dQuote{median}, \dQuote{lower 80}, \dQuote{upper 80}, \dQuote{lower 95}, \dQuote{upper 95}
	\item period1: e.g. \dQuote{2005-2010}, or \dQuote{2010}: Given population measure for the first time period
	\item period2: e.g. \dQuote{2010-2015}, or \dQuote{2015}: Given population measure for the second time period
	\item \dots further time period columns	
}
If  \code{expression} is given, \code{expression.label} (by default the full expression) is written as the first line of the file starting with #. The file contains one line per country, and possibly sex and age.

Function \code{write.pop.trajectories} writes out all trajectories, either in long format (default) or, if \code{wide = TRUE} in wide format (years become columns).
}

\note{
If the \code{expression} argument is used, the same applies as for \code{\link{pop.map}} in terms of Performance and Caching.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{pop.predict}}, \code{\link{pop.map}}, \code{\link{pop.expressions}}
}
\examples{
outdir <- tempfile()
dir.create(outdir)
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir=sim.dir, write.to.cache=FALSE)

# proportion of 65+ years old to the whole population
write.pop.projection.summary(pred, expression="PXXX[14:27] / PXXX", file.suffix="age65plus", 
    output.dir=outdir, include.observed=TRUE, digits=2)
    
# various measures
write.pop.projection.summary(pred, what=c("pop", "popsexage", "popsex"),
    output.dir=outdir)

unlink(outdir, recursive=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }

