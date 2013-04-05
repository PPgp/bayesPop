\name{write.pop.projection.summary}
\alias{write.pop.projection.summary}

\title{
Writing Projection Summary Files
}
\description{
The function creates ASCII files containing projection summaries, such as the median, the lower and upper bound of the 80 and 90\% probability intervals, respectively. 
}
\usage{
write.pop.projection.summary(pop.pred, what = NULL, expression = NULL, 
    output.dir = NULL, \dots)
}

\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}}.}
  \item{what}{A character vector specifying what kind of projection to write. Total population is specified by \dQuote{pop}. Vital events are specified by \dQuote{births}, \dQuote{deaths}, \dQuote{sr} (survival rate) and \dQuote{fertility}. Each of these strings can have a suffix \dQuote{sex} and/or \dQuote{age} if sex- and/or age-specific measure is desired. For example, \dQuote{popage}, \dQuote{birthssexage}, \dQuote{deaths}, \dQuote{deathssex}, are all valid values. Note that for survival, only \dQuote{srsexage} is allowed and fertility does not allow the \dQuote{sex} suffix. If the argument is \code{NULL}, all valid combinations are used.  The argument is not used if \code{expression} is given. Note that vital events can be only used if the prediction object contains vital events, i.e. if it was generated with the \code{keep.vital.events} argument being \code{TRUE} (see \code{\link{pop.predict}}).}
  \item{expression}{Expression defining the measure to be written. If it is not \code{NULL}, argument \code{what} is ignored. For expression syntax see \code{\link{pop.expressions}}. The country components of the expression should be given as \dQuote{XXX}.}
  \item{output.dir}{Directory in which the resulting files will be stored. If \code{NULL} \code{pop.pred$output.directory} is used.}
  \item{\dots}{These are arguments used if \code{expression} is given: \code{file.suffix} defines the file suffix; \code{expression.label} (defaults to the expression) is put as the first line in the resulting file; logical \code{include.observed} determines if observed data should be included; integer \code{decimal} defines the number of decimal places in the resulting file.}
}
\details{
There is one file created per value of \code{what}, or \code{expression}, called \file{projection_summary_}\emph{suffix}\file{.csv}, where \emph{suffix} is either \code{what} or \code{alias} of the expression. It is a comma-separated table with the following columns:
\itemize{
	\item 	\dQuote{country_name}: country name 
	\item \dQuote{country_code}: country code 
	\item \dQuote{variant}: name of the variant, such as \dQuote{median}, \dQuote{lower 80}, \dQuote{upper 80}, \dQuote{lower 95}, \dQuote{upper 95}
	\item period1: e.g. \dQuote{2005-2010}: Given population measure for the first time period
	\item period2: e.g. \dQuote{2010-2015}: Given population measure for the second time period
	\item \dots further time period columns	
}
If  \code{expression} is given, the full expression is written as the first line of the file starting with \#.
The file contains one line per country, and possibly sex and age.
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
write.pop.projection.summary(pred, expression="PXXX[14:27] / PXXX", file.suffix='age65plus', 
    output.dir=outdir, include.observed=TRUE, decimal=2)
    
# various measures
write.pop.projection.summary(pred, what=c('pop', 'popsexage', 'popsex'),
    output.dir=outdir)

unlink(outdir, recursive=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }

