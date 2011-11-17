\name{PercentASFR}
\alias{PercentASFR}
\alias{PercentASFR2008}
\alias{PercentASFR2010}
\docType{data}
\title{
Datasets on Age-specific Distribution of Fertility Projections
}
\description{
Datasets giving the percentage of fertility over ages 15-50. 
}
\usage{
    data(PercentASFR2008)
    data(PercentASFR2010)
}
\format{
  A data frame containing the following variables.
  \describe{
    \item{\code{index}}{Numerical index of the records. Not used by the package.}
    \item{\code{country}}{Country name. Not used by the package.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{age}}{A character string representing an age interval. For each country there are 7 values: \code{15-19} \code{20-24} \code{25-29} \code{30-34} \code{35-39} \code{40-44} \code{45-49}.}
    \item{\code{2005-2010}, \code{2010-2015}, \dots}{Fertility percentage for the given time period. Periods between \code{present.year} and \code{end.year} of the \code{\link{pop.predict}} function are used.}
  }
}
\details{
There is one record for each country and each age. Thus, for each country there are 7 records. The country records must be groupped together with age records following the order above.
}
\source{
Data provided by the United Nations Population Division.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(PercentASFR2010)
str(PercentASFR2010)
}
\keyword{datasets}
