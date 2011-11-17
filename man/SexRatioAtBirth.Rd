\name{SexRatioAtBirth}
\alias{SexRatioAtBirth}
\alias{SexRatioAtBirth2008}
\alias{SexRatioAtBirth2010}
\docType{data}
\title{
Datasets on Sex Ratio at Birth
}
\description{
Datasets giving predictions about the sex ratio at birth.
}
\usage{
    data(SexRatioAtBirth2008)
    data(SexRatioAtBirth2010)
}
\format{
  A data frame containing the following variables:
  \describe{
    \item{\code{index}}{Numerical index of the records. Not used by the package.}
    \item{\code{country}}{Country name. Not used by the package.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{2005-2010}, \code{2010-2015}, \dots}{Sex ratio at birth (females/males) for the given time period. Periods between \code{present.year} and \code{end.year} of the \code{\link{pop.predict}} function are used.}
  }
}
\details{
There is one record per country.
}
\source{
Data provided by the United Nations Population Division.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(SexRatioAtBirth2008)
str(SexRatioAtBirth2008)
}
\keyword{datasets}
