\name{ASMR}
\alias{ASMR}
\alias{ASMRMale}
\alias{ASMRMale2008}
\alias{ASMRMale2010}
\alias{ASMRFemale}
\alias{ASMRFemale2008}
\alias{ASMRFemale2010}
\docType{data}
\title{
Observed Mortality Data 
}
\description{
Age and sex specific historical data on mortality. Each dataset is connected to a WPP year which is part of the dataset name.
}
\usage{
data(ASMRMale2008)
data(ASMRMale2010)

data(ASMRFemale2008)
data(ASMRFemale2010)
}
\format{
  Data frames containing the following variables:
  \describe{
    \item{\code{index}}{Numerical index of the records. Not used by the package.}
    \item{\code{country}}{Country name. Not used by the package.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{age}}{A character string representing an age interval (given by the starting age of the interval). For each country there are 22 values: \dQuote{0}, \dQuote{1}, \dQuote{5}, \dQuote{10}, \dQuote{15}, \dQuote{20}, \dQuote{25}, \dQuote{30}, \dQuote{35}, \dQuote{40}, \dQuote{45},  \dQuote{50}, \dQuote{55}, \dQuote{60}, \dQuote{65}, \dQuote{70}, \dQuote{75}, \dQuote{80}, \dQuote{85}, \dQuote{90}, \dQuote{95}, and \dQuote{100+} in that order.}
    \item{\code{1950-1955}, \code{1955-1960}, \dots}{Mortality rate for the given time period. Only periods between \code{start.year} and \code{present.year} of the \code{\link{pop.predict}} function are used.}

  }
}
\details{
There is one record for each country and each age. Thus, for each country there are 22 records. The country records must be groupped together with age records following the order above. 
}
\source{
Data provided by the United Nations Population Division.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(ASMRFemale2010)
str(ASMRFemale2010)
}
\keyword{datasets}
