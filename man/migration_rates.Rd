\name{migration_rates}
\alias{migration_rates}
\alias{migration_rates_wpp2012}
\alias{migration_rates_wpp2015}

\docType{data}
\title{
Datasets with Historical Data on Migration Rates
}

\description{
Datasets with migration rates for 5-year periods from 1950-2015 for all countries.}
\usage{
    data(migration_rates_wpp2015)
    data(migration_rates_wpp2012) 
}

\format{
  A data frame containing the following variables:
  \describe{
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{name}}{Country name.}
    \item{\code{1950-1955}, \code{1955-1960}, \dots}{Migration rates in various five-year time intervals.}
   }
}

%\source{
%Data provided by the United Nations Population Division.
%}

\examples{
data(land_area_wpp2015)
str(land_area_wpp2015)
}
\keyword{datasets}
