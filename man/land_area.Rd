\name{land_area}
\alias{land_area}
\alias{land_area_wpp2012}
\alias{land_area_wpp2015}

\docType{data}
\title{
Datasets on Countries Land Area
}

\description{
Datasets of land area in square kilometers for all countries.}
\usage{
    data(land_area_wpp2015)
    data(land_area_wpp2012) 
}

\format{
  A data frame containing the following variables:
  \describe{
  	\item{\code{name}}{Country name.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{land_area}}{Land area in square kilometers for each country.}
   }
}

\source{
Data provided by the United Nations Population Division.
}
\examples{
data(land_area_wpp2015)
str(land_area_wpp2015)
}
\keyword{datasets}
