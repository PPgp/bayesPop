\name{LOCATIONS}
\alias{LOCATIONS}
\docType{data}
\title{
United Nations Table of Locations
}
\description{
United Nations table of locations, including regions. It is stored in a TAB-separated ASCII file \dQuote{LOCATIONS.txt}. 
}
\usage{data(LOCATIONS)}
\format{
  A data frame with one observations per country or region. It contains the following 9 variables:
  \describe{
    \item{\code{Index}}{Numeric vector enumerating the observations. Not used by the package.}
    \item{\code{name}}{Name of country or region (following ISO 3166 official short names in English - see \cr \url{http://www.iso.org/iso/country_codes/iso_3166_code_lists/english_country_names_and_code_elements.htm} and United Nations Multilingual Terminology Database - see \url{http://unterm.un.org}).}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{reg_code}}{Code of the regions.}
    \item{\code{reg_name}}{Name of the regions.}
    \item{\code{area_code}}{Area code.}
    \item{\code{area_name}}{Area names, such as  \code{Africa}, \code{Asia}, \code{Europe}, \code{World}.}
    \item{\code{location_type}}{Code giving the type of the observation (0=World, 2=Major Area, 3=Region, 4=Country/Area, 5=Development group, 12=Special groupings). This column is used to determine countries for aggregation in the \code{\link{pop.aggregate}} function.}
    \item{\code{include_code}}{Specific code. Not used by the package.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Data provided by the United Nations Population Division
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(LOCATIONS)
str(LOCATIONS)
}
\keyword{datasets}
