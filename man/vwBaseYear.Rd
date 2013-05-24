\name{vwBaseYear}
\alias{vwBaseYear}
\alias{vwBaseYear2008}
\alias{vwBaseYear2010}
\alias{vwBaseYear2012}
\docType{data}
\title{
Datasets on Migration Base Year and Type
}
\description{
Datasets giving information on the baseyear and type of migration for each country.
}
\usage{
    data(vwBaseYear2008)
    data(vwBaseYear2010)   
}
\format{
  A data frame containing the following variables:
  \describe{
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{country}}{Country name. Not used by the package.}
    \item{\code{isSmall}}{UN internal code. Not used by the package.}
    \item{\code{ProjFirstYear}}{The base year of migration.}
    \item{\code{MigCode}}{Type of migration. Zero means migration is evenly distributed over each time interval. Code 9 means migration is captured at the end of each interval.}
  }
}
\details{
There is one record for each country.
}
\source{
Data provided by the United Nations Population Division.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(vwBaseYear2010)
str(vwBaseYear2010)
}
\keyword{datasets}
