\name{vwBaseYear}
\alias{vwBaseYear}
\alias{vwBaseYear2010}
\alias{vwBaseYear2012}
\alias{vwBaseYear2015}
\docType{data}
\title{
Datasets on Migration Base Year and Type, and Mortality and Fertility Age Patterns
}
\description{
Datasets giving information on the baseyear and type of migration for each country. The 2012 and 2015 datasets also give information on country's specifics regarding mortality and fertility age patterns.
}
\usage{
    data(vwBaseYear2015)
    data(vwBaseYear2012) 
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
    \item{\code{WPPAIDS}}{Dummy indicating if the country has generilized HIV/AIDS epidemics.}
    \item{\code{AgeMortalityType}}{Type of mortality age pattern. Only relevant for countries with the entry \dQuote{Model life tables}. In such a case, the \eqn{b_x} Lee-Carter parameter is not estimated from historical data. Instead is taken from the dataset \code{\link{MLTbx}} using a pattern given in the \code{AgeMortalityPattern} column.}
    \item{\code{AgeMortalityPattern}}{If \code{AgeMortalityType} is equal to \dQuote{Model life tables}, this value determines which \eqn{b_x} is selected from the \code{\link{MLTbx}} dataset. It must sorrespond to one of the rownames of \code{\link{MLTbx}}, e.g. \dQuote{CD East}, \dQuote{CD West}, \dQuote{UN Latin American}.}
    \item{\code{LatestAgeMortalityPattern}}{Dummy indicating if the \eqn{a_x} Lee-Carter parameter should be computed using observed data from the latest time period only.}
    \item{\code{SmoothLatestAgeMortalityPattern}}{If \code{LatestAgeMortalityPattern} is 1, this column indicates if the \eqn{a_x} should be smoothed.}
    \item{\code{PasfrNorm}}{Type of norm for computing age-specific fertility patern to which the country belongs to. Currently only \dQuote{GlobalNorm} is used.}
    \item{\code{PasfrGlobalNorm, PasfrFarEastAsianNorm, PasfrSouthAsianNorm}}{Dummies indicating which country to include to compute the specific norms.}
  }
}
\details{
There is one record for each country. See Sevcikova et al (2015) on how information from the various columns is used for projections.
}
\source{
Data provided by the United Nations Population Division.
}
\references{
H. Sevcikova, N. Li, V. Kantorova, P. Gerland and A. E. Raftery (2015). Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. arXiv:1503.05215. \url{http://arxiv.org/abs/1503.05215}
}
\examples{
data(vwBaseYear2015)
str(vwBaseYear2015)
}
\keyword{datasets}
