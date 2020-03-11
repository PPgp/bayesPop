\name{vwBaseYear}
\alias{vwBaseYear}
\alias{vwBaseYear2010}
\alias{vwBaseYear2012}
\alias{vwBaseYear2015}
\alias{vwBaseYear2017}
\alias{vwBaseYear2019}
\docType{data}
\title{
Datasets on Migration Base Year and Type, and Mortality and Fertility Age Patterns
}
\description{
Datasets giving information on the baseyear and type of migration for each country. The 2012, 2015, 2017 and 2019 datasets also give information on country's specifics regarding mortality and fertility age patterns.
}
\usage{
    data(vwBaseYear2019)
    data(vwBaseYear2017)
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
    \item{\code{AgeMortProjMethod1}}{Method for projecting age-specific mortality rates. It is one of \dQuote{LC} (modified Lee-Carter, uses function \code{\link[MortCast]{mortcast}}), \dQuote{PMD} (pattern mortality decline, uses function \code{\link[MortCast]{copmd}}), \dQuote{MLT} (model life tables, uses function \code{\link[MortCast]{mlt}}), \dQuote{LogQuad} (log quadratic method, uses function \code{\link[MortCast]{logquad}}), or \dQuote{HIVmortmod} (HIV model life tables as implemented in the \pkg{HIV.LifeTables} package which can be installed from the \href{https://github.com/PPgP/HIV.LifeTables}{PPgP/HIV.LifeTables} GitHub repo).}
    \item{\code{AgeMortProjMethod2}}{If the mortality rates are to be projected via a blend of two methods (see \code{\link[MortCast]{mortcast.blend}}), this column determines the second method. The options are the same as in the column \code{AgeMortProjMethod1}.}
    \item{\code{AgeMortProjPattern}}{If one of the \code{AgeMortProjMethodX} colums contains the \dQuote{MLT} method, this column determines the type of the life table (see the argument \code{type} in the \code{\link[MortCast]{mlt}} function).}
    \item{\code{AgeMortProjMethodWeights}}{If the mortality rates are to be projected via a blend of two methods, this column determines the weights in the first and the last year of the projection, respectively. It should be given as an R vector, e.g. \dQuote{c(1, 0.5)} (see the argument \code{weights} in \code{\link[MortCast]{mortcast.blend}}).}
    \item{\code{AgeMortProjAdjSR}}{Code determining how the \dQuote{PMD} method should be adjusted if it's used. 0 means no adjustment, 1 means the argument \code{sexratio.adjust} in \code{\link[MortCast]{copmd}} is set to \code{TRUE}, and code 3 means that the argument \code{adjust.sr.if.needed} in \code{\link[MortCast]{copmd}} is set to \code{TRUE}.}
    \item{\code{LatestAgeMortalityPattern}}{Indicator \eqn{n} for how many latest time periods should be used to compute  the \eqn{a_x} Lee-Carter parameter. If \eqn{n} is zero, all time periods are used. If \eqn{n} is one, only the latest time period is used. If \eqn{n} is negative, the latest \eqn{n} time periods are excluded.}
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
data(vwBaseYear2019)
str(vwBaseYear2019)
}
\keyword{datasets}
