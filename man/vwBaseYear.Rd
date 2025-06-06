\name{vwBaseYear}
\alias{vwBaseYear}
\alias{vwBaseYear2010}
\alias{vwBaseYear2012}
\alias{vwBaseYear2015}
\alias{vwBaseYear2017}
\alias{vwBaseYear2019}
\alias{vwBaseYear2022}
\alias{vwBaseYear2024}
\docType{data}
\title{
Datasets on Migration Base Year and Type, and Mortality and Fertility Age Patterns
}
\description{
Datasets giving information on the baseyear and type of migration for each country. The 2012, 2015, 2017, 2019, 2022 and 2024 datasets also give information on country's specifics regarding mortality, fertility and migration age patterns.
}
\usage{
    data(vwBaseYear2024)
    data(vwBaseYear2022)
    data(vwBaseYear2019)
    data(vwBaseYear2017)
    data(vwBaseYear2015)
    data(vwBaseYear2012) 
    data(vwBaseYear2010)  
}
\format{
  A data frame containing the following variables:
  \describe{
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{https://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{country}}{Country name. Not used by the package.}
    \item{\code{isSmall}}{UN internal code. Not used by the package.}
    \item{\code{ProjFirstYear}}{The base year of migration.}
    \item{\code{MigCode}}{Type of migration. Zero means migration is evenly distributed over each time interval. Code 9 means migration is captured at the end of each interval.}
    \item{\code{WPPAIDS}}{Dummy indicating if the country has generalized HIV/AIDS epidemics.}
    \item{\code{AgeMortalityType}}{Type of mortality age pattern. Only relevant for countries with the entry \dQuote{Model life tables}. In such a case, the \eqn{b_x} Lee-Carter parameter is not estimated from historical data. Instead is taken from the dataset \code{\link{MLTbx}} using a pattern given in the \code{AgeMortalityPattern} column.}
    \item{\code{AgeMortalityPattern}}{If \code{AgeMortalityType} is equal to \dQuote{Model life tables}, this value determines which \eqn{b_x} is selected from the \code{\link{MLTbx}} dataset. It must sorrespond to one of the rownames of \code{\link{MLTbx}}, e.g. \dQuote{CD East}, \dQuote{CD West}, \dQuote{UN Latin American}.}
    \item{\code{AgeMortProjMethod1}}{Method for projecting age-specific mortality rates. It is one of \dQuote{LC} (modified Lee-Carter, uses function \code{\link[MortCast]{mortcast}}), \dQuote{PMD} (pattern mortality decline, uses function \code{\link[MortCast]{copmd}}), \dQuote{modPMD} (modified pattern mortality decline, uses function \code{\link[MortCast]{copmd}(\dots use.modpmd = TRUE)}), \dQuote{MLT} (model life tables, uses function \code{\link[MortCast]{mlt}}), \dQuote{LogQuad} (log quadratic method, uses function \code{\link[MortCast]{logquad}}), or \dQuote{HIVmortmod} (HIV model life tables as implemented in the \pkg{HIV.LifeTables} package which can be installed from the \href{https://github.com/PPgP/HIV.LifeTables}{PPgP/HIV.LifeTables} GitHub repo).}
    \item{\code{AgeMortProjMethod2}}{If the mortality rates are to be projected via a blend of two methods (see \code{\link[MortCast]{mortcast.blend}}), this column determines the second method. The options are the same as in the column \code{AgeMortProjMethod1}.}
    \item{\code{AgeMortProjPattern}}{If one of the \code{AgeMortProjMethodX} colums contains the \dQuote{MLT} method, this column determines the type of the life table (see the argument \code{type} in the \code{\link[MortCast]{mlt}} function).}
    \item{\code{AgeMortProjMethodWeights}}{If the mortality rates are to be projected via a blend of two methods, this column determines the weights in the first and the last year of the projection, respectively. It should be given as an R vector, e.g. \dQuote{c(1, 0.5)} (see the argument \code{weights} in \code{\link[MortCast]{mortcast.blend}}).}
    \item{\code{AgeMortProjAdjSR}}{Code determining how the \dQuote{PMD} method should be adjusted if it's used. 0 means no adjustment, 1 means the argument \code{sexratio.adjust} in \code{\link[MortCast]{copmd}} is set to \code{TRUE}, and code 3 means that the argument \code{adjust.sr.if.needed} in \code{\link[MortCast]{copmd}} is set to \code{TRUE}.}
    \item{\code{LatestAgeMortalityPattern}, \code{LatestAgeMortalityPattern1}}{Indicator \eqn{n} for how many latest time periods of historical mortality rates should be averaged to compute  the \eqn{a_x} Lee-Carter and modPMD parameter. If \eqn{n} is zero, all time periods are used. If \eqn{n} is one, only the latest time period is used. If \eqn{n} is negative, the latest \eqn{n} time periods are excluded. This can have also a form of a vector where the first element is either a negative or a zero. If it is negative, the vector must have only two elements. In such a case, the first element (must be negative) determines how many latest time periods should be excluded, while the second element (must be positive) determines how many latest time periods to include after the exclusion. If the vector starts with a zero, the following numbers are interpreted as individual indices to the time periods starting from the latest time point. Here are a few examples, assuming the available mortality rates are on annual scale, from 1950 to 2023:
    \describe{
    \item{\dQuote{0}: }{ using all years from 1950 to 2023}
    \item{\dQuote{3}: }{ using 2023, 2022, 2021} 
    \item{\dQuote{-3}: }{ using 1950 - 2020}
    \item{\dQuote{c(-2, 3)}: }{ 2023 and 2022 are excluded; using 2021, 2020, 2019}
    \item{\dQuote{c(-2, 1, 3)}: }{ invalid specification - must have two elements if it starts with a negative}
    \item{\dQuote{c(0, 3)}: }{ interpreted as an individual index; thus, using 2021 only}
    \item{\dQuote{c(0, 1, 3, 4)}: }{ interpreted as individual indices; using 2023, 2021, 2020}
    }
    If the \code{LatestAgeMortalityPattern1} column is present, it should contain values related to an annual simulation (1x1) while the \code{LatestAgeMortalityPattern} column relates to a 5x5 simulation.}
    \item{\code{SmoothLatestAgeMortalityPattern}}{If \code{LatestAgeMortalityPattern} is not zero, this column indicates if the \eqn{a_x} should be smoothed.}
    \item{\code{SmoothDFLatestAgeMortalityPattern}, \code{SmoothDFLatestAgeMortalityPattern1}}{Degree of freedom for smoothing \eqn{a_x}. By default (value 0) a half of the number of age groups is taken. If the \code{SmoothDFLatestAgeMortalityPattern1} column is present, it should contain values related to a 1x1 simulation while the \code{SmoothDFLatestAgeMortalityPattern} column relates to a 5x5 simulation.}
    \item{\code{PasfrNorm}}{Type of norm for computing age-specific fertility pattern to which the country belongs to. Currently only \dQuote{GlobalNorm} is used.}
    \item{\code{PasfrGlobalNorm, PasfrFarEastAsianNorm, PasfrSouthAsianNorm}}{Dummies indicating which country to include to compute the specific norms.}
    \item{MigFDMb0, MigFDMb1, MigFDMmin, MigFDMsrin, MigFDMsrout}{Available in the 2024 dataset. These are parameters of the Flow Difference Method to generate age-specific net migration patterns (Sevcikova et. al, 2024). They correspond to the intercept, slope, minimum flow rate, female sex ratio for the in-flow and out-flow, respectively.}
  }
}
\details{
There is one record for each country. See Sevcikova et al (2016) on how information from the various columns is used for projections.
}
\source{
Data provided by the United Nations Population Division.
}
\references{
H. Sevcikova, N. Li, V. Kantorova, P. Gerland and A. E. Raftery (2016). Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. In: Dynamic Demographic Analysis, ed. Schoen R. (Springer), pp. 285-310. Earlier version in \href{https://arxiv.org/abs/1503.05215}{arXiv:1503.05215}.

H. Sevcikova, J. Raymer J., A. E. Raftery (2024). Forecasting Net Migration By Age: The Flow-Difference Approach. \href{https://arxiv.org/abs/2411.09878}{arXiv:2411.09878}.
}
\examples{
data(vwBaseYear2019)
str(vwBaseYear2019)
}
\keyword{datasets}
