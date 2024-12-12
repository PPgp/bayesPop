\name{rcFDM}
\alias{rc1FDM}
\alias{rc5FDM}
\alias{rcFDM}
\docType{data}

\title{
Datasets on Inflow and Outflow Migration Schedules for FDM Method
}
\description{
Age-specific schedules of the inflow and outflow migration distribution used as input for the FDM method. \code{rc1FDM} corresponds to 1-year ages, while \code{rc5FDM} corresponds to 5-year age groups.
}
\usage{
data(rc1FDM)
data(rc5FDM)
}
\format{
  A data frame where countries and ages are rows. It has four columns:
  \describe{
  \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{https://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
  \item{\code{age}}{Either single ages from 0 to 100 (\code{rc1FDM}) or 5-year age groups, such as \dQuote{0-4}, \dQuote{5-9}, \dots, \dQuote{100+} (\code{rc5FDM}).}
  }
}
\details{
These datasets are used as the default datasets in \code{\link{pop.predict}} if \code{mig.age.method} is either \dQuote{fdmp} or \dQuote{fdmnop} and the \code{inputs} item \dQuote{mig.fdm} is not given. Other default parameters of the FDM method are read from the \code{\link{vwBaseYear}} dataset.
}

\source{
Most of the values were provided by the United Nations Population Division.
}

\references{
H. Sevcikova, J. Raymer J., A. E. Raftery (2024). Forecasting Net Migration By Age: The Flow-Difference Approach. \href{https://arxiv.org/abs/2411.09878}{arXiv:2411.09878}.
}


\seealso{
\code{\link{vwBaseYear}}
}
\examples{
data(rc1FDM)
head(rc1FDM)
}

\keyword{ datasets }
