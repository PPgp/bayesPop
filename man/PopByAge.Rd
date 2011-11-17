\name{PopByAgeMale2008}
\alias{PopByAge}
\alias{PopByAgeMale}
\alias{PopByAgeFemale}
\alias{PopByAgeMale2008}
\alias{PopByAgeMale2010}
\alias{PopByAgeFemale2008}
\alias{PopByAgeFemale2010}

\docType{data}
\title{
Datasets of Initial Population
}
\description{
Male and female datasets of initial population (i.e. at present time). Each dataset is connected to a WPP year which is part of the dataset name.
}
\usage{
data(PopByAgeMale2008)
data(PopByAgeMale2010)

data(PopByAgeFemale2008)
data(PopByAgeFemale2010)
}
\format{
  Data frames containing the following 9 variables:
  \describe{
    \item{\code{index}}{Numerical index of the records. Not used by the package.}
    \item{\code{country}}{Country name. Not used by the package.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{age}}{A character string representing an age interval. For each country there are 21 values: \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dQuote{15-19}, \dQuote{20-24}, \dQuote{25-29}, \dQuote{30-34}, \dQuote{35-39}, \dQuote{40-44}, \dQuote{45-49},  \dQuote{50-54}, \dQuote{55-59}, \dQuote{60-64}, \dQuote{65-69}, \dQuote{70-74}, \dQuote{75-79}, \dQuote{80-84}, \dQuote{85-89}, \dQuote{90-94}, \dQuote{95-99}, and \dQuote{100+} in that order.}
    \item{\code{2005}, \code{2010}}{Value of the initial population. The header of this column must correspond to the present year.}
  }
}
\details{
There is one record for each country and each age. Thus, for each country there are 21 records. The country records must be groupped together with age records following the order above. 
}
\source{
Data provided by the United Nations Population Division
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(PopByAgeMale2008)
}
\keyword{datasets}
