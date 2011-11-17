\name{MigByAgeMale2008}
\alias{MigByAge}
\alias{MigByAgeMale}
\alias{MigByAgeMale2008}
\alias{MigByAgeMale2010}
\alias{MigByAgeFemale}
\alias{MigByAgeFemale2008}
\alias{MigByAgeFemale2010}
\docType{data}
\title{
Datasets on Migration
}
\description{
Datasets giving prediction of male and female age-specific net migration.
}
\usage{
data(MigByAgeMale2008)
data(MigByAgeMale2010)

data(MigByAgeFemale2008)
data(MigByAgeFemale2010)
}
\format{
  Data frames containing the following variables.
  \describe{
    \item{\code{index}}{Numerical index of the records. Not used by the package.}
    \item{\code{country}}{Country name. Not used by the package.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{age}}{A character string representing an age interval. For each country there are 21 values: \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dQuote{15-19}, \dQuote{20-24}, \dQuote{25-29}, \dQuote{30-34}, \dQuote{35-39}, \dQuote{40-44}, \dQuote{45-49},  \dQuote{50-54}, \dQuote{55-59}, \dQuote{60-64}, \dQuote{65-69}, \dQuote{70-74}, \dQuote{75-79}, \dQuote{80-84}, \dQuote{85-89}, \dQuote{90-94}, \dQuote{95-99}, and \dQuote{100+} in that order.}
    \item{\code{2005-2010}, \code{2010-2015}, \code{2015-2020}, \dots}{Net migration for the specific time period. It must contain one column per each projection period, i.e. between \code{present.year} and \code{end.year} of the \code{\link{pop.predict}} function.}
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
data(MigByAgeMale2010)
str(MigByAgeMale2010)
}
\keyword{datasets}
