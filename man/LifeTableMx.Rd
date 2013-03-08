\name{LifeTableMx}
\alias{LifeTableMx}
\alias{LifeTableMxCol}
\title{
Life Table Functions
}
\description{
Functions for obtaining life table quantities. 
}
\usage{
LifeTableMx(mx, sex = c("Male", "Female"), include01 = FALSE)

LifeTableMxCol(mx, colname, \dots)
}
\arguments{
  \item{mx}{Vector of age-specific mortality rates nmx. The elements correspond to 1m0, 4m1, 5m5, 5m10, \dots. It can have no more than 28 elements which corresponds to age up to 130. In the \code{LifeTableMxCol} function, this argument can be a two-dimensional matrix with first dimension being the age.}
  \item{sex}{For which sex is the life table.}
  \item{include01}{Logical controlling if the resulting life table contains the age group 0-1.}
  \item{colname}{Name of the column of the life table that should be returned.}
}
\details{
Function \code{LifeTableMx} returns a life table for one set of mortality rates. Function \code{LifeTableMxCol} returns one column of the life table for (possibly) multiple sets of mortality rates.
}
\value{
Function \code{LifeTableMx} returns a data frame with the following elements:
\item{age}{Age groups}
\item{Lx}{nLx, person-years lived between ages n and x+n.}
\item{lx}{lx, number left alive at age x.}
\item{qx}{nqx, probability of dying between ages x ad x+n.}
\item{ax}{nax, average person-years lived in the interval by those dying in the interval.}
}
\references{
Preston, P., Heuveline, P., Guillot, M. (2001): Demography. Blackwell Publishing Ltd.
}
\author{
Hana Sevcikova
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
