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
LifeTableMx(mx, sex = c("Male", "Female"), include01 = TRUE,
	radix = 1, open.age = 130)

LifeTableMxCol(mx, colname = c("Lx", "lx", "qx", "mx", "dx", "Tx", "sx", "ex", "ax"), \dots)
}
\arguments{
  \item{mx}{Vector of age-specific mortality rates nmx. The elements correspond to 1m0, 4m1, 5m5, 5m10, \dots. It can have no more than 28 elements which corresponds to age up to 130. In the \code{LifeTableMxCol} function, this argument can be a two-dimensional matrix with first dimension being the age.}
  \item{sex}{For which sex is the life table.}
  \item{include01}{Logical. If it is \code{FALSE} the first two age groups (0-1 and 1-4) are collapsed to one age group (0-4).}
  \item{radix}{Base of the life table.}
  \item{open.age}{Open age group. If smaller than the last age group of \code{mxm}, the life table is truncated.}
  \item{colname}{Name of the column of the life table that should be returned.}
  \item{\dots}{Arguments passed to underlying functions. Argument \code{age05} is a logical vector of size three, specifying if the age groups 0-1, 1-4 and 0-5 should be included. Default value of \code{c(FALSE, FALSE, TRUE)} includes the 0-5 age group only.}
}
\details{
Function \code{LifeTableMx} returns a life table for one set of mortality rates. Function \code{LifeTableMxCol} returns one column of the life table for (possibly) multiple sets of mortality rates.
}
\value{
Function \code{LifeTableMx} returns a data frame with the following elements:
\item{age}{Age groups}
\item{mx}{mx, the input vector of mortality rates.}
\item{qx}{nqx, probability of dying between ages x ad x+n.}
\item{lx}{lx, number left alive at age x.}
\item{dx}{ndx, cohort deaths between ages x ad x+n.}
\item{Lx}{nLx, person-years lived between ages x and x+n.}
\item{sx}{sx, survival rate at age x.}
\item{Tx}{Tx, person-years lived above age x.}
\item{ex}{e0x, expectation of life at age x.}
\item{ax}{nax, average person-years lived in the interval by those dying in the interval.}

Function \code{LifeTableMxCol} returns one given column of the life table, possibly as a matrix (if \code{mx} is a matrix).
}

\references{
Preston, P., Heuveline, P., Guillot, M. (2001): Demography. Blackwell Publishing Ltd.
}
\author{
Hana Sevcikova, Thomas Buettner, Nan Li, Patrick Gerland
}

\seealso{
\code{\link{pop.expressions}} for examples on retrieving some life table quantities.
}
\examples{
\dontrun{
sim.dir <- tempfile()
pred <- pop.predict(countries="Ecuador", output.dir=sim.dir, wpp.year=2015,
    present.year=2015, keep.vital.events=TRUE, fixed.mx=TRUE, fixed.pasfr=TRUE)
# get male mortality rates from 2020 for age groups 0-1, 1-4, 5-9, ...
mxm <- pop.byage.table(pred, expression="MEC_M{age.index01(27)}", year=2020)[,1]
print(LifeTableMx(mxm), digits=3)
# female LT with first two age categories collapsed 
mxf <- pop.byage.table(pred, expression="MEC_F{age.index01(27)}", year=2020)[,1]
print(LifeTableMx(mxf, sex="Female", include01=FALSE), digits=3)
unlink(sim.dir, recursive=TRUE)}
}

\keyword{ attribute }
