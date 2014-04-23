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
LifeTableMx(mx, sex = c("Male", "Female"))

LifeTableMxCol(mx, colname=c('Lx', 'lx', 'qx', 'mx'), \dots)
}
\arguments{
  \item{mx}{Vector of age-specific mortality rates nmx. The elements correspond to 1m0, 4m1, 5m5, 5m10, \dots. It can have no more than 28 elements which corresponds to age up to 130. In the \code{LifeTableMxCol} function, this argument can be a two-dimensional matrix with first dimension being the age.}
  \item{sex}{For which sex is the life table.}
  \item{colname}{Name of the column of the life table that should be returned.}
  \item{\dots}{Arguments passed to underlying functions. Argument \code{age05} is a logical vector of size three, specifying if the age groups 0-1, 1-4 and 0-5 should be included. Default value of \code{c(FALSE, FALSE, TRUE)} includes the 0-5 age group only.}
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

Function \code{LifeTableMxCol} returns one given column of the life table, possibly as a matrix (if \code{mx} is a matrix).

If \code{include01} is \code{TRUE} the number of rows on the data frame corresponds to the size of \code{mx}. In such a case, the first row corresponds to the age group 0-1 and the second row to the age group 0-5 (not 1-5!). If \code{include01} is \code{FALSE}, the first row is eliminated.
}

\references{
Preston, P., Heuveline, P., Guillot, M. (2001): Demography. Blackwell Publishing Ltd.
}
\author{
Hana Sevcikova, Nan Li
}

\seealso{
\code{\link{pop.expressions}} for examples on retrieving some life table quantities.
}
\examples{
sim.dir <- tempfile()
pred <- pop.predict(countries="Ecuador", output.dir=sim.dir,
    keep.vital.events=TRUE)
# get male mortality rates from 2000 for age groups 0-1, 1-4, 5-9, ...
mx <- pop.byage.table(pred, expression="MEC_M{c(-1,0,2:21)}")[,1]
LifeTableMx(mx)
unlink(sim.dir, recursive=TRUE)
}

\keyword{ attribute }
