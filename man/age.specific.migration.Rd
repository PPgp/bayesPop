\name{age.specific.migration}
\alias{age.specific.migration}

\title{
Reconstruction of Sex- and Age-specific Migration
}
\description{
Reconstructs the sex- and age-specific net migration datasets out of the total net migration using a residual method.
}
\usage{
age.specific.migration(wpp.year = 2015, years = seq(1955, 2100, by = 5), 
    countries = NULL, smooth = TRUE, rescale = TRUE, ages.to.zero = 18:21,
    write.to.disk = FALSE, directory = getwd(), file.prefix = "migration", 
    verbose = TRUE)
}

\arguments{
  \item{wpp.year}{Integer determining which \pkg{wpp} package should be used to get the necessary data from. That package is required to have a dataset on total net migration (called \code{\link[wpp2015]{migration}}). Currently, only the \code{\link{wpp2015}} package contains it (see Details).}
  \item{years}{Array of years that the reconstruction should be made for. This should be a subset of years for which the total net migration is available.}
  \item{countries}{Numerical country codes to do the reconstruction for. By default it is performed on all countries included in the \code{\link[wpp2015]{migration}} dataset where aggregations are excluded.}
  \item{smooth}{Logical controlling if smoothing of the reconstructed curves is required. Due to rounding issues the residual method often yields unrealistic zig-zags on migration curves by age. Smoothing usually improves their look.}
  \item{rescale}{Logical controlling if the resulting migration should be rescaled to match the total migration.}
  \item{ages.to.zero}{Indices of age groups where migration should be set to zero. Default is 85 and older.}
  \item{write.to.disk}{If \code{TRUE} results are written to disk.}
  \item{directory}{Directory where to write the results if \code{write.to.disk} is \code{TRUE}.}
  \item{file.prefix}{If \code{write.to.disk} is \code{TRUE} results are written into two text files with this prefix, a letter \dQuote{M} and \dQuote{F} determining the sex, and concluded by the \dQuote{.txt} suffix. By default \dQuote{migrationM.txt} and \dQuote{migrationF.txt}.}
  \item{verbose}{Logical controlling the amount of output messages.}
}
\details{
Unlike in \code{\link[wpp2012]{wpp2012}}, for the newest release of the WPP, the \code{\link{wpp2015}}, the UN Population Division did not publish the sex- and age-specific net migration counts, only the totals. However, since the sex- and age-schedules are needed for population projections, this function attempts to reconstruct those missing datasets. It uses the published population projections by age and sex, fertility and mortality projections from the \pkg{wpp} package. It computes the population projection without migration and sets the residual to the published population projection as the net migration. By default such numbers are then scaled so that the sum over sexes and ages corresponds to the total migration count. 

If \code{smooth} is \code{TRUE} a smoothing procedure is performed over ages where necessary. Also, for simplicity, we set migration of old ages to zero (default is 85+). Both is done before the scaling. If it is desired to obtain raw residuals without any additional processing, set \code{smooth=FALSE}, \code{rescale=FALSE}, \code{ages.to.zero=c()}. 
}
\value{
List of two data frames (\code{male} and \code{female}), each having the same structure as  \code{\link[wpp2012]{migrationM}}.
}

\author{
Hana Sevcikova
}

\note{The function is called automatically from \code{\link{pop.predict}} if no migration inputs is given. Thus, only users that need sex- and age-specific migration for other purposes will need to call this function explicitly.}
 
\section{Warning}{
Due to rounding issues and slight differences in the methodology, this function does not reproduce the unpublished UN datasets exactly. It is only an approximation! Especially, the first age groups might be more off than other ages.
}



\seealso{
\code{\link{pop.predict}}, \code{\link[wpp2015]{migration}} \code{\link[wpp2012]{migrationM}}
}
\examples{
\dontrun{
asmig <- age.specific.migration()
head(asmig$male)
head(asmig$female)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

