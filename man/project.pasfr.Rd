\name{project.pasfr}
\alias{project.pasfr}

\title{
Projections of Percent Age-Specific Fertily Rate
}
\description{
The projections of percent age-specific fertility rate (PASFR) is normally computed within the \code{\link{pop.predict}} function for each trajectory. This function allows to project PASFR outside of population projections for the median total fertility and export it.
}
\usage{
project.pasfr(inputs = NULL, present.year = 2015, end.year = 2100, 
    wpp.year = 2015, digits = 2, out.file.name = "percentASFR.txt")
}

\arguments{
  \item{inputs}{List of input data (file names) with the same meaning as in \code{\link{pop.predict}}. The relevant items here are: either \code{tfr.file} or \code{tfr.sim.dir} (TFR projections), \code{pasfr} (PASFR for observed time periods), \code{srb} (observed and predicted sex ratio at birth), \code{patterns} (PASFR patterns).  All entries are optional. By default the data is taken from the corresponding wpp package.}
  \item{present.year}{Year of the last observed data point.}
  \item{end.year}{End year of the projection.}
  \item{wpp.year}{Year for which WPP data is used if one of the \code{inputs} components is left out.}
  \item{digits}{Number of decimal places in the results.}
  \item{out.file.name}{Name of the resulting file. If \code{NULL} nothing is written.}
}
\details{
The median total fertility rate is derived from the inputs and PASFR corresponding to this median is projected using the method from Sevcikova et al (2015).
}
\value{
Returns invisible data frame with the projected PASFR. 
}

\references{
H. Sevcikova, N. Li, V. Kantorova, P. Gerland and A. E. Raftery (2015). Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. arXiv:1503.05215. \url{http://arxiv.org/abs/1503.05215}
}
\author{
Igor Ribeiro, Hana Sevcikova
}


\seealso{
\code{\link{pop.predict}}
}

\examples{
inputs <- list(tfr.sim.dir=file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output"))
pasfr <- project.pasfr(inputs, out.file.name=NULL)
head(pasfr)
}

\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
\keyword{ manip }% __ONLY ONE__ keyword per line