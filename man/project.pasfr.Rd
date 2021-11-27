\name{project.pasfr}
\alias{project.pasfr}
\alias{project.pasfr.traj}

\title{
Projections of Percent Age-Specific Fertily Rate
}
\description{
The projections of percent age-specific fertility rate (PASFR) is normally computed within the \code{\link{pop.predict}} function for each trajectory. This function allows to project PASFR outside of population projections for the median total fertility rate (TFR) or user-provided TFR, and export it.
}
\usage{
project.pasfr(inputs = NULL, present.year = 2020, end.year = 2100, 
    wpp.year = 2019, annual = FALSE, nr.est.points = if(annual) 15 else 3,
    digits = 2, out.file.name = "percentASFR.txt", verbose = FALSE)
    
project.pasfr.traj(inputs = NULL, countries = NULL, nr.traj = NULL, 
    present.year = 2020, end.year = 2100, wpp.year = 2019, 
    annual = FALSE, nr.est.points = if(annual) 15 else 3,
    digits = 2, out.file.name = "percentASFRtraj.txt", verbose = FALSE)
}

\arguments{
  \item{inputs}{List of input data (file names) with the same meaning as in \code{\link{pop.predict}}. The relevant items here are: either \code{tfr.file} or \code{tfr.sim.dir} (TFR estimates and projections), \code{pasfr} (PASFR for observed time periods), and \code{patterns} (PASFR patterns). All entries are optional. By default the data is taken from the corresponding wpp package. See Details below.}
  \item{present.year}{Year of the last observed data point.}
  \item{end.year}{End year of the projection.}
  \item{wpp.year}{Year for which WPP data is used if one of the \code{inputs} components is left out.}
  \item{annual}{Logical that should be \code{TRUE} if the provided data on TFR and PASFR are annual-based data.}
  \item{nr.est.points}{Number of time points to be used for estimating the continuation of the observed PASFR trend. By default it is 15 years, corresponding to three time points for 5-year data.}
  \item{digits}{Number of decimal places in the results.}
  \item{out.file.name}{Name of the resulting file. If \code{NULL} nothing is written.}
  \item{verbose}{Logical switching verbose messages on and off.}
  \item{countries}{Vector of numerical country codes. By default the function is applied to all countries.}
  \item{nr.traj}{Number of trajectories on which the function should be applied. By default all trajectories are taken. Otherwise they are thinned appropriately.}
}
\details{
If the input TFR is given as an ASCII file (in \code{inputs$tfr.file}), it can be either a csv (comma-separated) file in long format, with columns \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory} and \dQuote{TF}. Or, it can be a tab-separated (wide format) file with column \dQuote{country_code} and each year or time period as a separate column (see \code{\link[wpp2019]{tfr}}). In the latter case, an additional \code{inputs} entry \code{tfr.file.type = "w"} must be provided to specify the file is in the wide format, which is a case whe there is only one trajectory. Note that the TFR input should cover all projection time period as well as observed TFR as the function assesses the start of Phase III, which could be in the past.

If observed PASFR is given (in \code{inputs$pasfr}), it is a tab-separated file in wide format as in \code{\link[wpp2019]{percentASFR}}. Fertility age patterns can be controlled by country via the \code{inputs$patterns} entry, which is a dataset in the same format and meaning as \code{\link{vwBaseYear}}. 

In addition, if the present year differs by country, the \code{inputs} list accepts the entry \code{last.observed}, which is a tab-separated file with columns \dQuote{country_code} and \dQuote{last.observed}. It can contain the year of the last observed time period for each country. 

In the \code{project.pasfr} function, if the TFR input (given either as a long file or as a simulation directory), contains more than one trajectory, the median is derived over the trajectories for each time period. Then, PASFR corresponding to this median is projected using the method from Sevcikova et al (2015).

For \code{project.pasfr.traj}, the PASFR is projected for single trajectories of TFR.
}

\value{
Returns invisible data frame with the projected PASFR. 
}

\references{
H. Sevcikova, N. Li, V. Kantorova, P. Gerland and A. E. Raftery (2016). Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. In: Dynamic Demographic Analysis, ed. Schoen R. (Springer), pp. 285-310. Earlier version: \url{https://arxiv.org/abs/1503.05215}
}

\author{
Hana Sevcikova, Igor Ribeiro
}


\seealso{
\code{\link{pop.predict}}
}

\examples{
# using TFR in simulation directory
inputs <- list(tfr.sim.dir=file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output"))
pasfr <- project.pasfr(inputs, out.file.name = NULL)
head(pasfr)

\dontrun{
pasfr.traj <- project.pasfr.traj(inputs, out.file.name = NULL)
head(pasfr.traj)}

# using TFR in wide-format file
inputs2 <- list(tfr.file = file.path(find.package("wpp2019"), "data", "tfrprojMed.txt"),
    tfr.file.type = "w")
pasfr2 <- project.pasfr(inputs2, out.file.name = NULL)
head(pasfr2)
}

\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
\keyword{ manip }% __ONLY ONE__ keyword per line
