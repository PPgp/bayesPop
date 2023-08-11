\name{age.specific.migration}
\alias{age.specific.migration}
\alias{migration.totals2age}
\alias{rcastro.schedule}

\title{
Generate Sex- and Age-specific Migration
}
\description{
Creates sex- and age-specific net migration datasets out of the total net migration using either a residual method or given age schedules.
}

\usage{
age.specific.migration(wpp.year = 2019, years = seq(1955, 2100, by = 5), 
    countries = NULL, smooth = TRUE, rescale = TRUE, ages.to.zero = 18:21,
    write.to.disk = FALSE, directory = getwd(), file.prefix = "migration", 
    depratio = wpp.year == 2015, verbose = TRUE)
    
migration.totals2age(df, ages = NULL, annual = FALSE, time.periods = NULL, 
    schedule = NULL, scale = 1, method = "auto", sex = "M",
    id.col = "country_code", country_code = NULL, mig.is.rate = FALSE, 
    alt.schedule.file = NULL, \dots)
    
rcastro.schedule(annual = FALSE)
}

\arguments{
  \item{wpp.year}{Integer determining which \pkg{wpp} package should be used to get the necessary data from. That package is required to have a dataset on total net migration (called \code{\link[wpp2019]{migration}}). Currently, only packages \code{wpp2015}, \code{wpp2017} and \code{\link{wpp2019}} contain it (see Details).}
  \item{years}{Array of years that the reconstruction should be made for. This should be a subset of years for which the total net migration is available.}
  \item{countries}{Numerical country codes to do the reconstruction for. By default it is performed on all countries included in the \code{\link[wpp2019]{migration}} dataset where aggregations are excluded.}
  \item{smooth}{Logical controlling if smoothing of the reconstructed curves is required. Due to rounding issues the residual method often yields unrealistic zig-zags on migration curves by age. Smoothing usually improves their look.}
  \item{rescale}{Logical controlling if the resulting migration should be rescaled to match the total migration.}
  \item{ages.to.zero}{Indices of age groups where migration should be set to zero. Default is 85 and older.}
  \item{write.to.disk}{If \code{TRUE} results are written to disk.}
  \item{directory}{Directory where to write the results if \code{write.to.disk} is \code{TRUE}.}
  \item{file.prefix}{If \code{write.to.disk} is \code{TRUE} results are written into two text files with this prefix, a letter \dQuote{M} and \dQuote{F} determining the sex, and concluded by the \dQuote{.txt} suffix. By default \dQuote{migrationM.txt} and \dQuote{migrationF.txt}.}
  \item{depratio}{If it is \code{TRUE} it will use an internal dataset on migration dependency ratios to adjust the first three age groups. It can also be a name of a binary file containing such dataset.}
  \item{verbose}{Logical controlling the amount of output messages.}
  \item{df}{data.frame, marix or data.table containing total migration counts or rates. Columns correspond to time, rows correspond to locations. Column \dQuote{country_code} (or column identified by \code{id.col}) contains identifiers of the locations. Names of the time columns should be either single years if \code{annual} is \code{TRUE}, e.g. \dQuote{2018}, \dQuote{2019} etc.,  or five year time periods if \code{annual} is \code{FALSE}, e.g. \dQuote{2010-2015}, \dQuote{2015-2020} etc.}
  \item{ages}{Labels of age groups into which the total migration is to be disaggregated. If it is missing, default age groups are determined depending on the argument \code{annual}.}
  \item{annual}{Logical determining if the age groups are 5-year age groups (\code{FALSE}) or 1-year ages (\code{TRUE}) on which the choice of the default schedule is dependent, if \code{schedule} is missing. It also determines the expected syntax of the names of time columns in \code{df}.}
  \item{time.periods}{Character vector determining which columns should be considered in the \code{df} dataset. It should be a subset of column names in \code{df}. By default, all time columns in \code{df} are considered.}
  \item{schedule}{Numeric vector used for distributing total migration across ages. It should have the same length as \code{ages}. If not given, the argument \code{method} determines how the totals are distributed.}
  \item{scale}{The migration schedule is multiplied by this number. It can be used for example, if total migration needs to be distributed between sexes.}
  \item{method}{Method to use for the distribution of totals into age groups. The \dQuote{un} method, available only for \pkg{wpp2022} and only for projected migration, uses the UN migration schedules. The \dQuote{rc} method uses a simple Rogers-Castro disaggregation via the function \code{rcastro.schedule}. The \dQuote{user} method allows to provide expernal schedules in a file given in the argument \code{alt.schedule.file}. The \dQuote{auto} method tries to apply the \dQuote{un} method first and for countries or time periods where there is no match, it applies the \dQuote{rc} method.}
  \item{sex}{\dQuote{M} or \dQuote{F} determining the sex of this schedule. It is used for matching to the UN schedules, i.e. if the \dQuote{auto} or \dQuote{un} method is used.}
  \item{id.col}{Name of the unique identifier of the locations.}
  \item{country_code}{Code of the country if the data belongs to only country only. It is used for matching to the UN schedules, i.e. if the \dQuote{auto} or \dQuote{un} method is used.}
  \item{mig.is.rate}{Logical indicating if the data in \code{df} should be interpreted as rates. If \code{FALSE}, \code{df} represent counts.}
  \item{alt.schedule.file}{If \code{method} is \dQuote{user}, alternative migration age schedules can be provided in this binary file. It has the same structure as the internal package objects \code{mig1.schedule} (annual) or \code{mig5.schedule} (5-year). If there are different schedules for positive and negative totals, the negative schedules are in \code{mig1.neg.schedule} or \code{mig5.neg.schedule} objects.}
  \item{\dots}{Further arguments passed to the underlying functions.}
}
\details{
\subsection{Function \code{age.specific.migration}}{
Unlike in \code{\link[wpp2012]{wpp2012}}, for the three latest releases of the WPP, the \pkg{wpp2015}, \pkg{wpp2017}, and \pkg{\link{wpp2019}}, the UN Population Division did not publish the sex- and age-specific net migration counts, only the totals. However, since the sex- and age-schedules are needed for population projections, the \code{age.specific.migration} function attempts to reconstruct those missing datasets. It uses the published population projections by age and sex, fertility and mortality projections from the \pkg{wpp} package. It computes the population projection without migration and sets the residual to the published population projection as the net migration. By default such numbers are then scaled so that the sum over sexes and ages corresponds to the total migration count. 

If \code{smooth} is \code{TRUE} a smoothing procedure is performed over ages where necessary. Also, for simplicity, we set migration of old ages to zero (default is 85+). Both is done before the scaling. If it is desired to obtain raw residuals without any additional processing, set \code{smooth=FALSE}, \code{rescale=FALSE}, \code{ages.to.zero=c()}. 
}
\subsection{Function \code{migration.totals2age}}{
This function should be used when working with \pkg{wpp2022} data. It allows users to disagregate total migration counts or rates (for multiple time periods and multiple locations) into age-specific ones by either the UN schedule used in WPP2022 (\code{method = "un"}), a Rogers-Castro (\code{method = "rc"}), a user-defined schedule (\code{method = "user"}) or a combination of the first two (\code{method = "auto"}). In the latter case, the UN schedule is used for records that match in terms of country code and time periods. For the remaining records, a Rogers-Castro schedule is used. The latter schedule can be accessed via the function \code{rcastro.schedule} where the \code{annual} argument specifies if it is for 1-year or 5-year age groups. 
}
}

\value{
Function \code{age.specific.migration} returns a list of two data frames (\code{male} and \code{female}), each having the same structure as  \code{\link[wpp2012]{migrationM}}.

Function \code{migration.totals2age} returns a data.table with the disaggregated counts. 

Function \code{rcastro.schedule} returns a vector of proportions for each age group.
}

\author{
Hana Sevcikova
}

\note{The \code{age.specific.migration} function is called automatically from \code{\link{pop.predict}} if no migration input is given. Similarly, if migration is given as totals, i.e. via the entries \code{mig}, \code{migMt}, or \code{migFt}, the function \code{migration.totals2age} is called using either the UN schedules or the default Rogers-Castro schedule. 
Thus, only users that need sex- and age-specific migration for other purposes, or modify the defaults, will need to call these functions explicitly.}
 
\section{Warning}{
Due to rounding issues and slight differences in the methodology, the \code{age.specific.migration} function does not reproduce the unpublished UN datasets exactly. It is only an approximation! Especially, the first age groups might be more off than other ages.
}



\seealso{
\code{\link{pop.predict}}, \code{\link[wpp2019]{migration}} \code{\link[wpp2012]{migrationM}}
}
\examples{
\dontrun{
asmig <- age.specific.migration()
head(asmig$male)
head(asmig$female)}

# simple disaggregation for one location
totmig <- c(30, -50, -100)
names(totmig) <- 2018:2020
asmig.simple <- migration.totals2age(totmig, annual = TRUE, method = "rc")
head(asmig.simple)

\dontrun{
# disaggregate WPP 2019 migration for all countries, one sex
data(migration, package = "wpp2019")
# assuming equal sex migration ratio
asmig.all <- migration.totals2age(migration, scale = 0.5, method = "rc") 
# plot result for the US in 2095-2100
mig1sex.us <- subset(asmig.all, country_code == 840)[["2095-2100"]]
plot(ts(mig1sex.us))
# check that the sum is half of the original total
sum(mig1sex.us) == subset(migration, country_code == 840)[["2095-2100"]]/2}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ datagen }
