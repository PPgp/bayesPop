\name{age.specific.migration}
\alias{age.specific.migration}
\alias{migration.totals2age}
\alias{rcastro.schedule}

\title{
Generate Sex- and Age-specific Migration
}
\description{
Creates sex- and age-specific net migration datasets out of the total net migration using different methods. The \code{age.specific.migration} is a legacy function that distributes UN 5-year totals into ages using a residual method. The \code{migration.totals2age} distribute given totals using Rogers-Castro and the Flow Difference Method (FDM).
}

\usage{
age.specific.migration(wpp.year = 2019, years = seq(1955, 2100, by = 5), 
    countries = NULL, smooth = TRUE, rescale = TRUE, ages.to.zero = 18:21,
    write.to.disk = FALSE, directory = getwd(), file.prefix = "migration", 
    depratio = wpp.year == 2015, verbose = TRUE)
    
migration.totals2age(df, ages = NULL, annual = FALSE, time.periods = NULL, 
    scale = 1, method = "rc", sex = "M",
    id.col = "country_code", mig.is.rate = FALSE, 
    rc.data = NULL, pop = NULL, pop.glob = NULL, \dots)
    
rcastro.schedule(annual = FALSE)
}

\arguments{
  \item{wpp.year}{Integer determining which \pkg{wpp} package should be used to get the necessary data from. That package is required to have a dataset on total net migration (called \code{\link[wpp2019]{migration}}).}
  \item{years}{Array of years that the reconstruction should be made for. This should be a subset of years for which the total net migration is available.}
  \item{countries}{Numerical country codes to do the reconstruction for. By default it is performed on all countries included in the \code{\link[wpp2019]{migration}} dataset where aggregations are excluded.}
  \item{smooth}{Logical controlling if smoothing of the reconstructed curves is required. Due to rounding issues the residual method often yields unrealistic zig-zags on migration curves by age. Smoothing usually improves their look.}
  \item{rescale}{Logical controlling if the resulting migration should be rescaled to match the total migration.}
  \item{ages.to.zero}{Indices of age groups where migration should be set to zero. Default is 85 and older.}
  \item{write.to.disk}{If \code{TRUE} results are written to disk.}
  \item{directory}{Directory where to write the results if \code{write.to.disk} is \code{TRUE}.}
  \item{file.prefix}{If \code{write.to.disk} is \code{TRUE} results are written into two text files with this prefix, a letter \dQuote{M} and \dQuote{F} determining the sex, and concluded by the \dQuote{.txt} suffix. By default \dQuote{migrationM.txt} and \dQuote{migrationF.txt}.}
  \item{depratio}{If it is \code{TRUE} it will use an internal dataset on migration dependency ratios to adjust the first three age groups. It can also be a name of a binary file containing such dataset. The default dataset is only available for 2015.}
  \item{verbose}{Logical controlling the amount of output messages.}
  \item{df}{data.frame, marix or data.table containing total migration counts or rates. Columns correspond to time, rows correspond to locations. Column \dQuote{country_code} (or column identified by \code{id.col}) contains identifiers of the locations. Names of the time columns should be either single years if \code{annual} is \code{TRUE}, e.g. \dQuote{2018}, \dQuote{2019} etc.,  or five year time periods if \code{annual} is \code{FALSE}, e.g. \dQuote{2010-2015}, \dQuote{2015-2020} etc.}
  \item{ages}{Labels of age groups into which the total migration is to be disaggregated. If it is missing, default age groups are determined depending on the argument \code{annual}.}
  \item{annual}{Logical determining if the age groups are 5-year age groups (\code{FALSE}) or 1-year ages (\code{TRUE}) on which the choice of the default schedule is dependent, if \code{schedule} is missing. It also determines the expected syntax of the names of time columns in \code{df}.}
  \item{time.periods}{Character vector determining which columns should be considered in the \code{df} dataset. It should be a subset of column names in \code{df}. By default, all time columns in \code{df} are considered.}
  \item{scale}{The migration schedule is multiplied by this number. It can be used for example, if total migration needs to be distributed between sexes.}
  \item{method}{Method to use for the distribution of totals into age groups. The \dQuote{rc} method uses either a basic Rogers-Castro disaggregation via the function \code{rcastro.schedule}, or a schedule given in the \code{rc.data} argument. The \dQuote{fdmp} and \dQuote{fdmnop} methods use the Flow Difference Method, where \dQuote{fdmp} weights the flows by population.}
  \item{sex}{\dQuote{M} or \dQuote{F} determining the sex of this schedule. It only impacts the FDM methods.}
  \item{id.col}{Name of the unique identifier of the locations.}
  \item{mig.is.rate}{Logical indicating if the data in \code{df} should be interpreted as rates. If \code{FALSE}, \code{df} represent counts.}
  \item{rc.data}{data.table containing either a family of Rogers-Castro proportions if \code{method = "rc"}, or various inputs for the FDM methods if \code{method} is either \dQuote{fdmp} or \dQuote{fdmnop}. 
  
  For the \dQuote{rc} method, mandatory columns are \dQuote{age} and \dQuote{prop}. Optionally, it can have a column \dQuote{mig_sign} with values \dQuote{Inmigration} and \dQuote{Emigration} (distinguishing schedules to be applied for positive and negative migration, respectively) and a column \dQuote{sex} with values \dQuote{Female} and \dQuote{Male}. The format corresponds to the dataset \code{DemoTools::mig_un_families}, subset to a single family. 
  
  For the FDM methods, it has columns contained in the \code{\link{rcFDM}} dataset, as well as columns \dQuote{beta0} (intercept), \dQuote{beta1} (slope), \dQuote{min} (minimum rate), \dQuote{in_sex_factor} (inflow female proportion), and \dQuote{out_sex_factor} (outflow female proportion), used in the FDM methods. These columns correspond to columns \dQuote{MigFDMb0}, \dQuote{MigFDMb1}, \dQuote{MigFDMmin}, \dQuote{MigFDMsrin} and \dQuote{MigFDMsrout}, respectively, in the \code{\link{vwBaseYear}} dataset.}
  \item{pop}{data.table with population counts needed for the FDM methods. It should have a location identifier column of the same name as \code{id.col}, further columns \dQuote{year}, \dQuote{age}, and \dQuote{pop}.}
  \item{pop.glob}{data.table with global population needed for the weighted FDM method (\dQuote{fdmp}). It should have columns \dQuote{year}, \dQuote{age}, and \dQuote{pop}.}
  \item{\dots}{Further arguments passed to the underlying functions.}
}
\details{
\subsection{Function \code{age.specific.migration}}{
Unlike in \code{\link[wpp2012]{wpp2012}}, for the four releases of the WPP between 2015 and 2022, the \pkg{wpp2015}, \pkg{wpp2017}, \pkg{\link[wpp2019]{wpp2019}}, and \href{https://github.com/PPgp/wpp2022}{\pkg{wpp2022}}, the UN Population Division did not publish the sex- and age-specific net migration counts, only the totals. However, since the sex- and age-schedules are needed for population projections, the \code{age.specific.migration} function attempts to reconstruct those missing datasets. It uses the published population projections by age and sex, fertility and mortality projections from the \pkg{wpp} package. It computes the population projection without migration and sets the residual to the published population projection as the net migration. By default such numbers are then scaled so that the sum over sexes and ages corresponds to the total migration count.

If \code{smooth} is \code{TRUE} a smoothing procedure is performed over ages where necessary. Also, for simplicity, we set migration of old ages to zero (default is 85+). Both is done before the scaling. If it is desired to obtain raw residuals without any additional processing, set \code{smooth=FALSE}, \code{rescale=FALSE}, \code{ages.to.zero=c()}. 

This function works only for 5-year data. 
}
\subsection{Function \code{migration.totals2age}}{
This function should be used when working with annual data or data from \href{https://github.com/PPgp/wpp2022}{\pkg{wpp2022}} and \href{https://github.com/PPgp/wpp2024}{\pkg{wpp2024}}. It allows users to disagregate total migration counts or rates (for multiple time periods and multiple locations) into age-specific ones by either a schedule similar to the one used by the UN in WPP2024 (\code{method = "fdmnop"}), a Rogers-Castro (\code{method = "rc"}), or by FDM weighted by population (\code{method = "fdmp"}) as described in Sevcikova et al (2024). The FDM method needs additional info passed via the arguments \code{rc.data}, \code{pop} and \code{pop.glob}. The default Rogers-Castro schedule can be accessed via the function \code{rcastro.schedule} where the \code{annual} argument specifies if it is for 1-year or 5-year age groups. Alternatively, an external schedule can be given via the \code{rc.data} argument, where one can distinguish between schedules for each sex, as well as for positive and negative net migration. It has the same structure as the dataset \code{DemoTools::mig_un_families}, but it should be a subset for a single family and converted to \code{data.table}.
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

\note{These functions are called automatically from \code{\link{pop.predict}} if needed, depending on the inputs. 
Thus, only users that need sex- and age-specific migration for other purposes, or modify the defaults, will need to call these functions explicitly. 

Further note that the \href{https://github.com/PPgp/wpp2024}{\pkg{wpp2024}} package does contain the age-specific net migration for projected years (datasets \code{migprojAge1dt}, \code{migprojAge5dt}). Thus, if running \code{\link{pop.predict}} with \code{wpp.year = 2024} and the default migration totals, no disagregation is necessary for the projected time periods. The disaggregation is only triggerered for the past time periods, or in a case when user-specific net migration totals are used.}
 
\section{Warning}{
Due to rounding issues and slight differences in the methodology, the functions do not reproduce the unpublished UN datasets exactly. It is only an approximation! Especially, the first age groups might be more off than other ages.
}

\references{
H. Sevcikova, J. Raymer J., A. E. Raftery (2024). Forecasting Net Migration By Age: The Flow-Difference Approach. \href{https://arxiv.org/abs/2411.09878}{arXiv:2411.09878}.
}


\seealso{
\code{\link{pop.predict}}, \code{\link[wpp2019]{migration}} \code{\link[wpp2012]{migrationM}}, \code{\link{rcFDM}}, \code{\link{vwBaseYear}}
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
