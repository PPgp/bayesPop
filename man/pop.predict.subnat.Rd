\name{pop.predict.subnat}
\alias{pop.predict.subnat}

\title{
Subnational Probabilistic Population Projection
}
\description{
Generates trajectories of probabilistic population projection for subregions of a given country.
}
\usage{
pop.predict.subnat(end.year = 2060, start.year = 1950, present.year = 2015, 
        wpp.year = 2017, output.dir = file.path(getwd(), "bayesPop.output"), 
        locations = NULL, default.country = NULL, 
        inputs = list(
            popM = NULL, popF = NULL, 
            mxM = NULL, mxF = NULL, srb = NULL, 
            pasfr = NULL, patterns = NULL, 
            migM = NULL, migF = NULL, 
            e0F.file = NULL, e0M.file = NULL, tfr.file = NULL, 
            e0F.sim.dir = NULL, e0M.sim.dir = NULL, tfr.sim.dir = NULL, 
            migMtraj = NULL, migFtraj = NULL
        ), 
        nr.traj = 1000, keep.vital.events = FALSE, 
        fixed.mx = FALSE, fixed.pasfr = FALSE, 
        replace.output = FALSE, verbose = TRUE)
}

\arguments{
  \item{end.year}{End year of the projection.}
  \item{start.year}{First year of the historical data on mortality rates. It determines the length of the historical time series used in the Lee-Carter estimation.}
  \item{present.year}{Year for which initial population data is to be used.}
  \item{wpp.year}{Year for which WPP data is used. The function loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and uses its data (corresponding to the \code{default.country}) as default datasets if region-specific alternatives are not given (see more details below).}
  \item{output.dir}{Output directory of the projection.}
  \item{locations}{Name of a tab-delimited file that contains definitions of the subregions. It has a similar structure as \code{\link[wpp2017]{UNlocations}}, with mandatory columns \code{reg_code} (unique identifier of the subregions) and \code{name} (name of the subregions). Optionally, \code{location_type} should be set to 4 for subregions to be processed. Column \code{country_code} can be included with the numerical code of the corresponding country. A row with \code{location_type} of 0 determines the country that the subregions belong to and is used for extracting default "national" datasets if the argument \code{default.country} is missing. In such a case, the code of the default country is taken from its column \code{country_code}. This is a mandatory argument. }
  \item{default.country}{Numerical code of a country to which the subregions belong to. It is used for extracting default datasets from the  \pkg{wpp} package if some region-specific input datasets are missing. Alternatively, it can be also included in the \code{locations} file, see above. In either case, the code must exists in the \code{\link[wpp2017]{UNlocations}} dataset.}
  \item{inputs}{A list of file names where input data is stored. Unless otherwise noted, these are tab delimited ASCII files with a mandatory column \code{reg_code} giving the numerical identifier of the subregions. If an element of this list is \code{NULL}, usually a default dataset corresponding to \code{default.country} is extracted from the \pkg{wpp} package. Names of these default datasets are shown in brackets. This list contains the following elements:
  		\describe{
  			\item{popM, popF}{Initial male/female age-specific population (at time \code{present.year}). Mandatory items, no defaults. Must contain columns \code{reg_code} and \code{age} and be of the same structure as \code{\link[wpp2017]{popM}} from \pkg{wpp}.}
  			\item{mxM, mxF}{Historical data and (optionally) projections of male/female age-specific death rates [\code{\link[wpp2017]{mxM}}, \code{\link[wpp2017]{mxF}}] (see also argument \code{fixed.mx}).}
  			\item{srb}{Projection of sex ratio at birth. [\code{\link[wpp2017]{sexRatio}}]}
  			\item{pasfr}{Historical data and (optionally) projections of percentage age-specific fertility rate [\code{\link[wpp2017]{percentASFR}}] (see also argument \code{fixed.pasfr}).}
  			\item{patterns}{Information on region's specifics regarding migration type, base year of the migration, mortality and fertility age patterns as defined in [\code{\link{vwBaseYear}}]. In addition, it can contain columns defining migration shares between the subregions, see Details below.}
  			\item{migM, migF}{Projection of male/female age-specific migration as net counts on the same scale as initital population. It should have the same format as \code{\link[wpp2012]{migrationM}}. If not available, the migration schedules are constructed from total migration counts of the \code{default.country} derived from \code{\link[wpp2017]{migration}} using Rogers Castro for age distribution. Migration shares between subregions (including sex-specific shares) can be given in the \code{patterns} file, see above and Details below.}
  			\item{e0F.file}{Comma-delimited CSV file with projected female life expectancy. It has the same structure as the file \dQuote{ascii_trajectories.csv} generated using \code{bayesLife::\link[bayesLife]{convert.e0.trajectories}} (which currently works for country-level results only). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{e0}. If \code{e0F.file} is \code{NULL}, data from the corresponding \pkg{wpp} package (for \code{default.country}) is taken, namely the median projections as one trajectory and the low and high variants (if available) as second and third trajectory. Alternatively, this element can be the keyword \dQuote{median_} in which case only the median is taken.}
  			\item{e0M.file}{Comma-delimited CSV file containing projections of male life expectancy of the same format as \code{e0F.file}. As in the female case, if \code{e0M.file} is \code{NULL}, data for \code{default.country} from the corresponding \pkg{wpp} package is taken.}
  			\item{tfr.file}{Comma-delimited CSV file with results of total fertility rate (generated using \pkg{\link{bayesTFR}}, function \code{\link{convert.tfr.trajectories}}, file \dQuote{ascii_trajectories.csv}). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{TF}. If this element is not \code{NULL}, the argument \code{tfr.sim.dir} is ignored. If both \code{tfr.file} and \code{tfr.sim.dir} are \code{NULL}, data for \code{default.country} from the corresponding \pkg{wpp} package is taken (median and the low and high variants as three trajectories). Alternatively, this argument can be the keyword \dQuote{median_} in which case only the wpp median is taken.}
  			\item{e0F.sim.dir}{Simulation directory with results of female life expectancy. Since \pkg{\link{bayesLife}} does not support subnational projections yet, this element should not be used. Instead use \code{e0F.file} if region-specific e0 projections are avaialable. Alternatively, it can be set to the keyword \dQuote{median_} which has the same effect as when \code{e0F.file} is \dQuote{median_}.}
  			\item{e0M.sim.dir}{This is analogous to \code{e0F.sim.dir}, here for male life expectancy. Use \code{e0M.file} instead of this item.}
  			\item{tfr.sim.dir}{Simulation directory with projections of total fertility rate (generated using \code{bayesTFR::\link[bayesTFR]{tfr.predict.subnat}}). It is only used if \code{tfr.file} is \code{NULL}.}
  			\item{migMtraj, migFtraj}{Comma-delimited CSV file with male/female age-specific migration trajectories. If present, it replaces deterministic projections given by the \code{migM} and \code{migF} items. It has a similar format as e.g. \code{e0M.file} with columns \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, \dQuote{Age} and \dQuote{Migration}. The \dQuote{Age} column must have values \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dots, \dQuote{95-99}, \dQuote{100+}.}
  		}
}
  \item{nr.traj, keep.vital.events, fixed.mx, fixed.pasfr, replace.output, verbose}{These arguments have the same meaning as in \code{\link{pop.predict}}.}
}

\details{
Population projection for subnational units (regions) is performed by applying the cohort component method to subnational datasets on projected fertility (TFR), mortality and net migration, starting from given sex- and age-specific population counts. The only required inputs are the initial sex- and age-specific population counts in each region (\code{popM} and \code{popF} elements of the \code{inputs} argument) and a file with a set of locations (argument \code{locations}). If no other input datasets are given, those datasets are replaced by the corresponding "national" values, taken from the corresponding \pkg{wpp} package. The argument \code{default.country} determines the country for those default "national" values. The default country can be also included in the \code{locations} file as a record with \code{location.type} being set to 0.

The TFR component can be given as a set of trajectories generated using the \code{\link[bayesTFR]{tfr.predict.subnat}} function of the \pkg{bayesTFR} package (\code{tfr.sim.dir} element). Alternatively, trajectories can be given in an ASCII file (\code{tfr.file}). Having a set of subnational TFR trajectories, the cohort component method is applied to each of them to yield a distribution of future subnational population.

Net migration can either be given as disaggregated sex- and age-specific datasets \code{migM} and \code{migF}. Alternatively, it can be given as shares between regions as columns in the \code{patterns} dataset. These are: \code{inmigrationM_share}, \code{inmigrationF_share}, \code{outmigrationM_share}, \code{outmigrationF_share}. The sex specification and/or direction specification (in/out) can be omitted, e.g. it can be simply \code{migration_share}. The function extracts the values of net migration projection on the national level and distributes it to regions according to the given shares. For positive (national) values, it uses the in-migration shares; for negative values it uses the out-migration shares. If the in/out prefix is omitted in the column names, the given migartion shares are used for both, positive and negative net migration projection. By default, if no \code{migM} and \code{migF} neither region-specific shares are given, the distribution between regions is proportional to the size of population. The age-specific schedules follow by default the Rogers-Castro age schedules. Note that when handling migration using shares as described here, it only affects the distribution of international migration into regions. It does not take into account between-region migration.

The package contains example datasets for Canada. Use these as templates for your own data. See Example below.
}

\value{
Object of class \code{\link{bayesPop.prediction}} containing the subnational projections. Note that this object can be use in the various \pkg{bayesPop} functions exactly the same way as an object with national projections. However, the meaning of the argument \code{country} in many of these functions (e.g. in \code{\link{pop.trajectories.plot}}) changes to an identification of the region (either as a numerical code or name as defined in the \code{locations} file).
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}

\section{Acknowledgment}{
We are greatful to Patrice Dion from Statistics Canada for providing us with example data. Note that the example datasets included in the package are not official STATCAN data - they only serve the purpose of illustration and templates.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pop.predict}}, \code{\link[bayesTFR]{tfr.predict.subnat}}, \code{\link{pop.aggregate.subnat}}
}
\examples{
\dontrun{
# Subnational projections for Canada
#########
data.dir <- file.path(find.package("bayesPop"), "extdata")

# Use national data for tfr and e0
###
sim.dir <- tempfile()
pred <- pop.predict.subnat(output.dir = sim.dir,
            locations = file.path(data.dir, "CANlocations.txt"),
            inputs = list(popM = file.path(data.dir, "CANpopM.txt"),
                          popF = file.path(data.dir, "CANpopF.txt"),
                          tfr.file = "median_"
                        ),
            verbose = TRUE)
pop.trajectories.plot(pred, "Alberta", sum.over.ages = TRUE)
unlink(sim.dir, recursive=TRUE)

# Use subnational TFR simulation
###
# Subnational TFR projections for Canada (from ?tfr.predict.subnat)
my.subtfr.file <- file.path(find.package("bayesTFR"), 'extdata', 'subnational_tfr_template.txt')
tfr.nat.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
tfr.reg.dir <- tempfile()
tfr.preds <- tfr.predict.subnat(124, my.tfr.file = my.subtfr.file,
    sim.dir = tfr.nat.dir, output.dir = tfr.reg.dir, start.year = 2013)
 
# Pop projections
sim.dir <- tempfile()
pred <- pop.predict.subnat(output.dir = sim.dir,
            locations = file.path(data.dir, "CANlocations.txt"),
            inputs = list(popM = file.path(data.dir, "CANpopM.txt"),
                          popF = file.path(data.dir, "CANpopF.txt"),
                          patterns = file.path(data.dir, "CANpatterns.txt"),
                          tfr.sim.dir = file.path(tfr.reg.dir, "subnat", "c124")
                        ),
            verbose = TRUE)
pop.trajectories.plot(pred, "Alberta", sum.over.ages = TRUE)
pop.pyramid(pred, "Manitoba", year = 2050)
get.countries.table(pred)

# Aggregate to country level
aggr <- pop.aggregate.subnat(pred, regions = 124, 
            locations = file.path(data.dir, "CANlocations.txt"))
pop.trajectories.plot(aggr, "Canada", sum.over.ages = TRUE)

unlink(sim.dir, recursive = TRUE)
unlink(tfr.reg.dir, recursive = TRUE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")

