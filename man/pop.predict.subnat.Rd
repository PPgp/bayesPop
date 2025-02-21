\name{pop.predict.subnat}
\alias{pop.predict.subnat}

\title{
Subnational Probabilistic Population Projection
}
\description{
Generates trajectories of probabilistic population projection for subregions of a given country.
}
\usage{
pop.predict.subnat(end.year = 2060, start.year = 1950, present.year = 2020, 
        wpp.year = 2019, output.dir = file.path(getwd(), "bayesPop.output"), 
        locations = NULL, default.country = NULL, annual = FALSE,
        inputs = list(
            popM = NULL, popF = NULL, 
            mxM = NULL, mxF = NULL, srb = NULL, 
            pasfr = NULL, patterns = NULL, 
            migM = NULL, migF = NULL, 
            migMt = NULL, migFt = NULL, mig = NULL, mig.fdm = NULL,
            e0F.file = NULL, e0M.file = NULL, tfr.file = NULL, 
            e0F.sim.dir = NULL, e0M.sim.dir = NULL, tfr.sim.dir = NULL, 
            migMtraj = NULL, migFtraj = NULL, migtraj = NULL,
            migFDMtraj = NULL, GQpopM = NULL, GQpopF = NULL, 
            average.annual = NULL
        ), 
        nr.traj = 1000, keep.vital.events = FALSE, 
        fixed.mx = FALSE, fixed.pasfr = FALSE, lc.for.all = TRUE,
         mig.is.rate = FALSE, mig.age.method = c("rc", "fdmp", "fdmnop"),
         mig.rc.fam = NULL, pasfr.ignore.phase2 = FALSE, 
         replace.output = FALSE, verbose = TRUE)
}

\arguments{
  \item{end.year}{End year of the projection.}
  \item{start.year}{First year of the historical data on mortality rates. It determines the length of the historical time series used in the Lee-Carter estimation.}
  \item{present.year}{Year for which initial population data is to be used.}
  \item{wpp.year}{Year for which WPP data is used. The function loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and uses its data (corresponding to the \code{default.country}) as default datasets if region-specific alternatives are not given (see more details below).}
  \item{output.dir}{Output directory of the projection.}
  \item{locations}{Name of a tab-delimited file that contains definitions of the subregions. It has a similar structure as \code{\link[wpp2019]{UNlocations}}, with mandatory columns \code{reg_code} (unique identifier of the subregions) and \code{name} (name of the subregions). Optionally, \code{location_type} should be set to 4 for subregions to be processed. Column \code{country_code} can be included with the numerical code of the corresponding country. A row with \code{location_type} of 0 determines the country that the subregions belong to and is used for extracting default "national" datasets if the argument \code{default.country} is missing. In such a case, the code of the default country is taken from its column \code{country_code}. This is a mandatory argument. }
  \item{default.country}{Numerical code of a country to which the subregions belong to. It is used for extracting default datasets from the  \pkg{wpp} package if some region-specific input datasets are missing. Alternatively, it can be also included in the \code{locations} file, see above. In either case, the code must exists in the \code{\link[wpp2019]{UNlocations}} dataset.}
  \item{annual}{Logical. If \code{TRUE} it is assumed that this is 1x1 simulation, i.e. one year age groups and one year time periods.}
  \item{inputs}{A list of file names where input data is stored. Unless otherwise noted, these are tab delimited ASCII files with a mandatory column \code{reg_code} giving the numerical identifier of the subregions. If an element of this list is \code{NULL}, usually a default dataset corresponding to \code{default.country} is extracted from the \pkg{wpp} package. Names of these default datasets are shown in brackets. This list contains the following elements:
  		\describe{
  			\item{popM, popF}{Initial male/female age-specific population (at time \code{present.year}). Mandatory items, no defaults. Must contain columns \code{reg_code} and \code{age} and be of the same structure as \code{\link[wpp2019]{popM}} from \pkg{wpp}.}
  			\item{mxM, mxF}{Historical data and (optionally) projections of male/female age-specific death rates [\code{\link[wpp2019]{mxM}}, \code{\link[wpp2019]{mxF}}] (see also argument \code{fixed.mx}).}
  			\item{srb}{Projection of sex ratio at birth. [\code{\link[wpp2019]{sexRatio}}]}
  			\item{pasfr}{Historical data and (optionally) projections of percentage age-specific fertility rate [\code{\link[wpp2019]{percentASFR}}] (see also argument \code{fixed.pasfr}).}
  			\item{patterns}{Information on region's specifics regarding migration type, base year of the migration, mortality and fertility age patterns as defined in [\code{\link{vwBaseYear}}]. In addition, it can contain columns defining migration shares between the subregions, see Details below.}
  			\item{migM, migF, migMt, migFt, mig}{Projection and (optionally) historical data of net migration on the same scale as the initital population. There are three ways of defining this quantity, here in order of priority: 1. via \code{migM} and \code{migF} which should give male and female age-specific migration [\code{\link[wpp2012]{migrationM}}, \code{\link[wpp2012]{migrationF}}]; 2. via \code{migMt} and \code{migFt} which should give male and female total net migration; 3. via \code{mig} which should give the total net migration. For 2. and 3., the totals are disagregated into age-specific migration by applying a Rogers-Castro schedule. For 3., the totals are equally split between sexes. If all of these input items are missing, the migration schedules are constructed from total migration counts of the \code{default.country} derived from \code{\link[wpp2019]{migration}} using Rogers Castro for age distribution. Migration shares between subregions (including sex-specific shares) can be given in the \code{patterns} file, see above and Details below. If no shares are given, it is distributed by population shares.}
  			\item{mig.fdm}{If \code{mig.age.method} is \dQuote{fdmp} or \dQuote{fdmnop}, this file is used to disaggregate total in- and out-migration into ages, giving proportions of the migration in-flow and out-flow for each age. It should have columns \dQuote{reg_code}, \dQuote{age}, \dQuote{in} and \dQuote{out}, where the latter two should each sum to 1 for each location. By default Rogers-Castro curves are used, obtained via the  function \code{\link{rcastro.schedule}}.}
  			\item{e0F.file}{Comma-delimited CSV file with projected female life expectancy. It has the same structure as the file \dQuote{ascii_trajectories.csv} generated using \code{bayesLife::\link[bayesLife]{convert.e0.trajectories}} (which currently works for country-level results only). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{e0}. If \code{e0F.file} is \code{NULL}, data from the corresponding \pkg{wpp} package (for \code{default.country}) is taken, namely the median projections as one trajectory and the low and high variants (if available) as second and third trajectory. Alternatively, this element can be the keyword \dQuote{median_} in which case only the median is taken.}
  			\item{e0M.file}{Comma-delimited CSV file containing projections of male life expectancy of the same format as \code{e0F.file}. As in the female case, if \code{e0M.file} is \code{NULL}, data for \code{default.country} from the corresponding \pkg{wpp} package is taken.}
  			\item{tfr.file}{Comma-delimited CSV file with results of total fertility rate (generated using \pkg{\link[bayesTFR]{bayesTFR}}, function \code{\link[bayesTFR]{convert.tfr.trajectories}}, file \dQuote{ascii_trajectories.csv}). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{TF}. If this element is not \code{NULL}, the argument \code{tfr.sim.dir} is ignored. If both \code{tfr.file} and \code{tfr.sim.dir} are \code{NULL}, data for \code{default.country} from the corresponding \pkg{wpp} package is taken (median and the low and high variants as three trajectories). Alternatively, this argument can be the keyword \dQuote{median_} in which case only the wpp median is taken.}
  			\item{e0F.sim.dir}{Simulation directory with results of female life expectancy, generated using \code{bayesLife::\link[bayesLife]{e0.predict.subnat}}. It is only used if \code{e0F.file} is \code{NULL}. Alternatively, it can be set to the keyword \dQuote{median_} which has the same effect as when \code{e0F.file} is \dQuote{median_}.}
  			\item{e0M.sim.dir}{This is analogous to \code{e0F.sim.dir}, here for male life expectancy. Use \code{e0M.file} instead of this item.}
  			\item{tfr.sim.dir}{Simulation directory with projections of total fertility rate (generated using \code{bayesTFR::\link[bayesTFR]{tfr.predict.subnat}}). It is only used if \code{tfr.file} is \code{NULL}.}
  			\item{migMtraj, migFtraj, migtraj}{Comma-delimited CSV file with male/female age-specific migration trajectories, or total migration trajectories (\code{migtraj}). If present, it replaces deterministic projections given by the \code{mig*} items. It has a similar format as e.g. \code{e0M.file} with columns \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, \dQuote{Age} (except for \code{migtraj}) and \dQuote{Migration}. For a five-year simulation, the \dQuote{Age} column must have values \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dots, \dQuote{95-99}, \dQuote{100+}. In an annual simulation, age is given by a single number between 0 and 100.}
  			\item{migFDMtraj}{Comma-delimited CSV file with trajectories of in- and out-migration schedules used for the FDM migration method, i.e. if \code{mig.age.method} is \dQuote{fdmp} or \dQuote{fdmnop}. The values have te same meaning as in the \code{mig.fdm} input item, except that here multiple trajectories of such schedules can be provided. It should have columns \dQuote{LocID}, \dQuote{Age}, \dQuote{Trajectory}, \dQuote{Value}, and \dQuote{Parameter}. For \dQuote{Age}, the same rules apply as for \code{migMtraj} above. The \dQuote{Parameter} column should have values \dQuote{in} for in-migration, \dQuote{out} for out-migration and \dQuote{v} for values of the variance denominator \eqn{v} used in Equation 22 of \href{https://arxiv.org/abs/2411.09878}{Sevcikova et al (2024)}. For the \eqn{v} parameter, the \dQuote{Age} column should be left empty.}
  			\item{GQpopM, GQpopF}{Age-specific population counts (male and female) that should be excluded from application of the cohort-component method (CCM). It can be used for defining group quarters. These counts are removed from population before the CCM projection and added back afterwards. It is not used when computing vital events on observed data. The datasets should have columns \dQuote{reg_code}, \dQuote{age} and \dQuote{gq}. In such a case the \dQuote{gq} amount is applied to all years. If it is desired to destinguish the amount that is added back for individual years, the \dQuote{gq} column should be replaced by columns indicating the individual years, i.e. single years for an annual simulation and time periods (e.g. \dQuote{2020-2025}, \dQuote{2025-2030}) for a 5-year simulation. For a five-year simulation, the \dQuote{age} column should include values \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dots, \dQuote{95-99}, \dQuote{100+}. However, rows with zeros do not need to be included. In an annual simulation, age is given by a single number between 0 and 100.}
  			\item{average.annual}{Character string with values \dQuote{TFR}, \dQuote{e0M}, \dQuote{e0F}. If this is a 5-year simulation, but the inputs of TFR or/and e0 comes from an annual simulation, including the corresponding string here will cause that the TFR or/and e0 trajectories are converted into 5-year averages.}
  		}
}
  \item{nr.traj, keep.vital.events, fixed.mx, fixed.pasfr, lc.for.all, mig.is.rate, mig.age.method, mig.rc.fam, replace.output, verbose}{These arguments have the same meaning as in \code{\link{pop.predict}}.}
  \item{pasfr.ignore.phase2}{Logical. If \code{TRUE} the TFR for all locations is considered being in phase III when predicting PASFR.}
}

\details{
Population projection for subnational units (regions) is performed by applying the cohort component method to subnational datasets on projected fertility (TFR), mortality and net migration, starting from given sex- and age-specific population counts. The only required inputs are the initial sex- and age-specific population counts in each region (\code{popM} and \code{popF} elements of the \code{inputs} argument) and a file with a set of locations (argument \code{locations}). If no other input datasets are given, those datasets are replaced by the corresponding "national" values, taken from the corresponding \pkg{wpp} package. The argument \code{default.country} determines the country for those default "national" values. The default country can be also included in the \code{locations} file as a record with \code{location.type} being set to 0.

The TFR component can be given as a set of trajectories generated using the \code{\link[bayesTFR]{tfr.predict.subnat}} function of the \pkg{bayesTFR} package (\code{tfr.sim.dir} element). Alternatively, trajectories can be given in an ASCII file (\code{tfr.file}). 

Similarly, the $e_0$ component can be given as a set of trajectories using the \code{\link[bayesLife]{e0.predict.subnat}} function of the \pkg{bayesLife} package (\code{e0F.sim.dir} element). If male projections are generated jointly (i.e. \code{predict.jmale = TRUE}), set \code{e0M.sim.dir = "joint_"}. Alternatively, trajectories can be given in an ASCII files (\code{e0F.file}, \code{e0M.file}).

Having a set of subnational TFR and $e_0$ trajectories, the cohort component method is applied to each of them to yield a distribution of future subnational population.

Projection of net migration can either be given as disaggregated sex- and age-specific datasets (\code{migM} and \code{migF}), or as sex totals (\code{migMt} and \code{migFt}), or as totals (\code{mig}), or as sex- and age-specific trajectories (\code{migMtraj} and \code{migFtraj}), or as total trajectories (\code{migtraj}). Alternatively, it can be given as shares between regions as columns in the \code{patterns} dataset. These are: \code{inmigrationM_share}, \code{inmigrationF_share}, \code{outmigrationM_share}, \code{outmigrationF_share}. The sex specification and/or direction specification (in/out) can be omitted, e.g. it can be simply \code{migration_share}. The function extracts the values of net migration projection on the national level and distributes it to regions according to the given shares. For positive (national) values, it uses the in-migration shares; for negative values it uses the out-migration shares. If the in/out prefix is omitted in the column names, the given migartion shares are used for both, positive and negative net migration projection. By default, if no migration datasets neither region-specific shares are given, the distribution between regions is proportional to the size of population. The age-specific schedules follow by default the Rogers-Castro age schedules. Note that when handling migration using shares as described here, it only affects the distribution of international migration into regions. It does not take into account between-region migration.

The package contains example datasets for Canada. Use these as templates for your own data. See Example below.
}

\value{
Object of class \code{\link{bayesPop.prediction}} containing the subnational projections. Note that this object can be used in the various \pkg{bayesPop} functions exactly the same way as an object with national projections. However, the meaning of the argument \code{country} in many of these functions (e.g. in \code{\link{pop.trajectories.plot}}) changes to an identification of the region (either as a numerical code or name as defined in the \code{locations} file).
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}

\section{Acknowledgment}{
We are greatful to Patrice Dion from Statistics Canada for providing us with example data. Note that the example datasets included in the package are not official STATCAN data - they only serve the purpose of illustration and templates. Data for the time period 2015-2020 has been imputed by the author.
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

# Use subnational probabilistic TFR simulation
###
# Subnational TFR projections for Canada (from ?tfr.predict.subnat)
my.subtfr.file <- file.path(find.package("bayesTFR"), 'extdata', 'subnational_tfr_template.txt')
tfr.nat.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
tfr.reg.dir <- tempfile()
tfr.preds <- tfr.predict.subnat(124, my.tfr.file = my.subtfr.file,
    sim.dir = tfr.nat.dir, output.dir = tfr.reg.dir, start.year = 2013)
 
# Use subnational probabilistic e0
### 
# Subnational e0 projections for Canada (from ?e0.predict.subnat)
# (here using the same female and male data, just for illustration)
my.sube0.file <- file.path(find.package("bayesLife"), 'extdata', 'subnational_e0_template.txt')
e0.nat.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
e0.reg.dir <- tempfile()
e0.preds <- e0.predict.subnat(124, my.e0.file = my.sube0.file,
    sim.dir = e0.nat.dir, output.dir = e0.reg.dir, start.year = 2018,
    predict.jmale = TRUE, my.e0M.file = my.sube0.file)
 
# Population projections
sim.dir <- tempfile()
pred <- pop.predict.subnat(output.dir = sim.dir,
            locations = file.path(data.dir, "CANlocations.txt"),
            inputs = list(popM = file.path(data.dir, "CANpopM.txt"),
                          popF = file.path(data.dir, "CANpopF.txt"),
                          patterns = file.path(data.dir, "CANpatterns.txt"),
                          tfr.sim.dir = file.path(tfr.reg.dir, "subnat", "c124"),
                          e0F.sim.dir = file.path(e0.reg.dir, "subnat_ar1", "c124"),
                          e0M.sim.dir = "joint_"
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
unlink(e0.reg.dir, recursive = TRUE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")

