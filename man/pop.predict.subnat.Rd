\name{pop.predict.subnat}
\alias{pop.predict.subnat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subnational Probabilistic Population Projection
}
\description{
Generates trajectories of probabilistic population projection for subregions of a given country.
}
\usage{
pop.predict.subnat(end.year = 2050, start.year = 1950, present.year = 2015, 
              wpp.year = 2017, output.dir = file.path(getwd(), "bayesPop.output"), 
              locations = NULL, default.country = NULL, 
              inputs = list(popM = NULL, popF = NULL, 
                            mxM = NULL, mxF = NULL, srb = NULL, 
                            pasfr = NULL, patterns = NULL, 
                            migM = NULL, migF = NULL, 
                            e0F.file = NULL, e0M.file = NULL, tfr.file = NULL, 
                            e0F.sim.dir = NULL, e0M.sim.dir = NULL, tfr.sim.dir = NULL, 
                            migMtraj = NULL, migFtraj = NULL), 
              nr.traj = 1000, keep.vital.events = FALSE, 
              fixed.mx = FALSE, fixed.pasfr = FALSE, 
              replace.output = FALSE, verbose = TRUE)
}

\arguments{
  \item{end.year}{End year of the projection.}
  \item{start.year}{First year of the historical data.}
  \item{present.year}{Year for which initial population data is to be used.}
  \item{wpp.year}{Year for which WPP data is used. The functions loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and uses the various datasets as default for the country corresponding to the regions (see below).}
  \item{output.dir}{Output directory of the projection.}
  \item{locations}{Name of a tab-delimited file that contains definitions of the regions. It has a similar structure as \code{\link[wpp2017]{UNlocations}}, with mandatory columns \code{reg_code} (unique identifier of the regions) and \code{name} (name of the regions). Optionally \code{location_type} should be set to 4 for regions to be processed. Column \code{country_code} can be included with the numerical code of the corresponding country. This is a mandatory argument. }
  \item{default.country}{Numerical code of a country to which the regions belong to. It is used for extracting default datasets from the  \pkg{wpp} package if some region-specific input datasets are missing. Alternatively, it can be also included in the \code{locations} file, where it would have \code{location_type} of 0. In such a case, the value is taken from its column \code{country_code}. In either case, the code must exists in the \code{\link[wpp2017]{UNlocations}} dataset.}
  \item{inputs}{A list of file names where input data is stored. Unless otherwise noted, these are tab delimited ASCII files with a mandatory column \code{reg_code} giving the numerical identifier of the regions. If an element of this list is \code{NULL}, usually a default dataset corresponding to \code{default.country} is extracted from the \pkg{wpp} package. Names of these default datasets are shown in brackets. This list contains the following elements:
  		\describe{
  			\item{popM, popF}{Initial male/female age-specific population (at time \code{present.year}). Mandatory items, no defaults. Must contain columns \code{reg_code} and \code{age} and be of the same structure as \code{\link[wpp2017]{popM}}.}
  			\item{mxM, mxF}{Historical data and (optionally) projections of male/female age-specific death rates [\code{\link[wpp2017]{mxM}}, \code{\link[wpp2017]{mxF}}] (see also argument \code{fixed.mx}).}
  			\item{srb}{Projection of sex ratio at birth. [\code{\link[wpp2017]{sexRatio}}]}
  			\item{pasfr}{Historical data and (optionally) projections of percentage age-specific fertility rate [\code{\link[wpp2017]{percentASFR}}] (see also argument \code{fixed.pasfr}).}
  			\item{patterns}{Information on region's specifics regarding migration type, base year of the migration, mortality and fertility age patterns as defined in [\code{\link{vwBaseYear}}]. In addition, it can contain columns defining migration shares between the regions. These are: \code{inmigrationM_share}, \code{inmigrationF_share}, \code{outmigrationM_share}, \code{outmigrationF_share}. The sex specification and/or direction specification (in/out) can be omitted, e.g. it can be simply \code{migration_share}. If sex is omitted, Rogers Castro age schedules are used for migration.}
  			\item{migM, migF}{Projection of male/female age-specific migration as net counts on the same scale as initital population [\code{\link[wpp2012]{migrationM}}, \code{\link[wpp2012]{migrationF}}]. If not available, the migration schedules are constructed from total migration counts of the \code{default.country} derived from \code{\link[wpp2017]{migration}} using Rogers Castro for age distribution. Migration shares between region (including sex-specific shares) can be given in the \code{patterns} file, see above.}
  			\item{e0F.file}{Comma-delimited CSV file with female life expectancy. It has the same structure as the file \dQuote{ascii_trajectories.csv} generated using \code{bayesLife::\link[bayesLife]{convert.e0.trajectories}} (which currently works for country-level results only). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{e0}. If this element is not \code{NULL}, the argument \code{e0F.sim.dir} is ignored. If both \code{e0F.file} and \code{e0F.sim.dir} are \code{NULL}, data from the corresponding \pkg{wpp} package (for \code{default.country}) is taken, namely the median projections as one trajectory and the low and high variants (if available) as second and third trajectory. Alternatively, this element can be the keyword \dQuote{median_} in which case only the median is taken.}
  			\item{e0M.file}{Comma-delimited CSV file containing results of male life expectancy of the same format as \code{e0F.file}. If this element is not \code{NULL}, the argument \code{e0M.sim.dir} is ignored. As in the female case, if both \code{e0M.file} and \code{e0M.sim.dir} are \code{NULL}, data for \code{default.country} from the corresponding \pkg{wpp} package is taken.}
  			\item{tfr.file}{Comma-delimited CSV file with results of total fertility rate (generated using \pkg{\link{bayesTFR}}, function \code{\link{convert.tfr.trajectories}}, file \dQuote{ascii_trajectories.csv}). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{TF}. If this element is not \code{NULL}, the argument \code{tfr.sim.dir} is ignored. If both \code{tfr.file} and \code{tfr.sim.dir} are \code{NULL}, data for \code{default.country} from the corresponding \pkg{wpp} package is taken (median and the low and high variants as three trajectories). Alternatively, this argument can be the keyword \dQuote{median_} in which case only the wpp median is taken.}
  			\item{e0F.sim.dir}{Simulation directory with results of female life expectancy. Since \pkg{\link{bayesLife}} does not support subnational projections yet, this element should not be used. Instead use \code{e0F.file} or leave it blank. Alternatively, it can be set to the keyword \dQuote{median_} which has the same effect as when \code{e0F.file} is \dQuote{median_}. This item is only used if \code{e0F.file} is \code{NULL}.}
  			\item{e0M.sim.dir}{This is analogous to \code{e0F.sim.dir}, here for male life expectancy. Use \code{e0M.file} instead of this item.}
  			\item{tfr.sim.dir}{Simulation directory with results of total fertility rate (generated using \code{bayesTFR::\link[bayesTFR]{tfr.predict.subnat}}). It is only used if \code{tfr.file} is \code{NULL}.}
  			\item{migMtraj, migFtraj}{Comma-delimited CSV file with male/female age-specific migration trajectories. If present, it replaces deterministic projections given by the \code{migM} and \code{migF} items. It has a similar format as e.g. \code{e0M.file} with columns \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, \dQuote{Age} and \dQuote{Migration}. The \dQuote{Age} column must have values \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dots, \dQuote{95-99}, \dQuote{100+}.}
  		}
}
  \item{nr.traj, keep.vital.events, fixed.mx, fixed.pasfr, replace.output, verbose}{These arguments have the same meaning as in \code{\link{pop.predict}}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
Object of class \code{\link{bayesPop.prediction}} containing the subnational projections. 
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pop.predict}}, \code{\link[bayesTFR]{tfr.predict.subnat}}
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")

