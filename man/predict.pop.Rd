\name{pop.predict}
\alias{pop.predict}
\alias{bayesPop.prediction}

\title{
Probabilistic Population Projection
}
\description{
The function generates trajectories of probabilistic population projection for all countries for which input data is available, or any subset of them.
}
\usage{
pop.predict(end.year = 2100, start.year = 1950, present.year = 2020, 
    wpp.year = 2019, countries = NULL, 
    output.dir = file.path(getwd(), "bayesPop.output"),
    annual = FALSE,
    inputs = list(popM=NULL, popF=NULL, mxM=NULL, mxF=NULL, srb=NULL,
        pasfr=NULL, patterns=NULL, 
        migM=NULL, migF=NULL, migMt=NULL, migFt=NULL, mig=NULL,
        mig.fdm = NULL, e0F.file=NULL, e0M.file=NULL, tfr.file=NULL,
        e0F.sim.dir=NULL, e0M.sim.dir=NULL, tfr.sim.dir=NULL,
        migMtraj = NULL, migFtraj = NULL, migtraj = NULL,
        migFDMtraj = NULL, GQpopM = NULL, GQpopF = NULL, 
        average.annual = NULL), 
    nr.traj = 1000, keep.vital.events = FALSE, 
    fixed.mx = FALSE, fixed.pasfr = FALSE,
    lc.for.hiv = TRUE, lc.for.all = TRUE, mig.is.rate = FALSE,
    mig.age.method  = c("auto", "fdmp", "fdmnop", "rc"), mig.rc.fam = NULL,
    my.locations.file = NULL, replace.output = FALSE, verbose = TRUE, \dots)
}
\arguments{
	\item{end.year}{End year of the projection.}
  	\item{start.year}{First year of the historical data.}
  	\item{present.year}{Year for which initial population data is to be used.}
  	\item{wpp.year}{Year for which WPP data is used. The functions loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and uses the various datasets as default if the corresponding \code{inputs} element is missing (see below).}
  	\item{countries}{Array of country codes or country names for which a projection is generated. If it is \code{NULL}, all available countries are used. If it is \code{NA} and there is an existing projection in \code{output.dir} and \code{replace.output=FALSE}, then a projection is performed for all countries that are not included in the existing projection. Names of countries are matched to those in the \code{\link[wpp2019]{UNlocations}} dataset (or in the dataset loaded from \code{my.locations.file} if used).}
  	\item{output.dir}{Output directory of the projection. If there is an existing projection in \code{output.dir} and \code{replace.output=TRUE}, everything in the directory will be deleted.}
  	\item{annual}{Logical. If \code{TRUE} it is assumed that this is 1x1 simulation, i.e. one year age groups and one year time periods. Note that this is still an experimental feature!}
  	\item{inputs}{A list of file names where input data is stored. It contains the following elements (Unless otherwise noted, these are tab delimited ASCII files; Names of default datasets from the corresponding \pkg{wpp} package which are used if the corresponding element is \code{NULL} are shown in brackets):
  		\describe{
  			\item{popM, popF}{Initial male/female age-specific population (at time \code{present.year}) [\code{\link[wpp2019]{popM}}, \code{\link[wpp2019]{popF}}].}
  			\item{mxM, mxF}{Historical data and (optionally) projections of male/female age-specific death rates [\code{\link[wpp2019]{mxM}}, \code{\link[wpp2019]{mxF}}] (see also argument \code{fixed.mx}).}
  			\item{srb}{Projection of sex ratio at birth. [\code{\link[wpp2019]{sexRatio}}]}
  			\item{pasfr}{Historical data and (optionally) projections of percentage age-specific fertility rate [\code{\link[wpp2019]{percentASFR}}] (see also argument \code{fixed.pasfr}).}
  			\item{patterns, mig.type}{Migration type and base year of the migration. In addition, this dataset gives information on country's specifics regarding mortality and fertility age patterns as defined in [\code{\link{vwBaseYear}}]. \code{patterns} and \code{mig.type} have the same meaning and can be used interchangeably.}
  			\item{migM, migF, migMt, migFt, mig}{Projection and (optionally) historical data of net migration on the same scale as the initital population. There are three ways of defining this quantity, here in order of priority: 1. via \code{migM} and \code{migF} which should give male and female age-specific migration [\code{\link[wpp2012]{migrationM}}, \code{\link[wpp2012]{migrationF}}]; 2. via \code{migMt} and \code{migFt} which should give male and female total net migration; 3. via \code{mig} which should give the total net migration. For 2. and 3., the totals are disagregated into age-specific migration by applying a schedule defined by the \code{mig.age.method} argument. If all of these input items are missing, for \code{wpp.year = 2024} or 2012, the UN age schedules are used. For other WPP revisions, the migration schedules are reconstructed from total migration counts derived from \code{\link[wpp2019]{migration}} using either the \code{\link{age.specific.migration}} or the \code{\link{migration.totals2age}} function.}
  			\item{mig.fdm}{If \code{mig.age.method} is \dQuote{fdmp} or \dQuote{fdmnop}, this file is used to disaggregate total in- and out-migration into ages, giving proportions of the migration in-flow and out-flow for each age. It should have columns \dQuote{country_code}, \dQuote{age}, \dQuote{in} and \dQuote{out}, where the latter two should each sum to 1 for each location. By default the function uses the \code{\link{rc1FDM}} (annual) or \code{\link{rc5FDM}} (5-year) datasets. For locations where the unique identifier does not match the country code in these default datasets, Rogers-Castro curves are used, obtained via the  function \code{\link{rcastro.schedule}}.}
  			\item{e0F.file}{Comma-delimited CSV file with results of female life expectancy (generated using \pkg{\link[bayesLife]{bayesLife}}, function \code{\link[bayesLife]{convert.e0.trajectories}}, file \dQuote{ascii_trajectories.csv}). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{e0}. If this element is not \code{NULL}, the argument \code{e0F.sim.dir} is ignored. If both \code{e0F.file} and \code{e0F.sim.dir} are \code{NULL}, data from the corresponding \pkg{wpp} package is taken, namely the median projections as one trajectory and the low and high variants (if available) as second and third trajectory. For 5-year simulations, column \dQuote{Year} should be the middle year of the time period, e.g. 2023, 2028 etc. }
  			\item{e0M.file}{Comma-delimited CSV file containing results of male life expectancy (generated using \pkg{\link[bayesLife]{bayesLife}}, function \code{\link[bayesLife]{convert.e0.trajectories}}, file \dQuote{ascii_trajectories.csv}). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{e0}. If this element is not \code{NULL}, the argument \code{e0M.sim.dir} is ignored. As in the female case, if both \code{e0M.file} and \code{e0M.sim.dir} are \code{NULL}, data from the corresponding \pkg{wpp} package is taken.}
  			\item{tfr.file}{Comma-delimited CSV file with results of total fertility rate (generated using \pkg{\link[bayesTFR]{bayesTFR}}, function \code{\link[bayesTFR]{convert.tfr.trajectories}}, file \dQuote{ascii_trajectories.csv}). Required columns are \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, and \dQuote{TF}. If this element is not \code{NULL}, the argument \code{tfr.sim.dir} is ignored. If both \code{tfr.file} and \code{tfr.sim.dir} are \code{NULL}, data from the corresponding \pkg{wpp} package is taken (median and the low and high variants as three trajectories). Alternatively, this argument can be the keyword \dQuote{median_} in which case only the wpp median is taken.}
  			\item{e0F.sim.dir}{Simulation directory with results of female life expectancy (generated using \pkg{\link[bayesLife]{bayesLife}}). It is only used if \code{e0F.file} is \code{NULL}.}
  			\item{e0M.sim.dir}{Simulation directory with results of male life expectancy (generated using \pkg{\link[bayesLife]{bayesLife}}). Alternatively, it can be the string \dQuote{joint_}, in which case it is assumed that the male life expectancy was projected jointly from the female life expectancy (see \link[bayesLife]{joint.male.predict}) and thus contained in the \code{e0F.sim.dir} directory. The argument is only used if \code{e0M.file} is \code{NULL}.}
  			\item{tfr.sim.dir}{Simulation directory with results of total fertility rate (generated using \pkg{\link[bayesTFR]{bayesTFR}}). It is only used if \code{tfr.file} is \code{NULL}.}
  			\item{migMtraj, migFtraj, migtraj}{Comma-delimited CSV file with male/female age-specific migration trajectories, or total migration trajectories (\code{migtraj}). If present, it replaces deterministic projections given by the \code{mig*} items. It has a similar format as e.g. \code{e0M.file} with columns \dQuote{LocID}, \dQuote{Year}, \dQuote{Trajectory}, \dQuote{Age} (except for \code{migtraj}) and \dQuote{Migration}. For a five-year simulation, the \dQuote{Age} column must have values \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dots, \dQuote{95-99}, \dQuote{100+}, and the \dQuote{Year} column should be the middle year of the time period, e.g. 2023, 2028 etc. In an annual simulation, age is given by a single number between 0 and 100, and \dQuote{Year} contains all projected years.}
  			\item{migFDMtraj}{Comma-delimited CSV file with trajectories of in- and out-migration schedules used for the FDM migration method, i.e. if \code{mig.age.method} is \dQuote{fdmp} or \dQuote{fdmnop}. The values have te same meaning as in the \code{mig.fdm} input item, except that here multiple trajectories of such schedules can be provided. It should have columns \dQuote{LocID}, \dQuote{Age}, \dQuote{Trajectory}, \dQuote{Value}, and \dQuote{Parameter}. For \dQuote{Age}, the same rules apply as for \code{migMtraj} above. The \dQuote{Parameter} column should have values \dQuote{in} for in-migration, \dQuote{out} for out-migration and \dQuote{v} for values of the variance denominator \eqn{v} used in Equation 22 of \href{https://arxiv.org/abs/2411.09878}{Sevcikova et al (2024)}. For the \eqn{v} parameter, the \dQuote{Age} column should be left empty.}
  			\item{GQpopM, GQpopF}{Age-specific population counts (male and female) that should be excluded from application of the cohort component method (CCM). It can be used for defining group quarters. These counts are removed from population before the CCM projection and added back afterwards. It is not used when computing vital events on observed data. The datasets should have columns \dQuote{country_code}, \dQuote{age} and \dQuote{gq}. In such a case the \dQuote{gq} amount is applied to all years. If it is desired to destinguish the amount that is added back for individual years, the \dQuote{gq} column should be replaced by columns indicating the individual years, i.e. single years for an annual simulation and time periods (e.g. \dQuote{2020-2025}, \dQuote{2025-2030}) for a 5-year simulation. For a five-year simulation, the \dQuote{age} column should include values \dQuote{0-4}, \dQuote{5-9}, \dQuote{10-14}, \dots, \dQuote{95-99}, \dQuote{100+}. However, rows with zeros do not need to be included. In an annual simulation, age is given by a single number between 0 and 100.}
  		\item{average.annual}{Character string with values \dQuote{TFR}, \dQuote{e0M}, \dQuote{e0F}. If this is a 5-year simulation, but the inputs of TFR or/and e0 comes from an annual simulation, including the corresponding string here will cause that the TFR or/and e0 trajectories are converted into 5-year averages.}
  	}
  	}
  	\item{nr.traj}{Number of trajectories to be generated. If this number is smaller than the number of available trajectories of the probabilistic components (TFR, life expectancy and migration), the trajectories are equidistantly thinned. 
  		If all of those components contain less trajectories than \code{nr.traj}, the value is adjusted to the maximum of available trajectories of the components. For those that have less trajectories than the adjusted number, the available trajectories are re-sampled, so that all components have the same number of trajectories.}
  	\item{keep.vital.events}{Logical. If \code{TRUE} age- and sex-specific vital events of births and deaths as well as other objects are stored in the prediction object, see Details.}
  	\item{fixed.mx}{Logical. If \code{TRUE}, it is assumed the dataset of death rates (mxM and mxF) include data for projection years and they are then used instead of the life expectancy.}
  	\item{fixed.pasfr}{Logical. If \code{TRUE}, it is assumed the dataset on percent age-specific fertility rate (percentASFR) include data for projection years and they are then used instead of computing it on the fly.}
  	\item{lc.for.hiv}{Logical controlling if the modified Lee-Carter method should be used 
  	for projection of mortality rates for countries with HIV epidemics. If \code{FALSE}, the function \code{hiv.mortmod} from the HIV.LifeTables package is used.}
  	\item{lc.for.all}{Logical controlling if the modified Lee-Carter method should be used 
  	for projection of mortality rates for all countries. If \code{FALSE}, the corresponding method is determined by the columns \dQuote{AgeMortProjMethod1} and \dQuote{AgeMortProjMethod2} of the \code{\link{vwBaseYear}} dataset.}
  	\item{mig.is.rate}{Logical determining if migration data are to be interpreted as net migration rates (\code{TRUE}) or counts (\code{FALSE}, default). It can also be a vector of two logicals, where the first element refers to observed data and the second element refers to predictions. A value of \code{c(FALSE, TRUE)} could for example be used if observed data in \code{inputs$mig} are counts, and migration trajectories in \code{inputs$migtraj} are rates.}
  	\item{mig.age.method}{If migration is given as totals, this argument determines a method to disaggregate into age-specific migration. 
  	
  	The \dQuote{rc} method uses a simple Rogers-Castro disaggregation, via the function \code{\link{rcastro.schedule}}. An alternative schedule can be passed via the \code{mig.rc.fam} argument.
  	
  	Values \dQuote{fdmp} and \dQuote{fdmnop} trigger the Flow Difference Method (Sevcikova et al, 2024), where \dQuote{fdmp} weights the flows by population, while \dQuote{fdmnop} is an unweighted version. They both split the total net migration into total in- and out-migration and then disaggregate these flows separately. These two FDM methods use additional inputs in the \code{inputs$rc.fdm} and/or \code{inputs$migFDMtraj} components. 
  	
  	The \dQuote{auto} method (default) uses \dQuote{rc} if sex-specific migration totals are given, i.e. in \code{inputs$migFt} and \code{inputs$migMt}. If \code{annual} is \code{FALSE} and \code{wpp.year} is 2015, 2017 or 2019, then the residual method using the function \code{\link{age.specific.migration}} is used. Otherwise the \dQuote{fdmp} method is applied.}
  	\item{mig.rc.fam}{Data frame providing a single family of Rogers-Castro parameters to be used if \code{mig.age.method} is set to \dQuote{rc}. Mandatory columns are \dQuote{age} and \dQuote{prop}. Optionally, it can have a column \dQuote{mig_sign} with values \dQuote{Inmigration} and \dQuote{Emigration} (distinguishing schedules to be applied for positive and negative migration, respectively) and a column \dQuote{sex} with values \dQuote{Female} and \dQuote{Male}. The format corresponds to the dataset \code{DemoTools::mig_un_families}, subset to a single family. If this argument is \code{NULL} and \code{mig.age.method = "rc"}, the function \code{\link{rcastro.schedule}} with equal sex ratio is used to distribute total migration into ages.}
  	\item{my.locations.file}{Name of a tab-delimited ascii file with a set of all locations for which a projection is generated. Use this argument if you are projecting for a country/region that is not included in the standard \code{\link[wpp2019]{UNlocations}} dataset. It must have the same structure.}
  	\item{replace.output}{Logical. If \code{TRUE}, everything in the directory \code{output.dir} is deleted prior to the prediction.} 
  	\item{verbose}{Logical controlling the amount of output messages.} 
  	\item{\dots}{Additional arguments passed to the underlying function. These can be \code{parallel} and \code{nr.nodes} for parallel processing and the number of nodes, respectively, as well as further arguments passed for creating a parallel cluster.}
}
\details{
The population projection is computed using the cohort component method and is based on an algorithm used by the United Nation Population Division (see also  Sevcikova et al (2016b) in the References below). For each country, one projection is calculated for each trajectory of male and female life expectancy, TFR and possibly migration. This results in a set of trajectories of population projection which forms its posterior distribution. The trajectories of life expectancy and TFR can be given either in its binary form generated by the packages \pkg{\link[bayesLife]{bayesLife}} and \pkg{\link[bayesTFR]{bayesTFR}}, respectively (as directories \code{e0M.sim.dir}, \code{e0F.sim.dir}, \code{tfr.sim.dir} of the \code{inputs} argument), or they can be given as ASCII tables in csv format, see above. The number of trajectories for male and female life expectancy must match, as does for male and female migration.

The projection is generated sequentially location by location. Results are stored in a sub-directory of \code{output.dir} called \file{prediction}. There is one binary file per location, called \file{totpop_country\eqn{x}.rda}, where \eqn{x} is the country code. It contains six objects: \code{totp}, \code{totpf}, \code{totpm} (trajectories of total population, age-specific female and age-specific male, respectively),  \code{totp.hch}, \code{totpf.hch}, \code{totpm.hch} (the UN half-child variant for total population, age-specific female and age-specific male, respectively). Optionally, if \code{keep.vital.events} is \code{TRUE}, there is an additional file per country, called \file{vital_events_country\eqn{x}.rda}, containing the following objects: \code{btm}, \code{btf} (trajectories for births by age of mothers for male and female child, respectively), \code{deathsm}, \code{deathsf} (trajectories for age-specific male and female deaths, respectively), \code{asfert} (trajectories of age-specific fertility), \code{mxm}, \code{mxf} (trajectories of male and female age-specific mortality rates), \code{migm}, \code{migf} (if used, these are trajectories of male and female age-specific migration), \code{btm.hch}, \code{btf.hch}, \code{deathsm.hch}, \code{deathsf.hch}, \code{asfert.hch}, \code{mxm.hch}, \code{mxf.hch} (the UN half-child variant for age- and sex-specific births, deaths, fertility rates and mortality rates).  An object of class \code{bayesPop.prediction} is stored in the same directory in a file \file{prediction.rda}. It is updated every time a country projection is finished.

See \code{\link{pop.trajectories}} for extracting trajectories.

To access a previously stored prediction object, use \code{\link{get.pop.prediction}}.
}
\value{
Object of class \code{bayesPop.prediction} with the following elements:
\item{base.directory}{Full path to the base directory \code{output.dir}.}
\item{output.directory}{Sub-directory relative to \code{base.directory} with the projections.}
\item{nr.traj}{The actual number of trajectories of the projections.}
\item{quantiles}{Three-dimensional array of projection quantiles (countries x number of quantiles x projection periods). The second dimension corresponds to the following quantiles: \eqn{0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975}.}
\item{traj.mean.sd}{Three-dimensional array of projection mean and standard deviation (countries x 2 x projection periods). First and second matrix of the second dimension, respectively, is the mean and standard deviation, respectively.}
\item{quantilesM, quantilesF}{Quantiles of male and female projection, respectively. Same structure as \code{quantiles}.}
\item{traj.mean.sdM, traj.mean.sdF}{Same as \code{traj.mean.sd} corresponding to male and female projection, respectively.}
\item{quantilesMage, quantilesFage}{Four-dimensional array of age-specific quantiles of male and female projection, respectively (countries x age groups x number of quantiles x projection periods). The same quantiles are used as in \code{quantiles}.}
\item{quantilesPropMage, quantilesPropFage}{Array of age-specific quantiles of male and female projection, respectively, divided by the total population. The dimensions are the same as in \code{quantilesMage}.}
\item{estim.years}{Vector of time for which historical data was used in the projections.}
\item{proj.years}{Vector of projection time periods starting with the present period.}
\item{wpp.year}{The wpp year used.}
\item{inputs}{List of input data used for the projection.}
\item{function.inputs}{Content of the \code{inputs} argument passed to the function.}
\item{countries}{Matrix of countries for which projection exists. It contains two columns: \code{code}, \code{name}.}
\item{ages}{Vector of age groups.}
\item{annual}{If \code{TRUE}, this object corresponds to a 1x1 prediction, otherwise 5x5.}
\item{cache}{This component is added by \code{\link{get.pop.prediction}} and modified and used by \code{\link{pop.map}} and \code{\link{write.pop.projection.summary}}. It is an environment for caching and re-using results of expressions.}
\item{write.to.cache}{Logical determining if \code{cache} should be modified.}
\item{is.aggregation}{Logical determining if this object is a result of \code{pop.predict} or \code{\link{pop.aggregate}}.}
}
\references{
H. Sevcikova, A. E. Raftery (2016a). bayesPop: Probabilistic
  Population Projections. Journal of Statistical Software, 75(5), 1-29.
  doi:10.18637/jss.v075.i05

A. E. Raftery,  N. Li, H. Sevcikova , P. Gerland, G. K. Heilig (2012). Bayesian probabilistic population projections for all countries. Proceedings of the National Academy of Sciences 109:13915-13921.

P. Gerland, A. E. Raftery, H. Sevcikova, N. Li, D. Gu, T. Spoorenberg, L. Alkema, B. K. Fosdick, J. L. Chunn, N. Lalic, G. Bay, T. Buettner, G. K. Heilig,  J. Wilmoth (2014). World Population Stabilization Unlikely This Century. Science 346:234-237.

H. Sevcikova, N. Li, V. Kantorova, P. Gerland and A. E. Raftery (2016b). Age-Specific Mortality and Fertility Rates for Probabilistic Population Projections. In: Dynamic Demographic Analysis, ed. Schoen R. (Springer), pp. 285-310. Earlier version in \href{https://arxiv.org/abs/1503.05215}{arXiv:1503.05215}.

H. Sevcikova, J. Raymer J., A. E. Raftery (2024). Forecasting Net Migration By Age: The Flow-Difference Approach. \href{https://arxiv.org/abs/2411.09878}{arXiv:2411.09878}.


%Working Paper 150. Center for Statistics and the Social Sciences, University of Washington. \url{http://www.csss.washington.edu/Papers/2015/wp150.pdf} 
}
\author{
Hana Sevcikova, Thomas Buettner, based on code of Nan Li and helpful comments from Patrick Gerland
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pop.trajectories.plot}}, \code{\link{pop.pyramid}}, \code{\link{pop.trajectories}}, \code{\link{get.pop.prediction}}, \code{\link{age.specific.migration}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
# Countries can be given as a combination of numerical codes and names
pred <- pop.predict(countries=c("Netherlands", 218, "Madagascar"), nr.traj=3, 
           output.dir=sim.dir)
pop.trajectories.plot(pred, "Ecuador", sum.over.ages=TRUE)
unlink(sim.dir, recursive=TRUE)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }

