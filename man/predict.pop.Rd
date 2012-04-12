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
pop.predict(end.year = 2100, start.year = 1950, present.year = 2010, 
    wpp.year = 2010, countries = NULL, 
    output.dir = file.path(getwd(), "bayesPop.output"),
    inputs = list(popM=NULL, popF=NULL, mxM=NULL, mxF=NULL, srb=NULL,
        pasfr=NULL, mig.type=NULL, migM=NULL, migF=NULL,
        e0F.file=NULL, e0M.file=NULL, tfr.file=NULL,
        e0F.sim.dir=NULL, e0M.sim.dir=NULL, tfr.sim.dir=NULL), 
    nr.traj = 1000, replace.output = FALSE, verbose = TRUE)
}
\arguments{
	\item{end.year}{End year of the projection.}
  	\item{start.year}{First year of the historical data.}
  	\item{present.year}{Year for which initial population data is to be used.}
  	\item{wpp.year}{Year for which WPP data is available.}
  	\item{countries}{Array of country codes or country names for which a projection is generated. If it is \code{NULL}, all available countries are used. If it is \code{NA} and there is an existing projection in \code{output.dir} and \code{replace.output=FALSE}, then a projection is performed for all countries that are not included in the existing projection. Names of countries are matched to those in the \code{\link{LOCATIONS}} dataset.}
  	\item{output.dir}{Output directory of the projection. If there is an existing projection in \code{output.dir} and \code{replace.output=TRUE}, everything in the directory will be deleted.}
  	\item{inputs}{A list of file names where input data is stored. It contains the following elements (names of default datasets which are used if the corresponding element is \code{NULL} are shown in brackets):
  		\describe{
  			\item{popM}{Initial male age-specific population (at time \code{present.year}) [\code{\link{PopByAgeMale}}].}
  			\item{popF}{Initial female age-specific population (at time \code{present.year}) [\code{\link{PopByAgeFemale}}].}
  			\item{mxM}{Historical data of male age-specific death rates [\code{\link{ASMRMale}}].}
  			\item{mxF}{Historical data of female age-specific death rates [\code{\link{ASMRFemale}}].}
  			\item{srb}{Projection of sex ratio at birth. [\code{\link{SexRatioAtBirth}}]}
  			\item{pasfr}{Projection of percentage age-specific fertility rate. [\code{\link{PercentASFR}}]}
  			\item{mig.type}{Migration type and base year of the migration. [\code{\link{vwBaseYear}}]}
  			\item{migM}{Projection of male age-specific migration [\code{\link{MigByAgeMale}}].}
  			\item{migF}{Projection of female age-specific migration [\code{\link{MigByAgeFemale}}].}
  			\item{e0F.file}{CSV file with results of female life expectancy (generated using \pkg{\link{bayesLife}}, function \code{\link{convert.e0.trajectories}}, file \dQuote{ascii_trajectories.csv}). If it is not \code{NULL}, the argument \code{e0F.sim.dir} is ignored. If both \code{e0F.file} and \code{e0F.sim.dir} are \code{NULL}, file \dQuote{e0Fwpp2010.csv} is used located in the \dQuote{ex-data} directory of the package which contains 100 trajectories.}
  			\item{e0M.file}{CSV file containing results of male life expectancy (generated using \pkg{\link{bayesLife}}, function \code{\link{convert.e0.trajectories}}, file \dQuote{ascii_trajectories.csv}). If it is not \code{NULL}, the argument \code{e0M.sim.dir} is ignored. If both \code{e0M.file} and \code{e0M.sim.dir} are \code{NULL}, file \dQuote{e0Mwpp2010.csv} is used located in the \dQuote{ex-data} directory of the package which contains 100 trajectories.}
  			\item{tfr.file}{CSV file with results of total fertility rate (generated using \pkg{\link{bayesTFR}}, function \code{\link{convert.tfr.trajectories}}, file \dQuote{ascii_trajectories.csv}). If it is not \code{NULL}, the argument \code{tfr.sim.dir} is ignored. If both \code{tfr.file} and \code{tfr.sim.dir} are \code{NULL}, file \dQuote{TFRwpp2010.csv} is used located in the \dQuote{ex-data} directory of the package which contains 100 trajectories.}
  			\item{e0F.sim.dir}{Simulation directory with results of female life expectancy (generated using \pkg{\link{bayesLife}}). It is only used if \code{e0F.file} is \code{NULL}.}
  			\item{e0M.sim.dir}{Simulation directory with results of male life expectancy (generated using \pkg{\link{bayesLife}}). Alternatively, it can be the string \dQuote{joint_}, in which case it is assumed that the male life expectancy was projected jointly from the female life expectancy (see \link[bayesLife]{joint.male.predict}) and thus contained in the \code{e0F.sim.dir} directory. The argument is only used if \code{e0M.file} is \code{NULL}.}
  			\item{tfr.sim.dir}{Simulation directory with results of total fertility rate (generated using \pkg{\link{bayesTFR}}). It is only used if \code{tfr.file} is \code{NULL}.}
  		}
  	}
  	\item{nr.traj}{Number of trajectories to be generated. If any of the simulations in \code{e0M.file/e0M.sim.dir}, \code{e0F.file/e0F.sim.dir} and \code{tfr.file/tfr.sim.dir} contains less trajectories than \code{nr.traj}, the number is set to the minimum of those.}
  	\item{replace.output}{Logical. If \code{TRUE}, everything in the directory \code{output.dir} is deleted prior to the prediction.} 
  	\item{verbose}{Logical controlling the amount of output messages.} 
}
\details{
The population projection is computed using the Cohort Component method and is based on an algorithm used by the United Nation Population Division. For each country, one projection is calculated for each trajectory of male and female life expectancy and TFR. This results in a set of trajectories of population projection which forms its posterior distribution. The trajectories of life expectancy and TFR can be given either in its binary form generated by the packages \pkg{\link{bayesLife}} and \pkg{\link{bayesTFR}}, respectively (as directories \code{e0M.sim.dir}, \code{e0F.sim.dir}, \code{tfr.sim.dir} of the \code{inputs} argument), or they can be given as ASCII tables in csv format, see above.

The projection is generated sequentially country by country. Results are stored in a sub-directory of \code{output.dir} called \file{prediction}. There is one binary file per country, called \file{totpop_country\eqn{x}.rda}, where \eqn{x} is the country code. It contains six objects: \code{totp}, \code{totpf}, \code{totpm} (trajectories of total population, female and male, respectively),  \code{totp.hch}, \code{totpf.hch}, \code{totpm.hch} (the UN half-child variant for total population, female and male, respectively). Furthermore, an object of class \code{bayesPop.prediction} is stored in the same directory in a file \file{prediction.rda}. It is updated every time a country projection is finished.

To access a previously stored prediction object, use \code{\link{get.pop.prediction}}.
}
\value{
Object of class \code{bayesPop.prediction} with the following elements:
\item{output.directory}{Full path to the projections.}
\item{nr.traj}{The actual number of trajectories of the projections.}
\item{quantiles}{Three-dimensional array of projection quantiles (countries x number of quantiles x projection periods). The second dimension corresponds to the following quantiles: \eqn{0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1}.}
\item{traj.mean.sd}{Three-dimensional array of projection mean and standard deviation (countries x 2 x projection periods). First and second matrix of the second dimension, respectively, is the mean and standard deviation, respectively.}
\item{quantilesM, quantilesF}{Quantiles of male and female projection, respectively. Same structure as \code{quantiles}.}
\item{traj.mean.sdM, traj.mean.sdF}{Same as \code{traj.mean.sd} corresponding to male and female projection, respectively.}
\item{quantilesMage, quantilesFage}{Four-dimensional array of age-specific quantiles of male and female projection, respectively (countries x age groups x number of quantiles x projection periods). The same quantiles are used as in \code{quantiles}.}
\item{quantilesPropMage, quantilesPropFage}{Array of age-specific quantiles of male and female projection, respectively, divided by the total population. The dimensions are the same as in \code{quantilesMage}.}
\item{estim.years}{Vector of time for which historical data was used in the projections.}
\item{proj.years}{Vector of projection time periods starting with the present period.}
\item{inputs}{List of input data used for the projection.}
\item{countries}{Matrix of countries for which projection exists. It contains two columns: \code{code}, \code{name}.}
\item{ages}{Vector of age groups.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, based on code of Nan Li
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pop.trajectories.plot}}, \code{\link{pop.pyramid}}, \code{\link{get.pop.prediction}}
}
\examples{
sim.dir <- tempfile()
# Countries can be given as a combination of numerical codes and names
pred <- pop.predict(countries=c("Netherlands", 218, "Madagascar"), nr.traj=10, 
           output.dir=sim.dir)
pop.trajectories.plot(pred, "Ecuador", sum.over.ages=TRUE)
unlink(sim.dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }

