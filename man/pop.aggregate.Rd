\name{pop.aggregate}
\alias{pop.aggregate}
\alias{get.pop.aggregation}

\title{
Aggregation of Population Projections
}
\description{
Aggregation of existing countries' population projections into projections of given regions, and accessing such aggregations.
}
\usage{
pop.aggregate(pop.pred, regions, 
    method = c('independence', 'regional'), name = method,
    inputs = list(e0F.sim.dir = NULL, e0M.sim.dir = "joint_", tfr.sim.dir = NULL), 
    verbose = FALSE)
    
get.pop.aggregation(sim.dir, name = NULL)
}

\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}} containing country-specific population projections.}
  \item{regions}{Vector of numerical codes of regions. It should correspond to values in the column \dQuote{area_code} in the \code{\link{LOCATIONS}} dataset.}
  \item{method}{There are two methods for aggregating projections, \dQuote{independence} and \dQuote{regional}, see Details. }
  \item{name}{Name of the aggregation. It becomes a part of a directory name where aggregation results are stored.}
  \item{inputs}{This argument is only used when the \dQuote{regional} method is selected. It is a list of inputs of probabilistic components of the projection:
  	\describe{
  		\item{e0F.sim.dir}{Simulation directory with projections of female life expectancy (generated using \pkg{\link{bayesLife}}). It must contain projections for the given regions (see functions \code{\link[bayesLife]{run.e0.mcmc.extra}}, \code{\link[bayesLife]{e0.predict.extra}}). If it is not given, the same e0 directory is taken which was used for generating the \code{pop.pred} object, in which case the e0 projections  are re-loaded from disk.}
  		\item{e0M.sim.dir}{Simulation directory with projections of male life expectancy. By default (value \code{NULL} or \dQuote{joint_}) the function assumes a joint female-male projections of life expectancy and thus tries to load the male projections from the female projection object created using the \code{e0F.sim.dir} argument.}
  		\item{tfr.sim.dir}{Simulation directory with projections of total fertility rate (generated using \pkg{\link{bayesTFR}}). It must contain projections for the given regions (see functions \code{\link[bayesTFR]{run.tfr.mcmc.extra}}, \code{\link[bayesTFR]{tfr.predict.extra}}). If it is not given, the same TFR directory is taken which was used for generating the \code{pop.pred} object, in which case the TFR projections  are re-loaded from disk.}
}}
  \item{verbose}{Logical switching log messages on and off.}
  \item{sim.dir}{Simulation directory where aggregation is stored. It is the same directory used for creating the \code{pop.pred} object.}
}
\details{
The dataset \code{\link{LOCATIONS}} is used to determine countries to be aggregated, in particular the field \dQuote{location_type} of the area given in the \code{regions} argument. There are three location types over which one can aggregate: Type 0 means aggregating all countries of the world, type 2 is aggregating over continents, and type 3 is aggregating over regions within continents. For type 2, countries are matched using the \dQuote{area_code} column; for type 3 the matching is done using the \dQuote{reg_code} column of the \code{\link{LOCATIONS}} dataset. 
	
There are two methods available for generating aggregations of population projection:
\describe{
\item{Independence Method}{
	Aggregations are created by summing trajectories over countries of the given region.
}
\item{Regional Method}{
The aggregation is generated using the same algorithm as population projections for single countries (function \code{\link{pop.predict}}), but it operates on aggregated input components. These are created as follows. Here \eqn{c} denotes countries over which we aggregate a region \eqn{R}, \eqn{s \in \{m, f\}}, \eqn{a}, \eqn{t} denote sex, age category and time, respectively. \eqn{t=P} denotes the present year of the prediction. \eqn{N_{s,a,t}^c} and \eqn{M_{s,a,t}^c}, respectively, denotes the historical population count and the Bayesian predictive median of population, respectively, of sex \eqn{s}, in age category \eqn{a} at time \eqn{t} for country \eqn{c} (refer to the links in parentheses for description of the data):
\describe{
	\item{Initial sex and age-specific population (\link[=PopByAgeMale]{popM, popF}):}{\eqn{N_{s,a,t=P}^R = \sum_c N_{s,a,t=P}^c}}
	\item{Sex and age-specific death rates (\link[=ASMRMale]{mxM, mxF}):}{\eqn{mx_{s,a,t}^R = \frac{\sum_c(mx_{s,a,t}^c \cdot N_{s,a,t})}{\sum_c N_{s,a,t}}}{mx_{s,a,t}^R = \sum_c(mx_{s,a,t}^c * N_{s,a,t}) / \sum_c N_{s,a,t}}}
	\item{Sex ratio at birth (\link[=SexRatioAtBirth]{srb}):}{\eqn{SRB_t^R = \frac{\sum_c M_{s=m,a=1,t}^c}{\sum_c M_{s=f,a=1,t}^c}}{SRB_t^R = \sum_c M_{s=m,a=1,t}^c / \sum_c M_{s=f,a=1,t}^c}}
	\item{Percentage age-specific fertility rate (\link[=PercentASFR]{pasfr}):}{\eqn{PASFR_{a,t}^R = \frac{\sum_c(PASFR_{a,t}^c \cdot M_{s=f,a,t})}{\sum_c M_{s=f,a,t}}}{PASFR_{a,t}^R = \sum_c(PASFR_{s=f,a,t}^c * M_{s=f,a,t}) / \sum_c M_{s=f,a,t}}}
	\item{Migration code and start year (\link[=vwBaseYear]{mig.type}):}{Aggregated migration code is the code of maximum counts over aggregated countries weighted by \eqn{N_{t=P}^c}. Migration start year is the maximum of start years over aggregated countries.}
	\item{Sex and age-specific migration (\link[=MigByAgeMale]{migM, migF}):}{\eqn{mig_{s,a,t}^R = \sum_c mig_{s,a,t}^c}}
	\item{Probabilistic projection of life expectancy:}{We assume an aggregation of life expectancy for the given regions was generated prior to this call, using the \code{\link[bayesLife]{run.e0.mcmc.extra}} and \code{\link[bayesLife]{e0.predict.extra}} functions of the \pkg{\link{bayesLife}} package.}
	\item{Probabilistic projection of total fertility rate:}{We assume an aggregation of total fertility for the given regions was generated prior to this call, using the \code{\link[bayesTFR]{run.tfr.mcmc.extra}} and \code{\link[bayesTFR]{tfr.predict.extra}} functions of the \pkg{\link{bayesTFR}} package.}	
}}
}
Aggregations can be performed for any area/region given in the \code{\link{LOCATIONS}} table. Their identifiers correspond to the column \dQuote{area_code}. In case of the \dQuote{regional} method, the same codes must be also included in predictions of life expectancy and TFR.

Results of the aggregations are stored in the same top directory as the \code{pop.pred} object, in a sudirectory called \sQuote{\code{aggregations_}\emph{name}}. They can be accessed using the function \code{get.pop.aggregation}. Note that multiple runs of this function with the same name will overwrite previous aggregations results of the same name.
}
\value{
Object of class \code{\link{bayesPop.prediction}} containing the aggregated results. In addition it contains elements \code{aggregation.method} giving the method used, and \code{aggregated.countries} which is a list of countries aggregated for each region. 
}

\author{
Hana Sevcikova, Adrian Raftery
}

\seealso{
\code{\link{pop.predict}}, \code{\link[bayesTFR]{tfr.predict.extra}}, \code{\link[bayesLife]{e0.predict.extra}}
}
\examples{
sim.dir <- tempfile()
pred <- pop.predict(countries=c(528,218,450), nr.traj = 10, output.dir=sim.dir)
aggr <- pop.aggregate(pred, 900) # aggregating World (i.e. all countries available in pred)
pop.trajectories.plot(aggr, 900, sum.over.ages=TRUE)
# countries over which we aggregated:
LOCATIONS[is.element(LOCATIONS[,'country_code'], aggr$aggregated.countries[['900']]),]
unlink(sim.dir, recursive=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }

