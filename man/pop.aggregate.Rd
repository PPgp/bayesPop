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
    input.type = c("country", "region"), name = input.type,
    inputs = list(e0F.sim.dir = NULL, e0M.sim.dir = "joint_", tfr.sim.dir = NULL),
    my.location.file = NULL, verbose = FALSE, \dots)
    
get.pop.aggregation(sim.dir = NULL, pop.pred = NULL, name = NULL, 
    write.to.cache = TRUE)
}

\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}} containing country-specific population projections.}
  \item{regions}{Vector of numerical codes of regions. It should correspond to values in the column \dQuote{country_code} in the \code{\link[wpp2012]{UNlocations}} dataset or in \code{my.location.file} (see below).}
  \item{input.type}{There are two methods for aggregating projections depending on the type of inputs, \dQuote{country}- and \dQuote{region}-based, see Details. }
  \item{name}{Name of the aggregation. It becomes a part of a directory name where aggregation results are stored.}
  \item{inputs}{This argument is only used when the \dQuote{region}-based method is selected. It is a list of inputs of probabilistic components of the projection:
  	\describe{
  		\item{e0F.sim.dir}{Simulation directory with projections of female life expectancy (generated using \pkg{\link{bayesLife}}). It must contain projections for the given regions (see functions \code{\link[bayesLife]{run.e0.mcmc.extra}}, \code{\link[bayesLife]{e0.predict.extra}}). If it is not given, the same e0 directory is taken which was used for generating the \code{pop.pred} object, in which case the e0 projections  are re-loaded from disk.}
  		\item{e0M.sim.dir}{Simulation directory with projections of male life expectancy. By default (value \code{NULL} or \dQuote{joint_}) the function assumes a joint female-male projections of life expectancy and thus tries to load the male projections from the female projection object created using the \code{e0F.sim.dir} argument.}
  		\item{tfr.sim.dir}{Simulation directory with projections of total fertility rate (generated using \pkg{\link{bayesTFR}}). It must contain projections for the given regions (see functions \code{\link[bayesTFR]{run.tfr.mcmc.extra}}, \code{\link[bayesTFR]{tfr.predict.extra}}). If it is not given, the same TFR directory is taken which was used for generating the \code{pop.pred} object, in which case the TFR projections  are re-loaded from disk.}
}}
\item{my.location.file}{User-defined location file that can contain other agreggation groups than the default UN location file. It should have the same structure as the \code{\link[wpp2012]{UNlocations}} dataset, see below.}
\item{verbose}{Logical switching log messages on and off.}
\item{sim.dir}{Simulation directory where aggregation is stored. It is the same directory used for creating the \code{pop.pred} object. Alternatively, \code{pop.pred} can be used. Either \code{sim.dir} or \code{pop.pred} must be given.}
\item{write.to.cache}{Logical controlling if functions operating on this object are allowed to write into its cache (see Details of \code{\link{get.pop.prediction}}).}
\item{\dots}{Additional arguments.}
}
\details{
The dataset \code{\link[wpp2012]{UNlocations}} or \code{my.location.file} is used to determine countries to be aggregated, in particular the field \dQuote{location_type} of the entries with \dQuote{country_code} given in the \code{regions} argument. One can aggregate over the following location types: Type 0 means aggregating all countries of the world (or in the file), type 2 is aggregating over continents, type 3 is aggregating over regions within continents, and any other integer (except 4) correponds to user-defined aggregations. Note that type 4 is reserved as a location type of countries and thus, all aggregations are performed over entries of this type. For type 2, countries are matched using the \dQuote{area_code} column; for type 3 the matching is done using the \dQuote{reg_code} column of the \code{\link[wpp2012]{UNlocations}} dataset. E.g., if \code{regions=908} (Europe) which has location type 2 in the default \code{\link[wpp2012]{UNlocations}} dataset, all countries are aggregated for which values of 908 are found in the \dQuote{area_code} column. If the location type is other than 0, 2, 3 and 4, there must be a column in the file called \dQuote{agcode_\eqn{x}} with \eqn{x} being the location type. This column is then used to match the countries to be aggregated. 

Consider the following example. Say we want to pair four countries (Germany [DE], France [FR], Netherlands [NL], Italy [IT]) in two different ways, so we have two overlapping grouppings, each of which has two groups (A,B): 
\enumerate{
	\item{group A = (DE, FR), group B = (NL, IT)}
	\item{group A = (DE, NL), group B = (FR, IT)}
}
Then,  \code{my.location.file} should have the following entries:

\tabular{rlrrr}{
	country_code \tab name \tab location_type \tab agcode_98 \tab agcode_99 \cr
	1001  \tab		groupping1_groupA \tab	98	\tab	-1 \tab -1 \cr			
	1002	  \tab		groupping1_groupB \tab	98	\tab	-1 \tab -1 \cr
	1003	  \tab		groupping2_groupA \tab	99  \tab	-1 \tab -1 \cr
	1004	 \tab		groupping2_groupB \tab	99 \tab	-1 \tab -1 \cr
	276	\tab	Germany	\tab		4		\tab	1001	 \tab 1003 \cr
	250	\tab	France	\tab		4	\tab		1001	\tab 1004 \cr
	258	\tab	Netherlands	 \tab	4	\tab	1002		\tab	1003 \cr
	380	\tab	Italy		\tab		4	\tab		1002	 \tab 1004 \cr
	1005 \tab		all \tab 0 \tab	-1 \tab -1 
}
The \dQuote{country_code} of the groups is user-specific, but it must be unique within the file. Values of \dQuote{country_code} for countries must match those in the prediction object. To run the aggregation for the four groups above we set \code{regions=1001:1004}. Having \dQuote{location_type} being 98 and 99, it is expected the file to have columns  \dQuote{agcode_98} and \dQuote{agcode_99} containing assignements to each of the two grouppings. Values in this columns corresponding to groups are not used and thus can have any value. For aggregating over all four countries, set \code{regions=1005} which has \dQuote{location_type} equal 0 and thus, it is aggregated over all entries with \dQuote{location_type} equals 4.


There are two methods available for generating aggregations of population projection:
\describe{
\item{Country-based Method}{
	Aggregations are created by summing trajectories over countries of the given region.
}
\item{Region-based Method}{
The aggregation is generated using the same algorithm as population projections for single countries (function \code{\link{pop.predict}}), but it operates on aggregated input components. These are created as follows. Here \eqn{c} denotes countries over which we aggregate a region \eqn{R}, \eqn{s \in \{m, f\}}{s in {m,f}}, \eqn{a}, and \eqn{t} denote sex, age category and time, respectively. \eqn{t=P} denotes the present year of the prediction. \eqn{N_{s,a,t}^c} and \eqn{M_{s,a,t}^c}, respectively, denotes the historical population count and the Bayesian predictive median of population, respectively, of sex \eqn{s}, in age category \eqn{a} at time \eqn{t} for country \eqn{c} (refer to the links in parentheses for description of the data):
\describe{
	\item{Initial sex and age-specific population (\link[wpp2012:popM]{popM, popF}):}{\eqn{N_{s,a,t=P}^R = \sum_c N_{s,a,t=P}^c}}
	\item{Sex and age-specific death rates (\link[wpp2012:mxM]{mxM, mxF}):}{\eqn{mx_{s,a,t}^R = \frac{\sum_c(mx_{s,a,t}^c \cdot N_{s,a,t})}{\sum_c N_{s,a,t}}}{mx_{s,a,t}^R = \sum_c(mx_{s,a,t}^c * N_{s,a,t}) / \sum_c N_{s,a,t}}}
	\item{Sex ratio at birth (\link[wpp2012:sexRatio]{srb}):}{\eqn{SRB_t^R = \frac{\sum_c M_{s=m,a=1,t}^c}{\sum_c M_{s=f,a=1,t}^c}}{SRB_t^R = \sum_c M_{s=m,a=1,t}^c / \sum_c M_{s=f,a=1,t}^c}}
	\item{Percentage age-specific fertility rate (\link[wpp2012:percentASFR]{pasfr}):}{\eqn{PASFR_{a,t}^R = \frac{\sum_c(PASFR_{a,t}^c \cdot M_{s=f,a,t})}{\sum_c M_{s=f,a,t}}}{PASFR_{a,t}^R = \sum_c(PASFR_{s=f,a,t}^c * M_{s=f,a,t}) / \sum_c M_{s=f,a,t}}}
	\item{Migration code and start year (\link[=vwBaseYear]{mig.type}):}{Aggregated migration code is the code of maximum counts over aggregated countries weighted by \eqn{N_{t=P}^c}. Migration start year is the maximum of start years over aggregated countries.}
	\item{Sex and age-specific migration (\link[wpp2012:migrationM]{migM, migF}):}{\eqn{mig_{s,a,t}^R = \sum_c mig_{s,a,t}^c}}
	\item{Probabilistic projection of life expectancy:}{We assume an aggregation of life expectancy for the given regions was generated prior to this call, using the \code{\link[bayesLife]{run.e0.mcmc.extra}} and \code{\link[bayesLife]{e0.predict.extra}} functions of the \pkg{\link{bayesLife}} package.}
	\item{Probabilistic projection of total fertility rate:}{We assume an aggregation of total fertility for the given regions was generated prior to this call, using the \code{\link[bayesTFR]{run.tfr.mcmc.extra}} and \code{\link[bayesTFR]{tfr.predict.extra}} functions of the \pkg{\link{bayesTFR}} package.}	
}}
}

Results of the aggregations are stored in the same top directory as the \code{pop.pred} object, in a sudirectory called \sQuote{\code{aggregations_}\emph{name}}. They can be accessed using the function \code{get.pop.aggregation}. Note that multiple runs of this function with the same name will overwrite previous aggregations results of the same name.
}
\value{
Object of class \code{\link{bayesPop.prediction}} containing the aggregated results. In addition it contains elements \code{aggregation.method} giving the \code{input.type} used, and \code{aggregated.countries} which is a list of countries aggregated for each region. 
}

\references{
H. Sevcikova, A. E. Raftery (2016). bayesPop: Probabilistic
  Population Projections. Journal of Statistical Software, 75(5), 1-29.
  doi:10.18637/jss.v075.i05
}

\author{
Hana Sevcikova, Adrian Raftery
}

\seealso{
\code{\link{pop.predict}}, \code{\link[bayesTFR]{tfr.predict.extra}}, \code{\link[bayesLife]{e0.predict.extra}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
pred <- pop.predict(countries=c(528,218,450), output.dir=sim.dir)
aggr <- pop.aggregate(pred, 900) # aggregating World (i.e. all countries available in pred)
pop.trajectories.plot(aggr, 900, sum.over.ages=TRUE)
# countries over which we aggregated:
subset(UNlocations, country_code \%in\% aggr$aggregated.countries[["900"]])
unlink(sim.dir, recursive=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }

