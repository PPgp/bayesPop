\name{pop.expressions}
\alias{pop.expressions}
\title{Expressions as used in Population Plot Functions}
\description{
	Documentation of expressions supported by functions \code{\link{pop.trajectories.plot}}, \code{\link{pop.trajectories.plotAll}}, \code{\link{pop.trajectories.table}}, \code{\link{pop.byage.plot}}, \code{\link{pop.byage.table}}, \code{\link{cohorts}}, \code{\link{pop.cohorts.plot}}, \code{\link{pop.map}}, \code{\link{pop.map.gvis}}, \code{\link{write.pop.projection.summary}}.
}
\details{
The functions above accept an argument \code{expression} which should define a population measure, i.e. a quantity that can be computed from population projections, observed population data or vital events. Such an expression is a collection of \emph{basic components} connected via usual arithmetic operators, such as \code{+},  \code{-}, \code{*}, \code{/}, \code{^}, \code{\%\%}, \code{\%/\%}, and combined using parentheses. In addition, standard R functions or predefined functions (see below) can be used within expressions.

A \bold{basic component} is a character string constituted of four parts, two of which are optional. They must be in the following order:
\enumerate{
\item Measure identification. One of the folowing upper-case characters:
	\itemize{
		\item \sQuote{P} - population, 
		\item \sQuote{D} - deaths, 
		\item \sQuote{B} - births, 
		\item \sQuote{S} - survival ratio,
		\item \sQuote{F} - fertility rate,
		\item \sQuote{R} - percent age-specific fertility,
		\item \sQuote{M} - mortality rate, 
		\item \sQuote{Q} - probability of dying,
		\item \sQuote{G} - net migration.
		}
	All but the \sQuote{P} and \sQuote{G} indicators are available only if the {\code{\link{pop.predict}}} function was run with \code{keep.vital.events=TRUE}.
\item Country part. One of the following:
	\itemize{
  		\item Numerical country code (as used in {\code{\link[wpp2012]{UNlocations}}}, see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}), 
  		\item two- or three-character ISO 3166 code, see \url{http://en.wikipedia.org/wiki/ISO_3166-1_alpha-2}, \url{http://en.wikipedia.org/wiki/ISO_3166-1_alpha-3}, 
  		\item characters \dQuote{XXX} which serves as a wildcard for a country code.
  	}
\item Sex part (optional): The country part can be followed by either \dQuote{_F} (for female) or \dQuote{_M} (for male).
\item Age part (optional): If used, the basic component is concluded by an age index given as an array. Such array is embraced by either brackets (\dQuote{[} and \dQuote{]}) or curly braces (\dQuote{\{} and \dQuote{\}}). The former invokes a summation of counts over given ages, the latter is used when no summation is desired. Note that if this part is missing, counts are automatically summed over all ages. To use all ages without summing, empty curly braces can be used. Age index one corresponds to age 0-4, index two corresponds to age 5-9 etc. Indicators \sQuote{S}, \sQuote{M} and \sQuote{Q} allow an index -1 which corresponds to age 0-1 and index 0 which corresponds to age 1-4. Use the pre-defined function \code{age.index01(\dots)} and \code{age.index05(\dots)} (see below) to define the right indices.
}
Not all combinations of the four parts above make sense. For example, \sQuote{F} and \sQuote{R} can be only combined with female sex, \sQuote{B}, \sQuote{F} and \sQuote{R} can be only combined with a subset of the age groups, namely child-bearing ages (indices 4 to 10). Or, there is no point in summing the rate indicators over sexes, i.e. using without the sex part, or over multiple age groups, i.e. using brackets.

Examples of basic components are \dQuote{P276}, \dQuote{D50_F[4:10]}, \dQuote{PXXX\{14:27\}}, \dQuote{SCZE_M\{\}}, \dQuote{QIE_M[-1]}.

When the expression is evaluated on a prediction object, each basic component is substituted by an array of four dimensions (using the \code{\link{get.pop}} function):
\enumerate{
\item Country dimension: Equals to one if a specific country code is given, or it equals the number of countries in the prediction object if a wildcard is used.
\item Age dimension: Equals to one if the third component above is missing or the age is defined within square brackets. If the age is defined within curly braces, this dimension corresponds to the length of the age array.
\item Time dimension: Depending on the time context of the expression, this dimension corresponds to either the number of projection periods or the number of observation periods.
\item Trajectory dimension: Corresponds to the number of trajectories in the prediction object, or one if the component is evaluated on observed data.
}

Depending on the context from which the expression is called, the trajectory dimension of the result of the expression can be reduced by computing given quantiles, and if only one country is evaluated, the first dimension is removed. In addition, with an exception of functions \code{\link{pop.byage.plot}}, \code{\link{pop.byage.table}},  \code{\link{cohorts}}, and \code{\link{pop.cohorts.plot}}, the expression should be constructed in a way that the age dimension is eliminated. This can be done for example by using brackets to define age, by using the \code{apply} function or one of the pre-defined functions described below. When using within \code{\link{pop.byage.plot}}, \code{\link{pop.byage.table}}, \code{\link{cohorts}}, or \code{\link{pop.cohorts.plot}}, the expression MUST include curly braces.
} 

\section{Pre-defined functions}{
\itemize{
	\item \code{gmedian(f, cat)} \cr
	It gives a median for grouped data with frequencies \code{f} and categories \code{cat}. This function is to be used in combination with \code{apply} or \code{pop.apply} (see below) along the age dimension. For example, \cr
	\dQuote{apply(P380\{\}, c(1,3,4), gmedian, cats=seq(0, by=5, length=28))} \cr
	is an expression for median age in Italy. (See \code{pop.apply} below for a simplified version.) 

\item \code{gmean(f, cat)} \cr
	Works like \code{gmedian} but gives the grouped mean.
	
\item \code{age.func(data, fun="*")} \cr
	This function applies \code{fun} to \code{data} and the corresponding age. The default case would multiply data by the corresponding age. As \code{gmedian}, it is to be used in combination with \code{apply} or \code{pop.apply}.

\item \code{drop.age(data)} \cr
	Drops the age dimension of the data. For example, if two basic components are combined where one is used within the \code{apply} function, the other will need to change its dimension in order to have conformable arrays. For example, \cr
	\dQuote{apply(age.func(P752\{\}), c(1,3,4), sum) / drop.age(P752)} \cr
	is an expression for the average age in Sweden. (See \code{pop.apply} below for a simplified version.)

\item \code{pop.apply(data, fun, \dots, split.along=c("None", "age", "traj", "country"))} \cr
	By default applies function \code{fun} to the age dimension of \code{data} and converts the result into the same format as returned by a basic component. This allows combining the \code{apply} function with other basic components without having to modify their dimensions. For example, \cr
	\dQuote{pop.apply(age.func(P752\{\}), fun=sum) / P752} gives the average age in Sweden, or \cr
	\dQuote{pop.apply(P380\{\}, gmedian, cats=seq(0, by=5, length=28))} gives the median age of Italy.
	If \code{slice.along} is not \sQuote{None}, it can be used as an \code{apply} function where the data is sliced along one axis. 
\item \code{pop.combine(data1, data2, fun, ..., split.along=c("age", "traj", "country"))} \cr
	Can be used if two basic components should be combined that result in different shapes. It tries to put data into the right format and calls \code{pop.apply}. For example, \cr
	\dQuote{pop.combine(PIND\{\}, PIND, '/')} give population by age per total population in India, or \cr
	\dQuote{pop.combine(BFR - DFR, GFR, '+', split.along='traj')} gives births minus deaths plus net migration in France. Here, \code{pop.combine} is necessary, because \sQuote{GFR} is a deterministic component and thus, has only one trajectory, whereas births and deaths are probabilistic.
\item \code{age.index01(end)} \cr
	Can be used with indicators \sQuote{S}, \sQuote{M} and \sQuote{Q} only. It returns an array of age group indices that include ages 0-1 and 1-4 and exclude 0-4. The last age index is \code{end}.
\item \code{age.index05(end)} \cr
	Returns an array of age group indices starting with group 0-4, 5-9 until the age group corresponding to index \code{end}.
}
There is also a help function available that generates an expression for the mean age of childbearing, see \code{\link{mac.expression}}.
}

\note{The expression parser is simple and far from being perfect. We recommend to leave spaces around the basic components.}

\author{Hana Sevcikova, Adrian Raftery}
\seealso{\code{\link{mac.expression}}, \code{\link{get.pop}}, \code{\link{pop.trajectories.plot}}, \code{\link{pop.map}}, \code{\link{write.pop.projection.summary}}.}
\examples{
sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
pred <- get.pop.prediction(sim.dir, write.to.cache=FALSE)

# median age of women in child-bearing ages in Netherlands and all countries - trajectories
pop.trajectories.plot(pred, nr.traj=50,
    expression="pop.apply(P528_F{4:10}, gmedian, cats= seq(15, by=5, length=8))")
pop.trajectories.plotAll(pred, nr.traj=50, 
    expression="pop.apply(PXXX_F{4:10}, gmedian, cats= seq(15, by=5, length=8))")

# mean age of women in child-bearing ages in Netherlands - table
pop.trajectories.table(pred, 
    expression="pop.apply(age.func(P528_F{4:10}), fun=sum) / P528_F[4:10]")
# for the mean age of childbearing, see ?mac.expression

# migration per capita by age
pop.byage.plot(pred, expression="GNL{} / PNL{}", year=2000)

\dontrun{
# potential support ratio - map
pop.map(pred, expression="PXXX[5:13] / PXXX[14:27]")}

# proportion of 0-4 years old to whole population - export to an ASCII file
dir <- tempfile()
write.pop.projection.summary(pred, expression="PXXX[1] / PXXX", output.dir=dir)
unlink(dir)

\dontrun{
# These are vital events only available if keep.vital.events=TRUE in pop.predict 
# log female mortality rate by age for Netherlands in 2050, including 0-1 and 1-4 age groups
pop.byage.plot(pred, expression="log(MNL_F{age.index01(27)})", year=2050)

# trajectories of male 1q0 and table of 5q0 for Netherlands
pop.trajectories.plot(pred, expression="QNLD_M[-1]")
pop.trajectories.table(pred, expression="QNLD_M[1]")}
}
\keyword{documentation}
