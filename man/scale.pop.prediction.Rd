\name{pop.scale.prediction}
\alias{pop.scale.prediction}
\alias{write.scaled.pop}

\title{
Population Scaling To Match Aggregated Targets
}
\description{
Functions scale given sex- and age-specific population counts for each location and year, so that the sum of the scaled counts over all locations matches the given target for each year, sex and age. Function \code{scale.pop.prediction} creates a new prediction object that contain the scaled trajectories. Function \code{write.scaled.pop} writes a file containing the scaled statistics for each location, year, sex and age.
}
\usage{
pop.scale.prediction(pop.pred, target.file, output.dir,
    target.code = NULL, variant.name = "mean", 
    target.id.column = "country_code", stat = "mean", 
    exclude.codes = NULL, verbose = TRUE)
    
write.scaled.pop(pop.pred, target.file, 
    output.file = "adjusted_population.txt", target.code = NULL, 
    variant.name = "mean", target.id.column = "country_code", 
    output.id.column = "reg_code", stat = "mean", 
    exclude.codes = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop.pred}{Object of class \code{\link{bayesPop.prediction}} containing population projections to be scaled.}
  \item{target.file}{Comma-delimited file containing the summary statistics of the aggregated targets. For example, if locations in \code{pop.pred} are states, this file would contain the summary statistics for the corresponding country. 
  It is expected to have the same format as a file created with \code{\link{write.pop.projection.summary}(what = 'popsexage')}. That is, it should have the following columns: \code{country_code} (unique identifier of the aggregated locations), \code{sex} (either \dQuote{male} or \dQuote{female}), \code{age} (values from 0 to 130), columns for each year or time period.  Optionally, column \code{variant} can be present.
}
  \item{output.dir}{Directory to write the scaled projection into.}
  \item{target.code}{If the \code{target.file} contains more than one aggregated location, this argument gives the code of the location to use.
  }
  \item{variant.name}{If the \code{target.file} contains the column \code{variant}, this argument gives the value of that column for selecting the summary statistics. If the \code{target.file} is created via the \code{\link{write.pop.projection.summary}} function, run it with the argument \code{include.means = TRUE} to generate a file with the \dQuote{mean} variant, which is the default value here.
}
  \item{target.id.column}{Name of the column in \code{target.file} that determines the unique identifiers of the aggregated locations.
}
  \item{stat}{Summary statistics to use for the scaling. By default, the scaling shifts are determined by comparing the means of the trajectories with the targets. 
}
 \item{exclude.codes}{Vector of location codes that should be excluded from the scaling.
 }
  \item{output.file}{File name to write the scaled sumary statistics into.
}
  \item{output.id.column}{Column name of the unique identifier of the locations in the output file. 
}
\item{verbose}{Logical switching processing info on and off.}
}
\details{
These functions can be used if it is necessary to scale the population distribution so that a summary statistics of the sum of the locations matches some given targets. For example, if states should sum up to a given country projection. Thus, its use is meant to be mainly for subnational projections.

The scaling is done for each sex and age, and the scaling proportions are determined by the population sizes within each sex and age group. 

The function \code{scale.pop.prediction} creates a new prediction object, saved into the directory given by \code{output.dir}. It contains adjusted population trajectories, as well as recomputed quantiles. Note however that only population quantities are affected by this scaling. There is no change to vital events or rates. 

Function \code{write.scaled.pop} derives the necessary shifts for the scaling and stores the adjusted summary statistics into a file. It does not create any new projection object.
}
\value{
Function \code{scale.pop.prediction} returns an object of class \code{\link{bayesPop.prediction}} containing the scaled predictions. Function \code{write.scaled.pop} returns the name of the file that it created.
}

\author{
Hana Sevcikova
}

\keyword{ manip }
\keyword{ print }

