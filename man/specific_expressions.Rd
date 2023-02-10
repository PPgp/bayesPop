\name{mac.expression}
\alias{mac.expression}
\alias{mac.expression1}
\alias{mac.expression5}
\title{
Expression Generator
}
\description{
Help functions to easily generate commonly used expressions. 
}
\usage{
mac.expression(country)
mac.expression1(country)
mac.expression5(country)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{country}{Country code as defined for \code{\link[=pop.expressions]{expressions}}.
}
}
\details{
\code{mac.expression} and \code{mac.expression1} generate expressions for the mean age of childbearing of the given country, for 5-year age groups and 1-year age groups, respectively. \code{mac.expression5} is a synonym for \code{mac.expression}.
Note that \code{\link{pop.predict}} has to be run with \code{keep.vital.events=TRUE} for this to work.
}
\value{
\code{mac.expression} returns a character string corresponding to the formula
{\eqn{(17.5*R_c(15-19) + 22.5*R_c(20-24) + ... + 47.5*R_c(45-49))/100}}
where \eqn{R_c(x)} denotes the country-specific percent age-specific fertility for the age group \eqn{x}. 

\code{mac.expression1} returns a character string corresponding to the formula
{\eqn{(10.5*R_c(10-11) + 11.5*R_c(11-12) + ... + 54.5*R_c(54-55))/100}}
}

\seealso{
\code{\link{pop.expressions}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
# Run pop.predict with storing vital events
pred <- pop.predict(countries=c("Germany", "France"), nr.traj=3, 
           keep.vital.events=TRUE, output.dir=sim.dir)
# plot the mean age of childbearing 
pop.trajectories.plot(pred, expression=mac.expression("FR"), cex.main = 0.7)
unlink(sim.dir, recursive=TRUE)}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
