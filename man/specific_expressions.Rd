\name{mac.expression}
\alias{mac.expression}
\title{
Expression Generator
}
\description{
Help functions to easily generate commonly used expressions. 
}
\usage{
mac.expression(country)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{country}{Country code as defined for \code{\link[=pop.expressions]{expressions}}.
}
}
\details{
\code{mac.expression} generates an expression for the mean age of childbearing of the given country. Note that \code{\link{pop.predict}} has to be run with \code{keep.vital.events=TRUE} for this to work.
}
\value{
\code{mac.expression} returns a character string corresponding to the formula
{\eqn{(17.5*R_c(15-19) + 22.5*R_c(20=24) + ... + 47.5*R_c(45-49))/100}}
where \eqn{R_c(x)} denotes the country-specific percent age-specific fertility for the age group \eqn{x}. 
}

\seealso{
\code{\link{pop.expressions}}
}
\examples{
# Run pop.predict with storing vital events
\dontrun{
sim.dir <- tempfile()
pred <- pop.predict(countries=c("Germany", "France"), nr.traj=3, 
           keep.vital.events=TRUE, output.dir=sim.dir)
# plot the mean age of childbearing 
pop.trajectories.plot(pred, expression=mac.expression("FR"))
unlink(sim.dir, recursive=TRUE)}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
