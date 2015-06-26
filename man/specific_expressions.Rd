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
  \item{country}{Country code as defined for \code{\link[expressions]{pop.expressions}}.
}
}
\details{
\code{mac.expression} generates an expression for the mean age of childbearing of the given country.
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
    sim.dir <- file.path(find.package("bayesPop"), "ex-data", "Pop")
    pred <- get.pop.prediction(sim.dir)
    # plot the mean age of childbearing 
    pop.trajectories.plot(pred, expression=mac.expression("NL"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
