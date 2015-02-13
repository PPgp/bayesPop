\name{MLTbx}
\alias{MLTbx}
\docType{data}
\title{
Dataset on Lee-Carter bx for Modeled Countries
}
\description{
Dataset with values of the Lee-Carter bx parameter for countries where mortality was obtained using model life tables.
}
\usage{
    data(MLTbx)
}
\format{
  A data frame with nine rows and 28 columns. Each row corresponds to one mortality age pattern as defined in the \code{\link{vwBaseYear}} dataset. Each column corresponds to an age group, starting with 0-1, 1-4, 5-9, 10-14, ... up to 125-129, 130+.}
\details{
	These values are used for countries for which the column \code{AgeMortalityType} in \code{\link{vwBaseYear}} is equal to \dQuote{Model life tables}. In such a case a row is selected that corresponds to the corresponding value of the column \code{AgeMortalityPattern} (also in  \code{\link{vwBaseYear}}). These values are then used instead of estimating the Lee-Carter \eqn{b_x} from the country's historical data.
}
\source{
Data provided by the United Nations Population Division.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\seealso{\code{\link{vwBaseYear}}}

\examples{
data(MLTbx)
str(MLTbx)
}
\keyword{datasets}
