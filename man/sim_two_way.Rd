\name{sim_two_way}
\alias{sim_two_way}
\docType{data}
\title{
Simulated dataset for a multivariate two-way design
}
\description{
Simulated dataset for a multivariate two-way design.
}
\usage{data("sim_two_way")}
\format{
  A data frame with 30 observations on the following 5 variables.
  \describe{
    \item{\code{y1}}{a numeric vector. response 1.}
    \item{\code{y2}}{a numeric vector. response 2.}
    \item{\code{y3}}{a numeric vector. response 3.}
    \item{\code{f1}}{a factor with levels \code{1} \code{2}}
    \item{\code{f2}}{a factor with levels \code{1} \code{2} \code{3}}
  }
}
\details{
Generated data for a 2x3 design. 
The responses are vectors of three measurements. The design is balanced with 5 observations in each of the 6 cells. Some of the true interactions are not zero in the model that generated the responses. The random errors were drawn from a multivariate t-distribution with 5 degrees of freedom.
}
\source{
Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.
}
\references{
Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.
}
\examples{
data(sim_two_way)
}
\keyword{datasets}
