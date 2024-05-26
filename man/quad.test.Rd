\name{quad.test}
\alias{quad.test}
\title{
Quadradic Test of Hypothesis in a Multivariate Regression Setting.
}
\description{
Quadradic Test of Hypothesis in a Multivariate Regression Setting.
Similar to the traditional Lawley-Hotelling developed for a robust rank-based analysis.
}
\usage{
quad.test(fit, M, K, ...)
}
\arguments{
  \item{fit}{result from mv_rfit }
  \item{M}{matrix.}
  \item{K}{matrix.}
  \item{\dots}{addtional arguments. currently not used.}
}
\details{
Test the hypothesis H0: M beta K = 0 versus HA: M beta K ne 0.
Test is described in Chapter 10 of Kloke and McKean 2024.
This test and its asymptotic theory were developed by Davis and McKean
(1993); see, also, Chapter 6 of Hettmansperger and McKean (2011).
}
\value{
\item{test.statistic}{test statistic}
\item{p.value}{p-value based on chi-squared distribution}
\item{M}{M matrix}
\item{K}{K matrix}
}
\references{
Davis, J. and McKean. (1993), Rank-based methods for multivariate linear models," \emph{Journal of the American Statistical Association}, 88, 245-251.

Hettmansperger, T. P. and McKean, J. W. (2011), \emph{Robust nonparametric statistical methods, 2nd ed.}, Boca Raton, FL: Chapman & Hall.

Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.
}
\author{
John Kloke, Joseph McKean
}
\seealso{
\code{\link{mv_rfit}}
}

