\name{twoway.mvrfit}
\alias{twoway.mvrfit}
\title{ Rank-based Multivariate Two-way Analysis of Variance }
\description{
Multivariate Two-way Analysis of Variance based on a rank-based fit.
}
\usage{
twoway.mvrfit(y, g, ...)
}
\arguments{
  \item{y}{matrix of responses.}
  \item{g}{n x 2 matrix of group membership}
  \item{\dots}{additional arguments. passed to fitting routine.}
}
\details{
Performs a multivariate two-way ANOVA as described in Kloke and McKean (2024).
}
\value{
\item{anova.table}{table of tests for interaction and main effects}
\item{fit}{results from call to mv_rfit}
\item{muhat}{estimated cell medians}
}
\references{
Kloke, J. and McKean, J.W. (2024), \emph{Nonparametric statistical methods using R, Second Edition}, Boca Raton, FL: Chapman-Hall.
}
\author{
John Kloke, Joseph McKean
}
\seealso{
\code{\link{mv_rfit}}, ~~~
}
