\name{lincomb.mv_rfit}
\alias{lincomb.mv_rfit}
\title{
Estimate and Confidence Interval for Linear Combination
}
\description{
Estimate and Confidence Interval for Linear Combination in a Multivariate Regression Setting.
}
\usage{
lincomb.mv_rfit(fit, h, l, conf.level = 0.95)
}
\arguments{
  \item{fit}{result from mv_rfit }
  \item{h}{vector} 
  \item{l}{vector}
  \item{conf.level}{scalar.confidence level}
}
\details{
Estimates confidence interval of the form h^T beta l
as described in Chapter 10 of Kloke and McKean (2024).
}
\value{
Table with estimate (Est), standard error (SE), lower confidence limit (LCL), and upper confidence limit (LCL) is returned.
}
\references{
Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.
}
\author{
John Kloke, Joseph McKean
}
\seealso{
\code{\link{mv_rfit}}
}
