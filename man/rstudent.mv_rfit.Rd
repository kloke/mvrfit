\name{rstudent.mv_rfit}
\alias{rstudent.mv_rfit}
\title{
Studentized residuals for a multivariate rank-based fit.
}
\description{
Studentized residuals for a multivariate rank-based fit.
}
\usage{
\method{rstudent}{mv_rfit}(model, ...)
}
\arguments{
  \item{model}{result from a mv_rfit.}
  \item{\dots}{additional arguments.  see details below.}
}
\details{
Studentized residuals from the individual rank-based fits for each of the responses are returned.
Default is to return as a matrix.  Option to return as a list using as_list set to TRUE.
}
\value{
Matrix or list of residuals.
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.

Kloke, J. and McKean, J.W. (2024), \emph{Nonparametric statistical methods using R, Second Edition}, Boca Raton, FL: Chapman-Hall.
}
\author{
John Kloke, Joseph McKean
}
