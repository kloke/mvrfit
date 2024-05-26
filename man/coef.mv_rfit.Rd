\name{coef.mv_rfit}
\alias{coef.mv_rfit}
\title{Return Multivariate Regression Coefficients}
\description{
Return Multivariate Regression Coefficients
}
\usage{
\method{coef}{mv_rfit}(object,...)
}
\arguments{
  \item{object}{ object of class mv_rfit }
  \item{\dots}{ additional arguments. see description below.}
}
\details{
additional arguments are passed to coef.
as_list (logical) may be used to return the result as a list rather than as a matrix (the default).  
}
\value{
Returns list or matrix of the fitted multivariate regression model coefficients.
}
\author{ John Kloke }

\seealso{
\code{\link{mv_rfit}}
}
