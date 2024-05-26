\name{mv_rfit}
\alias{mv_rfit}
\title{
Multivariate Rank-based Regression
}
\description{
An extention of Rfit for the multivariate (multiple responses per subject) case.
}
\usage{
mv_rfit(x, y,...)
}
\arguments{
  \item{x}{
n x p design matrix
}
  \item{y}{
n x k response matrix
}
  \item{\dots}{
additional arguments passed to rfit
}
}
\details{
Fits model 
Y = 1 alpha' + X Beta + E 
where
Y is an nxk matrix where the rows represent response vectors for subjects 
X is an nxp matrix design matrix
Beta is a pxk matrix of regression coefficients
E is an nxk matrix of error terms where the rows are assumed to be independent
1 is an nx1 vector on ones and alpha is px1 vector of intercept terms.

Each column of Y is fit versus the columns of X via calls to rfit.
}
\value{
\item{coefficients}{List of length k containg (p+1)x1 vectors of coefficients}
\item{residuals}{nxk matrix of residuals}
\item{fitted.values}{nxk matrix of fitted values}
\item{scores}{object of class scores used for the fit}
\item{x}{nxp design matrix}
\item{y}{nxk response matrix}
\item{tauhatVec}{kx1 vector of tauhats}
\item{taushatVec}{kx1 vector of taushats}
\item{dispVec}{kx1 vector of dispersion function evaluated at the full model fit}
\item{S}{kxk matrix of correlation of the scores}
\item{fits}{List of length k of the list returned by rfit for each of the fits}
}
\references{
Davis, J. and McKean. (1993), Rank-based methods for multivariate linear models," \emph{Journal of the American Statistical Association}, 88, 245-251.

Hettmansperger, T. P. and McKean, J. W. (2011), \emph{Robust nonparametric statistical methods, 2nd ed.}, Boca Raton, FL: Chapman & Hall.

Johnson, R. and Wichern, D. (2007), \emph{Applied multivariate statistical analysis, 6th ed.} Upper Saddle River, NJ: Pearson.

Kloke, J. and McKean, J.W. (2024), \emph{Nonparametric statistical methods using R, Second Edition}, Boca Raton, FL: Chapman-Hall.

}
\author{
John Kloke, Joseph McKean
}

\seealso{
\code{\link{rfit}}
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
\keyword{regression}
