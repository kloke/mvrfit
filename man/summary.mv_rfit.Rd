\name{summary.mv_rfit}
\alias{summary.mv_rfit}
\title{
Summarize Multivariate Rank-Based Linear Model Fits 
}
\description{
Provides a list of rfit summaries.  
The default is to return a list of tables of coefficients.
Verbose option allows for complete rfit summary to be returned.
}
\usage{
\method{summary}{mv_rfit}(object,...)
}
\arguments{
  \item{object}{
object of class mv_rfit
}
  \item{\dots}{
addtional arguments passed to Rfit::summary.
}
}
\details{
Uses lapply to repeated call the summary funtion in Rfit on the list of fitted models.
}
\value{
List of summaries.
}
\references{
Kloke, J. and McKean, J.W. (2024), \emph{Nonparametric statistical methods using R, Second Edition}, Boca Raton, FL: Chapman-Hall.
}

\author{
John Kloke
}

\seealso{
\code{\link{summary.rfit}}, \code{\link{mv_rfit}}
}
