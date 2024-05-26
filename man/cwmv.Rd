\name{cwmv}
\alias{cwmv}
\title{
Componentwise Multivariate
}
\description{
Componentwise multivariate analysis including rank-based and least squares methods.
}
\usage{
cwmv(x, method = "HL", H = diag(rep(1, ncol(x))), g0 = rep(0, ncol(x)), C = diag(rep(1, ncol(x))), test = FALSE, test.method = "Wald", conf.int = FALSE, conf.level = 0.95, conf.adjust = "none", ...)
}
\arguments{
  \item{x}{n x k matrix of responses}
  \item{method}{Estimation method. Default is Hodges-Lehmann (HL)}
  \item{H}{k x q dimensional hypothesis matrix (for testing)}
  \item{g0}{null value of Hypothesis}
  \item{C}{contrast matrix for estimation}
  \item{test}{logical. whether to perform the test of hypothesis. default is FALSE.}
  \item{test.method}{character scalar. which method of test to use.}
  \item{conf.int}{logical. whether to estimate the confidence intervals.  default is TRUE.}
  \item{conf.level}{scalar numeric. level confidence to use.  default is 0.95}
  \item{conf.adjust}{character scalar. }
  \item{\dots}{addtional arguments.  currently not used.}
}
\details{
Performs componentwise multivariate analysis on the columns of the data matrix x.  Options for method include Wilcoxon or Hodges-Lehmann (HL); sign or median (Median); mean or least squares (LS).  Confidence intervals for contrasts are provided by default (option to omit).  Optional Bonferroni (bonferroni) adjustment.
Tests of hypothesis my be requested including Wald and Gradient tests.
}
\value{
\item{F.statistic}{F statistic for the test of hypothesis}
\item{p.value}{p-value for the test of hypothesis}
\item{g0}{null value used in the test of hypothesis}
\item{test.method}{which method was used to perform the test of hypothesis}
\item{conf.int.tab}{matrix (table) of confidence intervals}
\item{conf.adjust}{method of confidence interval adjustment}

}
\references{
Hettmansperger, T. P. and McKean, J. W (2011), Robust nonparametric statistical methods, 2nd Edition, Boca Raton, FL: Chapman & Hall.

Kloke, J. and McKean, J.W. (2024), \emph{Nonparametric statistical methods using R, Second Edition}, Boca Raton, FL: Chapman-Hall.
}
\author{
John Kloke, Joseph McKean
}

\examples{
n <- 100
xmat <- cbind(rnorm(n,pi),rnorm(n,2*pi),rnorm(n,4*pi))
cwmv(xmat,method='HL',conf.int=TRUE)
cwmv(xmat,method='Median',conf.int=TRUE)
cwmv(xmat,method='LS',conf.int=TRUE)
}
