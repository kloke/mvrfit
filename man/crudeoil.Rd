\name{crudeoil}
\alias{crudeoil}
\docType{data}
\title{
Study of Crude Oil versus Zones of Sandstone}
\description{
Samples of crude oil were taken from three zones of sandstone: (1) Wilhelm; (2) Sub-mulina; and (3) Upper.   For each sampled item 5 variables were measured: vanadium; iron; beryllium; sathydro; and aromhydro.
}
\usage{data("crudeoil")}
\format{
	A data frame consisting of 56 obervations and 6 variables.
	\describe{
		\item{\code{vanadium}}{ Vanadium in percent ash}
		\item{\code{iron}}{ Iron in percent ash}
		\item{\code{beryllium}}{ Beryllium in percent ash}
		\item{\code{sathydro}}{ Saturated hydrocarbons in percent area}
		\item{\code{aromhydro}}{ Aromatic hydrocarbons in percent area}
		\item{\code{zone}}{ Indicator variable: (1) Wilhelm; (2) Sub-mulina; and (3) Upper}
	}
}
\details{
Johnson and Wichern (2007) discuss a study presented by Gerrild and Lantz (1969).  Samples of crude oil were taken from three zones of sandstone: (1) Wilhelm; (2) Sub-mulina; and (3) Upper.   For each sampled item 5 variables were measured: vanadium; iron; beryllium; sathydro; and aromhydro.
}
\source{
Gerrild, P.M. and Lantz, R.J. (1969), Chemical analysis of 75 crude oil samples from Pliocene sand units, Elk Hills Oil Field, California, \emph{U.S. Geological Survey Open-file Report}.
}
\references{
Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.

Johnson, R.A. and Wichern, D.W., (2007), \emph{Applied multivariate statistical analysis, 6th ed.}, Upper Saddle River, NJ: Pearson.
}
\examples{
data(crudeoil)
boxplot(vanadium~zone,data=crudeoil)
}
\keyword{datasets}
