\name{paspalum}
\alias{paspalum}
\docType{data}
\title{
Study of the Plant Paspalum }
\description{
Design is A cross B with A at 2 levels ((1) Control and (2) Inoculated with a fungal infection) and B at 4 levels of temperatures ((1) 14 degress C; (2) 18 degrees C.; (3) 22 degrees C.; and (4) 26 degrees C.).  
For each experimental unit, 3 observations were taken: X1 the fresh weight of roots in gm; X2 the maximum root length in mm; and X3 the fresh weight of tops in gm).
}
\usage{data("paspalum")}
\format{
	A data frame with 48 observations and 5 variables.
	\describe{
		\item{\code{X1}}{Fresh weight of roots in gm}
		\item{\code{X2}}{Maximum root length in mm}
		\item{\code{X3}}{Fresh weight of tops in gm}
		\item{\code{c1}}{Indicator for Factor A: (1) Control and (2) Inoculated}
		\item{\code{c2}}{Indicator for Factor B: Inoculation applied with temperatures (1) 14 degress C; (2) 18 degrees C.; (3) 22 degrees C.; and (4) 26 degrees C.}
		}
		}
\details{
Seber (1984) presents the results of a study on the growth of paspalum grass.
The design was A cross B with A at 2 levels ((1) Control and (2) Inoculated with a fungal infection) and B at 4 levels of temperatures ((1) 14 degress C; (2) 18 degrees C.; (3) 22 degrees C.; and (4) 26 degrees C.).  For each experimental unit, 3 observations were taken: X1 the fresh weight of roots in gm; X2 the maximum root length in mm; and X3 the fresh weight of tops in gm).
}
\source{
	Seber, G.A.F., (1984), \emph{Multivariate Observations}, John Wiley: NY.
}
\references{
Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.
}
\examples{
data(paspalum)
boxplot(X1~c1,data=paspalum)
}
\keyword{datasets}
