\name{anteaters}
\alias{anteaters}
\docType{data}
\title{Study of the Size of Anteaters}
\description{
Results from a study on anteaters concerning their skull size from three different areas of South America.
Three measurements (basal length (y1), occipitionasal length (y2), and greatest length of nasal (y3)$ of skull size were recorded.
The logs of the responses were recorded.
}
\usage{data("anteaters")}
\format{
	A data frame consisting of 13 observations and 4 variables.
	\describe{
		\item{\code{y1}}{Basal length}
		\item{\code{y2}}{Occipitionasal length}
		\item{\code{y3}}{Greatest length of nasal}
		\item{\code{area}}{ (1) Minas Graes, Brazil; (2) Matto Grosso, Brazil;
				(3) Santa Cruz, Bolivia}
	}
}

\details{
Seber (1984) discusses a study on anteaters concerning the skull
size of anteaters from three different areas of South America.
The original source of the data is Reeve (1941).
Three measurements (basal length (y1), occipitionasal length (y2), and greatest length of nasal (y3)$ of skull size were recorded.
The logs of the responses were recorded.
The areas sampled were: Minas Graes, Brazil;
Matto Grosso, Brazil; and Santa Cruz, Bolivia.

CAUTION: Small sample sizes!
}
\source{
  Reeve, E.C.R. (1941), A statistical analysis of taxonomic differences within the genus \emph{Tamandua} Gray (Xenarthra), \emph{Proc.\ Zool.\ Soc.\ London A}, 111, 279--302.
}
\references{
Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.

Seber, G.A.F., (1984), \emph{Multivariate Observations}, NY: John Wiley.
}
\examples{
data(anteaters)
boxplot(y1~area,data=anteaters)
}
\keyword{datasets}
