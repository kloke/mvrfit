\name{amidata}
\alias{amidata}
\docType{data}
\title{
Study of the Drug Amitriptyline}
\description{
Study on 17 patients who were admitted to the hospital after taking an overdose of the antidepressant amitriptyline.
The two response variables  recorded are  total TCAD plasma level and the amount of amitriptyline in TCAD.  Several predictors were also recorded: gender, PR (electrocardiogram measuremen), amount of antidepressants taken at the time of overdose, diastolic blood pressure, QRS wave measurements.
}
\usage{data("amidata")}
\format{
A data frame with 17 observations and 7 variables.
\describe{
	\item{\code{tot}}{Total TCAD, response}
	\item{\code{ami}}{Amount of amitriptyline in TCAD, response}
	\item{\code{gen}}{Gender, predictor}
	\item{\code{amt}}{Amount of antidepressants taken at the time of overdose, predictor}
	\item{\code{pr}}{Electrocardiogram measurement, predictor}
	\item{\code{diap}}{Diastolic blood pressure, predictor}
	\item{\code{qrs}}{QRS wave measurement, predictor}
}
}
\details{
        Johnson and Wichern (2007) present the results of a study on 17 patients who were admitted to the hospital after taking an overdose of the antidepressant amitriptyline.
        The two response variables  recorded are  total TCAD plasma level
        and the amount of amitriptyline in TCAD.
        Several predictors were also recorded: gender, PR (electrocardiogram measuremen), amount of antidepressants taken at the time of overdose, diastolic blood pressure, QRS wave measurements.

}
\source{
	Johnson, R.A. and Wichern, D.W., (2007), \emph{Applied multivariate statistical analysis, 6th ed.}, Upper Saddle River, NJ: Pearson.
}
\references{
Kloke, J. and McKean J.W. (2024), \emph{Nonparametric Statistical Methods using R, 2nd ed. }, Boca Raton, FL: Chapman-Hall.
}
\examples{
data(amidata)
opar_mfrow <- par()$mfrow
par(mfrow=c(1,2))
plot(tot~amt,data=amidata)
plot(ami~amt,data=amidata)
par(mfrow=opar_mfrow)
}
\keyword{datasets}
