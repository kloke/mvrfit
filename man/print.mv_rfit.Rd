\name{print.mv_rfit}
\alias{print.mv_rfit}
\title{
mv_rfit helper function
}
\description{
Print the output in a user-friendly manner using the internal R function \code{print}
}
\usage{
\method{print}{mv_rfit}(x, ...)
}
\arguments{
  \item{x}{
object of class mv_rfit
}
  \item{\dots}{
additional arguments passed along to print...
}
}
\author{
John Kloke
}

\seealso{
\code{\link{mv_rfit}}
}
\examples{
## The function initially defined as
function (x, ...) 
{
    cat("Call:\n")
    print(x$call)
    coef <- coef(x)
    cat("\nCoefficients:\n")
    print(coef, ...)
  }
}
