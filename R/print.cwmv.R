print.cwmv <- function(x, digits = max(5, .Options$digits - 2), ...) {

  method <- switch(x$method,
    HL='Hodges-Lehmann',
    Median='Median',
    LS='Least Squares')

  cat("\nComponentwise Multivariate Procedure\n")
  cat("\nEstimation method:", method)

  cat("\n\nEstimated Location Vector:\n")
  print(format(x$estimate,digits=digits),quote=FALSE)

  if(x$test) {

    cat("\n", x$test.method, " Test\n")
    y <- c(x$F.statistic, x$p.value)
    names(y) <- c("F-Statistic", "p-value")
    print(format(y, digits = digits), quote = FALSE)
    cat("Null Location Vector:\n")
    print(format(x$g0, digits = digits), quote = FALSE)

  }


  if(x$conf.int) {
    cat("\nWald Confidence Intervals\n")
    print(format(x$conf.int.tab, digits = digits), quote = FALSE)
    cat("CI adjust method:", x$conf.adjust)
    cat("\n")
  }


  
  


}
