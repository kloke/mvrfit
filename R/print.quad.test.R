print.mv_rfit.quad.test <- function(x,digits=max(5, .Options$digits - 2),...) {
  cat("\nMultivariate Rank-Based Regression Quadratic Test\n")

  y <- c(x$test.statistic, x$p.value)
  names(y) <- c("Test_Statistic", "p-value")

  print(format(y, digits = digits), quote = FALSE)


}
