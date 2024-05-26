print.twoway.mvrfit <- function(x,digits=max(5, .Options$digits - 2),...) {
  cat("\nMultivariate Rank-Based Regression Quadradic Tests\n")
  cat("Main Effects and Interaction\n")
  print(format(x$anova.table, digits = digits), quote = FALSE)
}
