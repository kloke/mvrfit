print.mv_rfit <-
function(x,...) {
  cat("Call:\n")
  print(x$call)
  coef <- coef(x)
  cat("\nCoefficients:\n")
  print(coef, ...)

}
