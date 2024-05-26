"coef.mv_rfit" <-
function(object,...) { 
  x <- object
  coef <- matrix(unlist(x$coef),ncol=length(x$fits))
  rownames(coef) <- names(x$coef[[1]])
  colnames(coef) <- colnames(x$y)
  return(coef)
}
