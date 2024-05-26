summary.mv_rfit <-
function(object,...) { 
#  if( verbose ) sf <- summary
#  else 
  sf <- function(x,...) return(summary(x,...)$coef)
  lapply(object$fits,sf,...)
}
