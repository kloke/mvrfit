"rstudent.mv_rfit" <- function(model,...) { 
  mc <- match.call()
  fit <- mc$model
  if(is.null(mc$as_list)) as_list <- FALSE
  else as_list <- eval(mc$as_list) == TRUE
  rstudent_l <- lapply(fit$fits,rstudent)
  if(as_list) return(rstudent_l)
  rstudent <- matrix(unlist(rstudent_l),ncol=length(rstudent_l))
  colnames(rstudent) <- names(rstudent_l)
  return(rstudent)
}
