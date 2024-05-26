lincomb.mv_rfit <- function(fit,h,l,conf.level=0.95) {

  betahatMat <- coef.mv_rfit(fit)[-1,]
  est <- t(h)%*%betahatMat%*%l
  se <- sqrt(t(h)%*%fit$xpxi%*%h * t(l)%*%fit$TST%*%l) 
  pp1 <- ncol(fit$x) + 1
  n <- nrow(fit$y)
  tcv <- qt( (1-conf.level)/2, n-pp1, lower.tail=FALSE )
   
  lcl <- est - tcv*se
  ucl <- est + tcv*se
  tab <- c(est,se,lcl,ucl)
  names(tab) <- c('Est','SE','LCL','UCL')

  return(tab)

}
