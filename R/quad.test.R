quad.test <- function(fit,M,K,...) {

  if( (!is.matrix(M)) | (!is.matrix(K)) ) stop("invalid data type on input")

  q <- nrow(M)
  s <- ncol(K)

  b <- coef(fit)[-1,]
  mbk <- M%*%b%*%K
  v1 <- solve( M %*% fit$xpxi %*% t(M) )
  v2 <- solve( t(K) %*% fit$TST %*%K )

  test.statistic <- sum(diag( ( t(mbk) %*% v1 %*% mbk %*% v2 ) ))

  p.value <- pchisq(test.statistic,q*s,lower.tail=FALSE)

  result <- list(test.statistic=test.statistic,p.value=p.value,M=M,K=K)

  class(result) <- 'mv_rfit.quad.test'

  result

}
