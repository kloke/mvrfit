# internal functions for cwmv

hle_loc <- function(x) {

# original code from wilcox.test in
# R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
# under GPL 2 or 3

  psums <- outer(x, x, "+")
  median( psums[!lower.tri(psums)] / 2 )

}

hle <- signedrank

npmvEst <- function(x,method='HL') {

  func <- switch(method,
    HL = hle,
    Median = median,
    LS = mean
  )

  apply(x,2,func)

}

npmvEstVar <- function(x,method='HL',est=npmvEst(x,method)) {

var_hle_0 <- function(x,est=npmvEst(x,'HL')) {

  x <- as.matrix(x)
  n <- nrow(x)

  r <- sweep(x,2,est)
  sr <- sign(r) * apply(abs(r),2,rank)
  Esr <- crossprod(sr)/((n+1)*(n+1)*n)*3
  diag(Esr) <- 1

  taus <- apply(x,2,gettau,0)
  tauMat <- tcrossprod(taus)

  S <- Esr * tauMat
  S/n

}

var_hle_1 <- function(x,est=npmvEst(x,'HL')) {

  x <- as.matrix(x)
  n <- nrow(x)

  r <- sweep(x,2,est)
  sr <- apply(r,2,rank)

  Esr <- 3*(crossprod(sr)/((n+1)*(n+1)*n)*4-1)
  diag(Esr) <- 1

  taus <- apply(x,2,gettau,0)
  tauMat <- tcrossprod(taus)

  S <- Esr * tauMat
  S/n

}

var_hle <- var_hle_1

var_median <- function(x,est=npmvEst(x,'Median')) {

  x <- as.matrix(x)
  n <- nrow(x)

  sr <- sign(sweep(x,2,est))
  Esr <- crossprod(sr)/n
  diag(Esr) <- 1

  taus <- apply(x,2,taustar,0)
  tauMat <- tcrossprod(taus)

  S <- Esr * tauMat
  S/n

}

var_ls <- function(x,...) {
  var(x)/nrow(x)
}

  func <- switch(method,
    HL = var_hle,
    Median = var_median,
    LS = var_ls
  )
 
  func(x,est)

}

gradient.test <- function(x, g0=rep(0,ncol(x)), sign_rank=TRUE, ...) { 

  x <- as.matrix(x)

  n <- nrow(x)
  k <- ncol(x)

  xc <- sweep(x,2,g0)
  sign_mat <- sign(xc)

  if( sign_rank ) {
    r_mat <- apply(abs(xc),2,rank)
    Svec <- apply( sign_mat*r_mat, 2, sum )/(n+1)
    rank_mat <- apply(xc,2,rank)
    Bmat <- crossprod(rank_mat)/(n*(n+1)*(n+1))*4-1
    diag(Bmat) <- 1/3
  } else {
    Svec <- apply( sign_mat, 2, sum )
    Bmat <- crossprod(sign_mat)/n
    diag(Bmat) <- 1
  }

  T <- t(Svec)%*%solve(Bmat)%*%Svec/n
  F.statistic <- (n-k)/(k*(n-1))*T
  p.value <- pf(F.statistic,k,n-k,lower.tail=FALSE)

  result <- list(F.statistic=F.statistic, p.value=p.value, g0=g0)

  # to do add class for gradient test
  result

}

print.gradient.test<-function(x,...) {
	with(x,cat("F-statistic = ",F.statistic,", p-value = ", p.value, "\n"))
}


