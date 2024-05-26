cwmv <- function(x,method='HL',H=diag(rep(1,ncol(x))),g0=rep(0,ncol(x)), C=diag(rep(1,ncol(x))),
                      test=FALSE, test.method="Wald", 
                      conf.int=FALSE, conf.level=0.95, conf.adjust='none', ...
                     ) {

# nonparametric multivariate analysis
# as discussed in Chapter in Kloke & McKean (2nd)

# componentwise Wald estimation and testing

# input:
#   x - n x k matrix of responses where the rows represent multivariate responses (i.e. experimental units) and the columns represent measurements (i.e. variables).
#   (optional) y - m x k matrix of responses for a second treatment group
#   method
#   H - k x q dimensional contrast matrix
#   g0 - null value of Hypothesis
#   C - 
#
# conf.adjust one of 'none', 'bonferroni'

# author: johndkloke@gmail.com
# maintainer: Rfit.developer@gmail.com

# details:
# nonparametric componentwise estimation and Wald inference for one-sample multivariate problems
# 
# One sample estimation and testing may be median, HL:HodgesLehmann (Wilcoxon), mean (LS:LeastSquares)
# tests the hypothesis of 
#   H0: H theta = g0
#   HA: H theta ne g0
# 

# recommends/depends:  

  n <- nrow(x)
  k <- nrow(H)

# check matrix rank
if(qr(H)$rank != k) stop("Hypothesis matrix H not full rank.")

  est <- npmvEst(x,method=method)
  Sn <- npmvEstVar(x,method=method,est=est) 

  result <- list(estimate=est,vcov=Sn,test=test,conf.int=conf.int,method=method)

  if( test ) {
    if( (method == 'LS') & test.method != 'Wald') {
      warning("setting test.method to 'Wald' for Least Squares (LS)")
      test.method <- 'Wald'
    }
    if( test.method == 'Wald' ) {

      a <- H%*%est-g0
      v <- H%*%Sn%*%t(H)

      T2.statistic <- t(a)%*%solve(v)%*%a
      F.statistic <- (n-k)/((n-1)*(k))*T2.statistic
      p.value <- pf(F.statistic,k,n-k,lower.tail=FALSE)

      result$F.statistic <- F.statistic
      result$p.value <- p.value
      result$g0 <- g0

#      result$test <- list(F.statistic=F.statistic, p.value=p.value, g0=g0)
    } else if( test.method == 'Gradient' ) {

      test_results <- gradient.test(x, g0=g0, sign_rank=(method=='HL'))
      result$F.statistic <- test_results$F.statistic
      result$p.value <- test_results$p.value
      result$g0 <- test_results$g0

    } else {
      stop( "test.method must be either 'Wald' or 'Gradient'" )
    }
    result$test.method <- test.method
  }
 

  if( conf.int ) {

    c.est <- C%*%est
    c.se <- sqrt(diag(C%*%Sn%*%t(C)))

    k.ci <- nrow(C)
    alphastar <- (1-conf.level)/(2*ifelse(conf.adjust=='bonferroni',k.ci,1))

    tcv <- qt(alphastar,nrow(x)-1,lower.tail=FALSE)

    CIs <- cbind(c.est-tcv*c.se,c.est+tcv*c.se)
    colnames(CIs) <- c('ci.lower','ci.upper')

    if(!is.null(colnames(x))) rownames(CIs) <- colnames(x)

    result$conf.int.tab <- CIs
    result$conf.adjust <- conf.adjust

  }

  class(result) <- 'cwmv'
  result

}

