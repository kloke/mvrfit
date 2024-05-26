twoway.mvrfit <- function(y,g,...) {

emat <- function(r) {
  d <- rep(1,r)
  return(cbind(d,diag(d)[,-1]))
}

g <- as.data.frame(g)
if(is.null(colnames(g))) {
  fn <- c('f1','f2') 
} else {
  fn <- colnames(g)
}

g1 <- factor(paste(paste(fn[1],g[,1],sep='_'), paste(fn[2],g[,2],sep='_'), sep='.'))

x <- model.matrix(~g1-1)

p <- ncol(x)
x1 <- x[,-1]
fit <- mv_rfit(x1,y,...)

muhat <- coef.mv_rfit(fit)
for(j in 2:p) {
  muhat[j,] <- muhat[1,] + muhat[j,]
}
rownames(muhat) <- levels(g1)

m <- ncol(y)
kmat <- diag(1,nrow=m,ncol=m)

listtests <- subsets(2)

atab <- matrix(NA,nrow=3,ncol=3)
colnames(atab) <- c("Test_Statistic","DF","p-value")
rownames(atab) <- c(fn,paste(fn,collapse='x'))

for(i in 1:3) {

  hmat <- khmat(g,listtests[i,])
  hmat2 <- hmat%*%emat(p)[,-1,drop=FALSE]
  ti <- quad.test(fit,hmat2,kmat)
  dfi <- nrow(hmat2)*ncol(kmat)
  atab[i,] <- c(ti[[1]],dfi,ti[[2]])

}

result <- list(anova.table=atab,fit=fit,muhat=muhat)
class(result) <- 'twoway.mvrfit'
result

}
