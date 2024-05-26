mv_rfit <- function(x,y,...) {

call <- match.call()

x.o <- x
y.o <- y

x <- as.data.frame(x)
y <- as.data.frame(y)

x.names <- names(x)
y.names <- names(y)

if(nrow(x) != nrow(y)) stop("incorrect data dimensions")

n <- nrow(y)

rfit1 <- function(y,...) Rfit::rfit(y~.,data=cbind(y,x),...)

fits <- lapply(as.list(as.data.frame(y)),rfit1,...)

# fits

l2m <- function(l,n,names) {
  m <- matrix(unlist(l),nrow=n)
  colnames(m) <- names
  return(m)
}

rtau <- function(x) return(x$tauhat)
rtaus <- function(x) return(x$taushat)
rdisp <- function(x) return(x$disp)

tauhatVec <- drop(l2m(lapply(fits,rtau),1,y.names))

Emat <- l2m(lapply(fits,residuals),n,y.names)
R <- apply(Emat,2,rank,ties.method='first')
gs1 <- function(x,scores) getScores(scores,x)
scores <- fits[[1]]$scores
S1 <- apply(R/(n+1),2,gs1,scores=scores)
S <- crossprod(S1,S1)/(n+1)

xc <- scale(x,scale=FALSE,center=TRUE)
xpxi <- chol2inv(chol(crossprod(xc,xc)))

result <- list(coefficients=lapply(fits,coefficients),
               residuals=Emat,
               fitted.values=l2m(lapply(fits,fitted.values),n,y.names),
               scores=scores,
               x=x.o,y=y.o,
               tauhatVec=tauhatVec,
               tauhatVec=drop(l2m(lapply(fits,rtau),1,y.names)),
               taushatVec=drop(l2m(lapply(fits,rtaus),1,y.names)),
               dispVec=drop(l2m(lapply(fits,rdisp),1,y.names)),
               S=S,xpxi=xpxi,TST=diag(tauhatVec)%*%S%*%diag(tauhatVec),
               fits=fits
              )


result$call <- call
class(result) <- list("mv_rfit")

result

}
