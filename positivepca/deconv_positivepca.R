
rm(list=ls())
library(corpcor)
diag2 <- function(x){
  if(length(x)==1){
    return(as.matrix(x))
  } else {
    return(diag(x))
  }
}
ginv<-function(A){
  A_svd<-fast.svd(A)
  if(length(A_svd$d)==1){
    A_inv<-A_svd$v%*%as.matrix(1/A_svd$d)%*%t(A_svd$u)
  }else{
    A_inv<-A_svd$v%*%diag(1/A_svd$d)%*%t(A_svd$u)
  }
  return(A_inv)
}
positive <- function(x){
  ifelse(x>0,x,0)
}
qpca <- function(A,scale=T,rank=0){
  if(scale){A <- scale(A)}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-10]
  r <- length(d)
  prop <- d^2; prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop)
  return(rlt)
}


#input
set.seed(12); raw <- A <- rmatrix(1000,100)
# A <- scale(cbind(USArrests,iris[1:50*3,1:4]))[,]
K <- 5
lambda = 100

#initialization
  # A.svd <- svd(qpca(A,F,K)$Z)
  A.svd <- svd(scale(A))
    # cor(A,A.svd$u %*% diag(A.svd$d) %*% t(A.svd$v))
  U <- A.svd$u[,1:K,drop=F]
  V <- A.svd$v[,1:K,drop=F]
  D <- diag2(A.svd$d[1:K])
    # cor(A,U %*% diag2(A.svd$d[1:K]) %*% t(V))
  
  X <- U %*% sqrt(D)
  Y <- sqrt(D) %*% t(V)
  Iy <- (Y<0)
    # cor(A,X%*%Y)

#Loops
loop <- list()
for(i in 1:100000){
  Y2 <- positive(ginv(t(X)%*% X) %*% t(X) %*% A - lambda * Iy)
  X2 <- positive(A %*% t(Y2) %*% ginv(Y2 %*% t(Y2)))
  loop[[i]] <- list(X2=X2,Y2=Y2,
                    Xf=matrixcalc::frobenius.norm(X2-X),
                    Yf=matrixcalc::frobenius.norm(Y2-Y))
  X <- X2; Y <- Y2  
  lambda <- lambda * 0.9
}
A2 <- X2 %*% Y2
plot.ts(xf <- sapply(loop,function(x){x$Xf}))
plot.ts(yf <- sapply(loop,function(x){x$Yf}))

#nmf

test.nmf <- nmf(A[,],3)
A.nmf <- basis(test.nmf) %*% coef(test.nmf)
diag(cor(A.nmf,A))
diag(cor(A2,A))

sum((A-A.nmf)^2)
sum((A-A2)^2)
test <-sapply(loop,function(x){sum((x$X2 %*% x$Y2 - A)^2)})
plot.ts(cbind(xf,yf,test))
summary(cbind(xf,yf,test))
