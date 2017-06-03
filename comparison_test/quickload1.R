
rm(list=ls())
setwd('/home/zhu/deconv/data')
load('processeddata.rda')

pca <- function(X,prop=1){
  X <- scale(X)[,]
  m = nrow(X)
  n = ncol(X)
  Xeigen <- svd(X)
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  score2 <- score[,1:which(value>=prop)[1],drop=F]
  list(score=score,score2=score2,value=value,mat=Xeigen$v)
}

A.adni <- adnidata$A2
A.rush <- rushdata$A2
X.ref <- refdata$X2
save(A.adni,A.rush,X.ref,pca,file='quickload.rda')
