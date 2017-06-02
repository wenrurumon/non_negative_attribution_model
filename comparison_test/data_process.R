
rm(list=ls())
setwd("/home/zhu/deconv/data")
load('loaddata.rda')

pca1 <- function(X,prop=1){
  X <- scale(X)[,]
  m = nrow(X)
  n = ncol(X)
  Xeigen <- svd(X)
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  score[,1]
}

#adnidata
  A <- adnidata$A
  g <- adnidata$gene
  Asplit <- lapply(1:length(g),function(i){
    print(i/length(g))
    temp <- matrix(rep(as.numeric(A[i,,drop=F]),length(g[[i]])),ncol=length(g[[i]]),dimnames=list(colnames(A),g[[i]]))
    temp
  })
  A <- do.call(cbind,Asplit)
  g <- unique(colnames(A))
  g <- g[g%in%gene_sel]
  A <- A[,colnames(A)%in%g]
  Asplit <- sapply(g,function(gi){
    temp <- A[,colnames(A)==gi,drop=F]
    pca1(temp)
  })
adnidata$A2 <- Asplit

#rushdata
  A <- rushdata$A
  g <- unique(rushdata$gene)
  g <- g[g%in%gene_sel]
  A <- apply(A,2,as.numeric)
