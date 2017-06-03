
rm(list=ls())
setwd("/home/zhu/deconv/data")
load('loaddata.rda')

pca1 <- function(X,prop=1){
  X[is.na(X)] <- 0
  X <- scale(X)[,]
  X[is.na(X)] <- 0
  m = nrow(X)
  n = ncol(X)
  Xeigen <- svd(X)
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  score[,1]
}

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
  Asplit <- sapply(gene_sel,function(gi){
    temp <- A[,colnames(A)==gi,drop=F]
    pca1(temp)
  })
rownames(Asplit) <- rownames(A)
adnidata$A2 <- Asplit

#rushdata
  A <- rushdata$A
  A <- A[,colnames(A)%in%gene_sel]
  A <- scale(apply(A,2,as.numeric))
  rownames(A) <- rownames(rushdata$A)
  i <- 0
  Asplit <- sapply(gene_sel,function(gi){
    print((i <<- i+1)/length(gene_sel))
    temp <- pca1(A[,colnames(A)==gi,drop=F])
  })
  rownames(Asplit) <- rownames(A)
  rushdata$A2 <- Asplit

#reference
  X <- apply(refdata$X,2,as.numeric)
  dimnames(X) <- dimnames(refdata$X)
  Xsplit <- lapply(unique(refdata$cell),function(ci){
    print(ci)
    pca(X[,refdata$cell==ci,drop=F],prop=0.6)
  })
