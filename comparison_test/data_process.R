
rm(list=ls())
setwd('/home/zhu/deconv/data')
library(data.table)

#macro

pca <- function(X,prop=1){
  X <- scale(X)[,]
  m = nrow(X)
  n = ncol(X)
  Xeigen <- svd(X)
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  score[,which(value>=prop)[1],drop=F]
}

#getdata

raw <- lapply(dir(),fread)

#adnidata
  raw1 <- as.matrix(raw[[1]])
  expr <- raw1[-1:-8,-1:-3]
  colnames(expr) <- raw1[2,-1:-3]
  rownames(expr) <- raw1[-1:-8,3]
  gene <- rownames(expr)
  ginfo <- expr[nchar(gene)>0,ncol(expr)]
  expr <- expr[nchar(gene)>0,-ncol(expr)]
  gene <- rownames(expr)
  gene <- strsplit(gene,' \\|\\| ')
  adnidata <- list(expr=expr,gene=gene)
