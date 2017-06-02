
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
raw <- lapply(dir()[!grepl('rda',dir())],fread)

#adnidata
  raw1 <- as.matrix(raw[[2]])
  expr <- raw1[-1:-8,-1:-3]
  colnames(expr) <- raw1[2,-1:-3]
  rownames(expr) <- raw1[-1:-8,3]
  gene <- rownames(expr)
  ginfo <- expr[nchar(gene)>0,ncol(expr)]
  expr <- expr[nchar(gene)>0,-ncol(expr)]
  gene <- rownames(expr)
  gene <- strsplit(gene,' \\|\\| ')
  adnidata <- list(A=expr,gene=gene,sample=colnames(expr))
  rm(expr,gene,ginfo,raw1)

#rushdata
  raw1 <- as.matrix(raw[[5]])
  expr <- raw1[,c(-ncol(raw1),-1)]
  rownames(expr) <- raw1[,ncol(raw1)]
  colnames(expr) <- raw[[6]]$genenm
  gene <- colnames(expr)
  rushdata <- list(A=expr,gene=colnames(expr),sample=rownames(expr))
  rm(expr,gene,raw1)

#reference
  raw1 <- as.matrix(raw[[1]])
  X <- raw1[,-1]
  rownames(X) <- raw1[,1]
  Xinfo <- raw[[4]]
  X <- X[,colnames(X)%in%Xinfo[[1]],drop=F]
  Xcell <- Xinfo[[4]][match(colnames(X),Xinfo[[1]])]
  refdata <- list(X=X,cell=Xcell,info=Xinfo,gene=rownames(X))
  rm(raw1,X,Xcell,Xinfo)

#overlap gene
  gene1 <- unique(unlist(adnidata$gene))
  gene2 <- unique(rushdata$gene)
  gene3 <- unique(refdata$gene)
  gene_sel <- names(which(table(c(gene1,gene2,gene3))==3))
  rm(gene1,gene2,gene3)

#save data
save(rushdata,adnidata,refdata,gene_sel,file='loaddata.rda')
  
