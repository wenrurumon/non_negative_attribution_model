
rm(list=ls())
setwd('E:\\uthealth\\deconv')
load("E:/uthealth/getpathway/gene39761.rda")

########################
# Function
########################

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
qpca2 <- function(A,prop1=0.9){
  rank1 <- qpca(A,T,0)$prop
  qpca(A,T,rank=which(rank1>=0.9)[1])
}
pca <- function(X){
  X <- as.matrix(X)
  m = nrow(X)
  n = ncol(X)
  X = scale(X)
  Xeigen <- svd(as.matrix(X))
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  mat <- Xeigen$v
  list(score=score,prop=value,mat=mat)
}

########################
# UCSD data
########################

raw.refuc <- do.call(rbind,strsplit(readLines('15.anno_gene_log2tpm_protCod.dat_nc'),'\t'))
gene1 <- c(T,raw.refuc[-1,1]%in%colnames(raw_exp))
raw.mapuc <- do.call(rbind,strsplit(readLines('reference_cell.csv'),','))
refuc <- scale(apply(raw.refuc[,-1],2,function(x){as.numeric(x[-1])}))
colnames(refuc) <- raw.refuc[1,-1]
gene1 <- table(raw.refuc[-1,1])

refuc1 <- refuc[match(names(which(gene1==1)),raw.refuc[-1,1]),]
rownames(refuc1) <- names(which(gene1==1))
refuc2 <- refuc[which(raw.refuc[-1,1] %in% names(which(gene1>1))),]
refuc2 <- apply(refuc2,2,function(x){
  tapply(x,raw.refuc[-1,1][which(raw.refuc[-1,1] %in% names(which(gene1>1)))],mean)
})
refuc <- rbind(refuc1,refuc2)
colnames(refuc) <- substr(paste(raw.mapuc[match(colnames(refuc),raw.mapuc[,1]),5]),1,1)
ref1 <- lapply(unique(colnames(refuc)),function(x){
  refuc[,colnames(refuc)==x]
})
gene1 <- names(gene1)

########################
# CM data
########################

setwd('E:\\uthealth\\deconv\\20171009')
raw <- lapply(c("cortex.csv","microglia.csv"),read.csv)
gene <- paste(raw[[1]][,1])%in%colnames(raw_exp)
raw <- lapply(raw,function(x){x[gene,]})
gene <- paste(raw[[1]][,1])
ref2 <- lapply(raw,function(x){
  scale(apply(x[,-1],2,function(y){
    tapply(y,gene,mean)
  }))
})
gene2 <- gene

########################
# All Ref
########################

ref <- c(ref1,ref2)
names(ref) <- c(unique(colnames(refuc)),'C','M')
gene <- gene1[gene1%in%gene2]
ref <- lapply(ref,function(x){
  x[match(gene,rownames(x)),]
})
expr <- t(raw_exp[,match(gene,colnames(raw_exp))])

save(ref,expr,file='data45cell.rda')



