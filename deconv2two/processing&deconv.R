
rm(list=ls())
load("E:/uthealth/getpathway/gene39761.rda")
setwd('E:\\uthealth\\deconv\\20171009')

###########################################
# Process data
###########################################

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

raw <- lapply(dir(pattern='.csv'),read.csv)
gene <- paste(raw[[1]][,1])%in%colnames(raw_exp)
raw <- lapply(raw,function(x){x[gene,]})
gene <- paste(raw[[1]][,1])

raw <- lapply(raw,function(x){
  apply(x[,-1],2,function(y){
    tapply(y,gene,mean)
  })
})
gene <- unique(gene)
raw.pca <- lapply(raw,function(x){qpca2(x,0.9)$X})
colnames(raw.pca[[1]]) <- 
  paste0(gsub('.csv','',dir(patter='.csv')[1]),1:ncol(raw.pca[[1]]))
colnames(raw.pca[[2]]) <- 
  paste0(gsub('.csv','',dir(patter='.csv')[2]),1:ncol(raw.pca[[2]]))

ref <- do.call(cbind,raw.pca); rownames(ref) <- gene
expr <- t(raw_exp[,match(gene,colnames(raw_exp))])

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

expr.pca <- pca(t(expr))
save(ref,expr,expr.pca,file='data4deconv.rda')

###########################################
# Setup Deconvolution
###########################################

rm(list=ls())
setwd('E:\\uthealth\\deconv\\20171009')
load('data4deconv.rda')
source('..\\CIBERSORT.R')

# mat1 <- expr.pca$mat[,1:which(expr.pca$prop>=0.9)[1],drop=F]
mat2 <- expr.pca$mat[,1:which(expr.pca$prop>=0.95)[1],drop=F]

# r1 <- t(scale(t(ref)) %*% mat1)
# e1 <- t(scale(t(expr)) %*% mat1)
# rownames(r1) <- rownames(e1) <- paste0('g',1:nrow(r1))
# rlt1 <- CIBERSORT2(r1,e1)

r2 <- t(scale(t(ref)) %*% mat2)
e2 <- t(scale(t(expr)) %*% mat2)
rownames(r2) <- rownames(e2) <- paste0('g',1:nrow(r2))
rlt2 <- CIBERSORT2(r2,e2)

dcv.rlt <- rlt2
save(dcv.rlt,file='rlt4causal.rda')
