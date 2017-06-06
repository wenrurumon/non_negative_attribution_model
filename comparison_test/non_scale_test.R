
rm(list=ls())
setwd("/home/zhu/deconv/data")

################
# Process data
################

load('loaddata.rda')

pca1 <- function(X,prop=1){
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
  A <- apply(A,2,as.numeric)
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
  X <- do.call(cbind,lapply(unique(refdata$cell),function(ci){
    print(ci)
    out <- pca(X[,refdata$cell==ci,drop=F],prop=0.75)$score2
    colnames(out) <- paste(ci,1:ncol(out))
    out
  }))
  Xsplit <- lapply(gene_sel,function(gi){
    out <- X[rownames(X)==gi,,drop=F]
    if(nrow(out)>1){
      out <- t(pca1(t(out)))
    }
    out
  })
  X2 <- t(do.call(rbind,Xsplit))
  colnames(X2) <- gene_sel
  refdata$X2 <- X2

################
#quickload
################

A.adni <- adnidata$A2
A.rush <- rushdata$A2
X.ref <- refdata$X2
#save(A.adni,A.rush,X.ref,pca,file='quickload2.rda')

################
#pca mat
################

A.adni <- scale(A.adni); A.adni[is.na(A.adni)] <- 0
A.rush <- scale(A.rush); A.rush[is.na(A.rush)] <- 0
pca.adni <- pca(A.adni)
pca.rush <- pca(A.rush)
mat.adni <- pca.adni$mat[,1:which(pca.adni$value>0.8)[1],drop=F]
mat.rush <- pca.rush$mat[,1:which(pca.rush$value>0.8)[1],drop=F]

################
#STF Method
################

setwd('/home/zhu/deconv/data')
source('/home/zhu/deconv/CIBERSORT.R')

A <- A.rush %*% mat.rush
X <- X.ref %*% mat.rush
gcode <- paste0('g',1:ncol(A))
colnames(A) <- colnames(X) <- gcode
A <- data.matrix(t(A))
X <- data.matrix(t(X))
Yrs <- CIBERSORT2(X,A)

A <- A.adni %*% mat.adni
X <- X.ref %*% mat.adni
gcode <- paste0('g',1:ncol(A))
colnames(A) <- colnames(X) <- gcode
A <- data.matrix(t(A))
X <- data.matrix(t(X))
Yas <- CIBERSORT2(X,A)

################
#QNMF Method
################

source('/home/zhu/deconv/qnmf.R')
A <- A.rush %*% mat.rush
X <- X.ref %*% mat.rush
gcode <- paste0('g',1:ncol(A))
colnames(A) <- colnames(X) <- gcode
A.rush <- data.matrix(t(A))
X.rush <- data.matrix(t(X))

A <- A.adni %*% mat.adni
X <- X.ref %*% mat.adni
gcode <- paste0('g',1:ncol(A))
colnames(A) <- colnames(X) <- gcode
A.adni <- data.matrix(t(A))
X.adni <- data.matrix(t(X))

Yrn <- qnmf(A=A.rush,K=NULL,X=X.rush,a=0.5,lambda=0,maxitn=10000,deconv=T)$Y
Yan <- qnmf(A=A.adni,K=NULL,X=X.adni,a=0.5,lambda=0,maxitn=10000,deconv=T)$Y

###############
# Summary
###############

Yas2 <- Yas[,1:153]
Yan2 <- t(Yan)
Yrs2 <- Yrs[,1:153]
Yrn2 <- t(Yrn)
save(Yas2,Yan2,Yrs2,Yrn2,file='rlt4compare3.rda')


