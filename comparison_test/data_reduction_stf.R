
rm(list=ls())
setwd('/home/zhu/deconv/data')
source('/home/zhu/deconv/CIBERSORT.R')
load('quickload.rda')
 
pca.adni <- pca(A.adni)
pca.rush <- pca(A.rush)
mat.adni <- pca.adni$mat[,1:which(pca.adni$value>0.8)[1],drop=F]
mat.rush <- pca.rush$mat[,1:which(pca.rush$value>0.8)[1],drop=F]

#######################

A <- A.rush %*% mat.rush
X <- X.ref %*% mat.rush
gcode <- paste0('g',1:ncol(A))
colnames(A) <- colnames(X) <- gcode
A <- data.matrix(t(A))
X <- data.matrix(t(X))
rlt1.rush <- CIBERSORT2(X,A)
save(rlt1.rush,file='rlt1_rush.rda')

###########################

A <- A.adni %*% mat.adni
X <- X.ref %*% mat.adni
gcode <- paste0('g',1:ncol(A))
colnames(A) <- colnames(X) <- gcode
A <- data.matrix(t(A))
X <- data.matrix(t(X))
rlt1.adni <- CIBERSORT2(X,A)
save(rlt1.adni,file='rlt1_adni.rda')
