
rm(list=ls())
source('/home/zhu/deconv/qnmf.R')
setwd('/home/zhu/deconv/data')
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
A.rush <- data.matrix(t(A))
X.rush <- data.matrix(t(X))

A <- A.adni %*% mat.adni
X <- X.ref %*% mat.adni
gcode <- paste0('g',1:ncol(A))
colnames(A) <- colnames(X) <- gcode
A.adni <- data.matrix(t(A))
X.adni <- data.matrix(t(X))

#######################

rlt2.rush <- qnmf(A=A.rush,K=NULL,X=X.rush,a=0.5,lambda=0,maxitn=10000,deconv=T)$Y
rlt2.adni <- qnmf(A=A.adni,K=NULL,X=X.adni,a=0.5,lambda=0,maxitn=10000,deconv=T)$Y
load('rlt1_adni.rda'); load('rlt1_rush.rda')

##################################

Ystf.adni <- rlt1.adni[,1:58]
Ystf.rush <- rlt1.rush[,1:58]
Ynmf.adni <- t(rlt2.adni)
Ynmf.rush <- t(rlt2.rush)

cell <- colnames(Ystf.adni)
cell <- substr(cell,1,regexpr(' ',cell)-1)

agg2cell <- function(x){
  t(apply(x,1,function(xi){tapply(xi,cell,sum)}))
}

Ystf.adni <- agg2cell(Ystf.adni)
Ystf.rush <- agg2cell(Ystf.rush)
Ynmf.adni <- agg2cell(Ynmf.adni)
Ynmf.rush <- agg2cell(Ynmf.rush)

summary(sapply(1:nrow(Ystf.adni),function(i){cor.test(Ystf.adni[i,],Ynmf.adni[i,])$p.value}))
summary(sapply(1:nrow(Ystf.rush),function(i){cor.test(Ystf.rush[i,],Ynmf.rush[i,])$p.value}))
summary(sapply(1:ncol(Ystf.adni),function(i){cor.test(Ystf.adni[,i],Ynmf.adni[,i])$p.value}))
summary(sapply(1:ncol(Ystf.rush),function(i){cor.test(Ystf.rush[,i],Ynmf.rush[,i])$p.value}))

save(Ystf.adni,Ynmf.adni,Ystf.rush,Ynmf.rush,file='rlt4compare.rda')
