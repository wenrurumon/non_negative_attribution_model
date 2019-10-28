
rm(list=ls())
library(data.table)
library(dplyr)
source('/Users/wenrurumon/Documents/uthealth/deconv/ae/model.R')

#PosLM

PosLM <- function(Yi,X){
  Yi <- minmax(Yi)
  X <- apply(X,2,minmax)
  colnames(X) <- 1:ncol(X)
  Xi <- X
  while(T){
    Mi <- coef(summary(lm(Yi~Xi)))[-1,]
    sel <- (which.min(Mi[,3]))
    Xi <- Xi[,-sel,drop=F]
    if(Mi[sel,1] > 0){break}
  }
  outi <- rep(0,ncol(X))
  outi[as.numeric(gsub('Xi','',rownames(Mi)))] <- Mi[,1]
  outi <- colSums(X) * outi
  outi/sum(outi)
}

#Dummy data

x <- apply(matrix(rnorm(1000),100,10),2,minmax)
p <- apply(matrix(rnorm(200),20,10),2,minmax)
dimnames(x) <- list(paste0('sample',1:100),paste0('cell',1:10))
p[p<0.5] <- 0
p <- apply(p,1,function(x){x/sum(x)})
y <- x %*% p
y <- apply(y,2,minmax)
dimnames(y) <- list(paste0('sample',1:100),paste0('gene',1:20))
sum(x);sum(y)

#STF

out_stf <- deconv(x,y)$coef
diag(cor(out_stf,p))

#LM

out_lm <- apply(y,2,function(y){PosLM(y,X=x)})

#Validate

mse(t(dedeconv(t(x),out_stf,t(y))),y)
mse(t(dedeconv(t(x),out_lm,t(y))),y)

