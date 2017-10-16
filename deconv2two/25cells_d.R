
rm(list=ls())
library(assist)
library(dHSIC)
library(fANCOVA)
library(psych)
library(bnlearn)
library(infotheo)
library(descr)

setwd('/home/zhu/deconv/25cells')
source("contiANM_1.3.R")
source("contidANM.R")
source("discreteANM_1.1.R")
load("mdata4causal.rda" )

ttest <- function(x,y){
  if(var(x)==0){
    return(NA)
  }else{
    return(t.test(scale(x)~y)$p.value)
  }
}

causal.ad <- function(x){
  r1 <- c(unlist(permcontidANM(scale(x),ad,100))[1:2],
          P_ttest=ttest(x,ad))
  x <- x * mprop
  r2 <- c(unlist(permcontidANM(scale(x),ad,100))[1:2],
          P_ttest=ttest(x,ad))
  c(r1,r2)
}

idx <- 1:ncol(raw_exp)
idx <- cut(idx,20)
arg <- as.numeric(commandArgs(trailingOnly=TRUE))
is <- which(idx==unique(idx)[arg])

rlt <- lapply(is,function(i){
  print(i)
  x <- raw_exp[,i]
  try(causal.ad(x))
})
names(rlt) <- colnames(raw_exp)[is]

save(rlt,file=paste0('deconv2m_',arg,'.rda'))
