
rm(list=ls())
library(assist)
library(dHSIC)
library(fANCOVA)
library(psych)
library(bnlearn)
library(infotheo)
library(descr)

setwd('/zfsauton/project/highmark/aska/deconv')
source("contiANM_1.3.R")
source("contidANM.R")
source("discreteANM_1.1.R")
load("mdata4causal.rda" )
load('deconv2m.rda')
rlt <- do.call(rbind,rlt)

ttest <- function(x,y){
  if(var(x)==0){
    return(NA)
  }else{
    return(t.test(scale(x)~y)$p.value)
  }
}

causal.ad <- function(x){
  r1 <- c(unlist(permcontidANM(scale(x),ad,10000))[1:2],
          P_ttest=ttest(x,ad))
  x <- x * mprop
  r2 <- c(unlist(permcontidANM(scale(x),ad,10000))[1:2],
          P_ttest=ttest(x,ad))
  c(r1,r2)
}

sel1 <- rlt[,2]<=0.05 | rlt[,5]<=0.05
sel2 <- rlt[,3]<=(0.05/40120) | rlt[,6]<=(0.05/40120)
sel2[is.na(sel2)] <- FALSE
sel <- sel1|sel2
sel <- which(sel)
rm(sel1,sel2)

selcut <- as.numeric(cut(1:length(sel),40))

arg <- as.numeric(commandArgs(trailingOnly=TRUE))
if(length(arg)==0){
  arg<-10
  is <- sel[which(selcut==arg)][1:2]
} else {
  is <- sel[which(selcut==arg)]
}

rlt <- lapply(is,function(i){
	print(paste(arg,match(i,is)))
	x <- raw_exp[,colnames(raw_exp)==names(is)[match(i,is)]]
	try(causal.ad(x))
})
names(rlt) <- names(is)

save(rlt,file=paste0('deconv2m_',arg,'_r2.rda'))
