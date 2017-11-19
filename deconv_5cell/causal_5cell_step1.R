##############################
# Data Processing
##############################

rm(list=ls())

args <- as.numeric(commandArgs(trailingOnly=TRUE))
if(length(args)==0){args <- 1}

library(assist)
library(dHSIC)
library(fANCOVA)
library(psych)
library(bnlearn)
library(infotheo)
library(descr)
load("gene39761.rda")
rownames(disease) <- disease[,1]
disease[,5] <- ifelse(disease[,5]%in%c(4,5),1,0); 
disease <- disease[,-1]

load('data4deconv_5cell.rda')
load('data45cell.rda')
load('rlt4causal_5cell.rda')
source('discreteANM_1.1.R')
source('contiANM_1.3.R')
source('contidANM.R')

deconv.rlt <- deconv.rlt[,0:2-ncol(deconv.rlt)]
deconv.rlt <- t(apply(deconv.rlt,1,function(x){
  tapply(x,colnames(deconv.rlt),sum)
}))
deconv.rlt <- cbind(R=1,deconv.rlt)
ad <- disease[,3]

###################################
# Causal to 5
###################################

x <- apply(raw_exp,2,function(x){
  list(x * deconv.rlt)
})
x <- do.call(c,x)
sel <- as.numeric(cut(1:length(x),25))
x <- x[sel==args]

ctest <- function(x,y=ad,nop=100){
  set.seed(args)
  rlt <- apply(x,2,function(x){
    rlt1 <- unlist(permcontidANM(x,y,nop))[1:2]
    rlt2 <- try(t.test(x~y)$p.value)
    if(!is.numeric(rlt2)){rlt2 <- NA}
    rlt <- c(rlt1,p_ttest=rlt2)
  })  
  rlt
}

rlt <- lapply(x,ctest)
save(rlt,file=paste0('rlt_step1_',args,'.rda'))


