
rm(list=ls())
setwd('E:\\uthealth\\deconv\\20171009')

##############################
# Deconvolution
##############################

load('data4deconv_5cell.rda')
source('..\\CIBERSORT.R')

r1 <- ref2deconv
e1 <- expr2deconv
deconv.rlt <- CIBERSORT2(r1,e1)
save(deconv.rlt,file='rlt4causal_5cell.rda')

##############################
# Causal
##############################

rm(list=ls())
library(assist)
library(dHSIC)
library(fANCOVA)
library(psych)
library(bnlearn)
library(infotheo)
library(descr)
load("E:/uthealth/getpathway/gene39761.rda")
rownames(disease) <- disease[,1]
disease[,5] <- ifelse(disease[,5]%in%c(4,5),1,0); disease <- disease[,-1]
disease <- disease[,-1]

setwd('E:\\uthealth\\deconv\\20171009')
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

ad <- disease[,3]
rlt.ttest <- apply(deconv.rlt,2,function(x){
  t.test(x~ad)$p.value
})
set.seed(12345);rlt.causal <- apply(deconv.rlt,2,function(x){
  unlist(permcontidANM(x,ad,1000))
}); rlt.causal


