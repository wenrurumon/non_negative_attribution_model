##############################
# Data Processing
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
disease[,5] <- ifelse(disease[,5]%in%c(4,5),1,0); 
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
deconv.rlt <- cbind(R=1,deconv.rlt)

ad <- disease[,3]
# rlt.ttest <- apply(deconv.rlt,2,function(x){
#   t.test(x~ad)$p.value
# })
genesel <- colnames(raw_exp)[1:10]
expr_sel <- raw_exp[colnames(raw_exp)%in%genesel]

###################################
# Causal to 5
###################################

x <- apply(expr_sel,2,function(x){
  list(x * deconv.rlt)
})
x <- do.call(c,x)

ctest <- function(x,y=ad,nop=100){
  rlt <- apply(x,2,function(x){
    rlt <- unlist(permcontidANM(x,y,nop))
    rlt[!grepl('fct_',names(rlt))]
  })  
  rlt
}
system.time(rlt <- lapply(x,ctest,nop=100))





