##############################
# Data Processing
##############################

rm(list=ls())

rf <- 1
th <- 0.01
nop <- 100

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

###################################
# Gene Selection
###################################

# genesel <- colnames(raw_exp)[1:10]

load("E:/uthealth/deconv/20171009/rlt_m/rlt_r2.rda")
rlt_bk <- rlt
rlt <- rlt[,c(2,5)]
rlt <- apply(rlt,1,min)

genesel <- names(rlt[rlt<=th])
rlt_bk <- rlt_bk[rlt<=th,]

length(genesel)
expr_sel <- raw_exp[,colnames(raw_exp)%in%genesel,drop=F]

###################################
# Causal to 5
###################################

x <- apply(expr_sel,2,function(x){
  list(x * deconv.rlt)
})
x <- do.call(c,x)

i <- 0
ctest <- function(x,y=ad,nop=100){
  print(i<<-i+1)
  rlt <- apply(x,2,function(x){
    rlt1 <- unlist(permcontidANM(x,y,nop))[1:2]
    rlt2 <- try(t.test(x~y)$p.value)
    if(!is.numeric(rlt2)){rlt2 <- NA}
    rlt <- c(rlt1,p_ttest=rlt2)
  })  
  rlt
}

system.time(rlt <- lapply(x,ctest,nop=nop))


# save(rlt,file='gene1227_5cell.rda')

###################################
# output causal to 5
###################################

x <- rlt[[1]]
sel <- sapply(rlt,function(x){min(x[2,])})
names(rlt[sel<=0.05])
rlt[names(rlt)%in%c('CRNDE','KCNC3')]



