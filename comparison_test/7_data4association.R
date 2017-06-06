
rm(list=ls())
setwd('/home/zhu/rushdata/expression_clustering')
load('disease.rda')
setwd('/home/zhu/deconv')
load('data/rlt4compare.rda')

pen <- function(x,lambda=0){
  (x-lambda) * ((x-lambda)>0)
}

#rushdata

pid <- rownames(disease)
pid <- pid[pid%in%rownames(Ynmf.rush)]
dr <- disease[rownames(disease)%in%pid,4]
Yrn <- Ynmf.rush[match(pid,rownames(Ynmf.rush)),]
Yrs <- Ystf.rush[match(pid,rownames(Ystf.rush)),]

#adnidata

adnid <- read.csv("/home/zhu/deconv/data/ADNIMERGE.csv")
adnid <- unique(adnid[,c(2,8)])
pid <- paste(adnid[,1])
da <- 1-as.numeric(adnid[,2]=="CN")[match(rownames(Ynmf.adni),pid)]
Yan <- Ynmf.adni
Yas <- Ystf.adni

save(pen,Yan,Yas,Yrn,Yrs,dr,da,file='data/rlt4compare2.rda')
