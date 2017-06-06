
rm(list=ls())
setwd("/home/zhu/deconv")
load('data/rlt4compare2.rda')
load('data/rlt4compare3.rda')

asso <- function(d,y,pen=0){
  y <- (y-pen) * ((y-pen)>0)
  apply(y,2,function(x){
    t.test(x~d)$p.value
  })
}

assob <- function(d,ys,yn){
  rlt = list()
  rlt[[1]] <- asso(d,ys)
  rlt[[2]] <- asso(d,yn)
  for(i in 3:5){
    rlt[[i]] <- asso(d,yn,pen=(0.05)*(i-2))
  }
  do.call(cbind,rlt)
}

testr <- assob(dr,Yrn,Yrs)
testa <- assob(da,Yan,Yas)
testr2 <- assob(dr,Yrn2,Yrs2)
testa2 <- assob(da,Yan2,Yas2)
