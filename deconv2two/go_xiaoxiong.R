rm(list=ls())
setwd('E:\\uthealth\\deconv\\20171009\\onserver')
source("contiANM_1.3.R")
source("contidANM.R")
source("discreteANM_1.1.R")
load("mdata4causal.rda" )

causal.ad <- function(x){
  r1 <- c(unlist(permcontidANM(scale(x),ad,100))[1:2],
          P_ttest=t.test(scale(x)~ad)$p.value)
  x <- x * mprop
  r2 <- c(unlist(permcontidANM(scale(x),ad,100))[1:2],
          P_ttest=t.test(scale(x)~ad)$p.value)
  c(r1,r2)
}

rlt <- sapply(1:nrow(raw_exp),function(i){
  print(i)
  x <- raw_exp[,i]
  causal.ad(x)
})

save(rlt,file='rlt1.rda')
