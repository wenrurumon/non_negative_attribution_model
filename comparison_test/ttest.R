
rm(list=ls())
setwd('/home/zhu/rushdata/expression_clustering')
load('disease.rda')
setwd('/home/zhu/deconv')
load('data/rlt4compare.rda')

#rushdata

pid <- rownames(disease)
pid <- pid[pid%in%rownames(Ynmf.rush)]
d <- disease[rownames(disease)%in%pid,4]
Y1 <- Ynmf.rush[match(pid,rownames(Ynmf.rush)),]
Y2 <- Ystf.rush[match(pid,rownames(Ystf.rush)),]
pen <- function(x,lambda=0){
  (x-lambda) * ((x-lambda)>0)
}
test1 <- apply(Y1,2,function(x) t.test(x~d)$p.value)
test2 <- apply(Y2,2,function(x) t.test(x~d)$p.value)
test3 <- apply(Y2,2,function(x) t.test(pen(x,lambda=0.05)~d)$p.value)
