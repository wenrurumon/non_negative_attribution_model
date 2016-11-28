
rm(list=ls())
setwd('C:/Users/zhu2/Documents/deconv/sj2016')
source('source.R')
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)
########################
#Data Processing
########################
setwd('C:/Users/zhu2/Documents/trail/sample/')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
setwd('C:/Users/zhu2/Documents/deconv/sj2016')
load("decomp_forkl.rda")
load("pathlist.rda")
load("geneinpath.rda")
geneinpath <- geneinpath[names(geneinpath)%in%pathlist[grep('Sig|Neurodegenerative',pathlist[,1]),2]]
expr <- lapply(decomp_raw,function(x){
  lapply(geneinpath,function(genes){
    x <- x[,colnames(x)%in%genes,drop=F]
    x <- qpca(x,which(pca(x)$prop>=0.9)[1])
    list(scores = x$X[,1:which(x$prop>=.9)[1],drop=F],
         prop = x$prop[1:which(x$prop>=.9)[1]])
    })
  })
########################
#Sample network for Y
########################
# 
# Y <- expr[[1]]
# Y.sig <- Y[names(Y)%in%pathlist[grep('Sig',pathlist[,1]),2]]
# Y.d <- Y[names(Y)%in%pathlist[grep('Neurodegenerative',pathlist[,1]),2]]
# 
# Y.model <- Y.d
test <- function(Y.model){
  print('model')
  Y.sem <- sem_grplasso2(
    Y=do.call(cbind,lapply(Y.model,function(x)x[[1]])),
    Y.group=rep(1:length(Y.model),sapply(Y.model,function(x){ncol(x[[1]])})),
    Y.prop=unlist(lapply(Y.model,function(x){diff(c(0,x[[2]]))})),
    lambda1=.3,lambda2=.1,times=10,stability=.8)[[2]]>=.8
  Y.cnif <- CNIF(data=do.call(cbind,lapply(Y.model,function(x){x[[1]]})),init.adj=Y.sem,max_parent=3)
  gc()
  Y.out <- mat.sds(adj.group(Y.cnif*unlist(lapply(Y.model,function(x){diff(c(0,x[[2]]))})),rep(1:length(Y.model),sapply(Y.model,function(x){ncol(x[[1]])}))))
  dimnames(Y.out) <- list(names(Y.model),names(Y.model))
  plotnet(Y.out>0,'directed')
  return(list(data=Y.model,adj=Y.out))
}
j <- 0
Y.models_sig <- lapply(
  expr,function(Y){
    print(j<<-j+1)
    Y.sig <- Y[names(Y)%in%pathlist[grep('Sig',pathlist[,1]),2]]
    Y.d <- Y[names(Y)%in%pathlist[grep('Neurodegenerative',pathlist[,1]),2]]
    signet <- test(Y.sig)
    dnet <- test(Y.d)
    list(signet=test(Y.sig),dnet=test(Y.d))
})
