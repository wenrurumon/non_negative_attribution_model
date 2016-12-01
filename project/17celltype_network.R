
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
geneinpath <- geneinpath[names(geneinpath)%in%pathlist[grep('Signal |Neurodegenerative',pathlist[,1]),2]]
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
test <- function(Y.model,max_parent=3){
  Y=do.call(cbind,lapply(Y.model,function(x)x[[1]]))
  Y.group=rep(1:length(Y.model),sapply(Y.model,function(x){ncol(x[[1]])}))
  Y.prop=unlist(lapply(Y.model,function(x){diff(c(0,x[[2]]))}))
  lambda1=.4
  lambda2=.3;
  times=100;
  stability=.8
  Y.sem <- sem_grplasso2(Y,Y.group,Y.prop,lambda1,lambda2,times,stability)[[2]]>=0.8
  Y.cnif <- CNIF(data=do.call(cbind,lapply(Y.model,function(x){x[[1]]})),init.adj=Y.sem,max_parent=max_parent)
  gc()
  Y.out <- mat.sds(adj.group(Y.cnif*unlist(lapply(Y.model,function(x){diff(c(0,x[[2]]))})),rep(1:length(Y.model),sapply(Y.model,function(x){ncol(x[[1]])}))))
  dimnames(Y.out) <- list(names(Y.model),names(Y.model))
  plotnet(Y.out>0,'directed')
  return(list(data=list(Y=Y,Y.group=Y.group,Y.prop=Y.prop,Y.cnif=Y.cnif),adj=Y.out))
}
test2 <- function(x){
  temp <- try(test(x))
  if(!is.list(temp)){temp <- try(test(x,max_parent=2))}
  return(temp)
}

j <- 0
Y.models <- lapply(
  expr,function(Y){
    print(j<<-j+1)
    Y.sig <- Y[names(Y)%in%pathlist[grep('Sig',pathlist[,1]),2]]
    for(i in 1:length(Y.sig)){colnames(Y.sig[[i]]$scores) <- paste0(gsub(" ",'_',names(Y.sig)[i]),1:ncol(Y.sig[[i]]$scores))}
    Y.d <- Y[names(Y)%in%pathlist[grep('Neurodegenerative',pathlist[,1]),2]]
    for(i in 1:length(Y.d)){colnames(Y.d[[i]]$scores) <- paste0(gsub(" ",'_',names(Y.d)[i]),1:ncol(Y.d[[i]]$scores))}
    signet <- test2(Y.sig)
    dnet <- test2(Y.d)
    sigdnet <- sparse_2sem(Y=dnet[[1]]$Y,X=signet[[1]]$Y,lambda=0.1,Y.fixed=dnet[[1]]$Y.cnif,times=100)[[1]]
    rlt <- matrix(0,ncol(sigdnet),ncol(sigdnet))
    rlt[1:nrow(sigdnet),] <- sigdnet
    rlt[-1:-ncol(dnet[[1]]$Y),-1:-ncol(dnet[[1]]$Y)] <- signet$data$Y.cnif
    newprop <- c(dnet$data$Y.prop,signet$data$Y.prop)
    newgroup <- c(dnet$data$Y.group,signet$data$Y.group+max(dnet$data$Y.group))
    return(list(data=cbind(dnet[[1]]$Y,signet[[1]]$Y),adj=rlt,Y.group=newgroup,Y.prop=newprop))
  })
dname <- c(names(expr[[1]])[names(expr[[1]])%in%pathlist[grep('Neurodegenerative',pathlist[,1]),2]],
           names(expr[[1]])[names(expr[[1]])%in%pathlist[grep('Sig',pathlist[,1]),2]])
save(Y.models,file='Y_models.rda')

#########################
# Phenotype
#########################

load('Y_models.rda')
load('phenonet.rda')
phe[260,] <- colMeans(phe[260:261,]); phe <- phe[-261,]; rownames(phe) <- unique(pheid)
j <- 0
YX.models <- lapply(Y.models,function(x){
  print(j<<-j+1)
  phe_model <- sparse_2sem(Y=phe,X=x$data,Y.fixed=sem_final_phe[[1]],lambda=.01,times=100)
  par(mfrow=c(2,1))
  rlt <- matrix(0,ncol(phe_model$eq_stability),ncol(phe_model$eq_stability))
  rlt[1:nrow(phe_model$eq_stability),] <- phe_model$eq_stability>=.8
  rlt[-1:-nrow(phe_model$eq_stability),-1:-nrow(phe_model$eq_stability)] <- x$adj
  rlt <- list(adj=rlt,group=c(1:20,x$Y.group+20),prop=c(rep(1,20),x$Y.prop))
  rlt$adj2 = mat.sds(adj.group(rlt$adj*rlt$prop,rlt$group))
  dimnames(rlt$adj2) <- list(c(colnames(phe),dname),c(colnames(phe),dname))
  plot(colSums(phe_model$eq_stability>=.8));plotnet(rlt$adj2>.1,'directed')
  return(rlt)
})
save(Y.models,YX.models,file='Y_models.rda')

#######################
# Summary
#######################

out <- YX.models
graph2adj <- function(g){
  temp <- matrix(0,0,4)
  rownames(g)[1:20] <- paste0('pheno::',rownames(g)[1:20])
  colnames(g)[1:20] <- paste0('pheno::',colnames(g)[1:20])
  for(i in 1:nrow(g)){
    for(j in 1:ncol(g)){
      if(g[i,j]>0){
        tempi <- c(paste(rownames(g)[i],colnames(g)[j],sep=','),rownames(g)[i],colnames(g)[j],g[i,j])
        temp <- rbind(temp,tempi)
      }
    }
  }
  return(temp)
}
out <- lapply(out,function(x){graph2adj(x$adj2)})
key <- unique(paste(unlist(lapply(out,function(x){x[,1]}))))

write.csv(cbind(do.call(rbind,strsplit(key,',')),sapply(out,function(out){as.numeric(out[match(key,out[,1]),4])})),'temp.csv')
