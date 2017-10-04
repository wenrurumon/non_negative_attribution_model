
rm(list=ls())
library(SIMLR)
library(data.table)
library(igraph)

#setwd('e:\\uthealth\\deconv')
setwd('/zfsauton/project/highmark/aska/sc2')
load("data4simlr.rda")
#browseVignettes("SIMLR")

length(unique(data$clust$celltype))
system.time(test <- SIMLR(X = data$ref, c = 17, cores.ratio = 0))
save(test,file='test_rlt.rda')

table(test$y$cluster,data$clust$celltype[1:length(test$y$cluster)])
