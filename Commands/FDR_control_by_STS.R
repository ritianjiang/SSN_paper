setwd("~/data/CPM_methods_2/Networks/All_p1/")
library(qvalue)
library(mutoss)
library(data.table)
tissues<-list.dirs(recursive=F)

qRes<-data.frame(sample=list.files(recursive = T,full.names = T,include.dirs = F),
                 q005=0,p001=0)
for (index in 1:nrow(qRes)){
  x<-fread(qRes$sample[index],header=T)
  colnames(x)[4]<-"p_value"
  qvs<-adaptiveSTS(x$p_value,alpha = 0.05,lambda = 0.5,silent = T)
  
  #y<-fread(paste0("~/data/CPM_methods/Networks/",files[i]),header=T)
  sum(qvs$adjPValues<0.05 & x$p_value<0.001)->nums
  # length(intersect(paste0(x[qvs$qvalue<0.05,]$Gene1,"_",x[qvs$qvalue<0.05,]$Gene2),
  #           paste0(y$Gene1,"_",y$Gene2)))->nums
  qRes[index,]$p001<-sum(x$p_value<0.001); 
  qRes[index,]$q005<-nums;
}
all(qRes$q005 == qRes$p001)

