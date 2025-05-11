#Retain q-value < 0.05 

setwd("~/data/CPM_methods_2/Networks/All_p1/")
library(qvalue)
library(data.table)
tissues<-list.dirs(recursive=F)

qRes<-data.table(sample=rep("Nx",15394),
                 q005=rep(0,15394),
                 p001=rep(0,15394),tissue=rep("Nx",12000))
index=1

for (d in tissues){
    files<-list.files(path=d,full.names=T)
    for(i in 1:length(files)){
        x<-fread(files[i],header=T)
        colnames(x)[4]<-"p_value"
        qvs<-qvalue(x$p_value,lambda=seq(0.01,0.45,0.02))
        
        #y<-fread(paste0("~/data/CPM_methods/Networks/",files[i]),header=T)
        sum(qvs$qvalues<0.05 & x$p_value<0.001)->nums
        # length(intersect(paste0(x[qvs$qvalue<0.05,]$Gene1,"_",x[qvs$qvalue<0.05,]$Gene2),
        #           paste0(y$Gene1,"_",y$Gene2)))->nums
        qRes[index,]$sample<-files[i];qRes[index,]$tissue<-d
        qRes[index,]$p001<-sum(x$p_value<0.001); 
        qRes[index,]$q005<-sum(qvs$qvalues<0.05 & x$p_value<0.001);
        index=index+1
        fwrite(x[qvs$qvalues<0.05 & x$p_value<0.001,],col.names = T,quote=F,sep="\t",
                file=paste0("../FDRed/",files[i]))
        if(index %in% c(100,500,1000,2000,4000,6000,8000,10000)){cat(index);cat(" ")}
    }
}
