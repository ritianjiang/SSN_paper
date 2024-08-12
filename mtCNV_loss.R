setwd("/home/owht/Data/KIZ/data/pCPM/Result/Summary/Loss_sig_metrics")
library(ggplot2)
library(stringr)
library(magrittr)
library(tidyverse)
library(ggpubr)

load("All_sample_edge_change.RData")
mtCN<-read.table("../mtCN+GTEx_pnas.2402291121.txt",sep="\t",header=T)
mtCN$tissue<-str_replace(mtCN$Tissue.Site.Detail,pattern = " - ",replacement = "#")
mtCN$tissue<-str_split(mtCN$tissue,pattern = "#",simplify = T)[,1]
mtCN$ind<-str_split(mtCN$Sample.ID.for.data.sharing.and.public.release,
                    pattern = "-",simplify = T)[,1:2] %>% apply(1,str_c,collapse="-")

All_sample_edge_change$ind<-str_split(All_sample_edge_change$ID,
                                      pattern = "-",simplify = T)[,1:2] %>% 
  apply(1,str_c,collapse="-")

final<-merge(All_sample_edge_change,mtCN,by.x=c("ind","tissue"),by.y=c("ind","tissue"))

mt_res<-data.frame(tissue = unique(final$tissue),cor=0,p=-1)
for(i in 1:20){
  tmp_loss<-final[final$tissue == mt_res[i,]$tissue,]
  fit_loss<-cor.test(tmp_loss$counts,tmp_loss$mtCN.Per.Cell,method = "spearman");
  mt_res[i,]$cor<-fit_loss$estimate;
  mt_res[i,]$p<-fit_loss$p.value
}
mt_res$sig<-ifelse(test = mt_res$p<0.05,yes = "Sig.",no = "N.S.")

ggbarplot(mt_res,x = "tissue",y="cor",
          fill="sig",palette = c("grey","#FC4E07"),
          ylab = "Correlation with MT copy number")+
  coord_flip()
