setwd("/home/owht/Data/KIZ/data/pCPM/Result/Summary/Loss_sig_metrics")
library(ggplot2)
library(stringr)
library(magrittr)
library(tidyverse)
library(ggpubr)

load("All_sample_edge_change.RData")
#All_sample_edge_change$geneage<-0
stem<-read.csv("../Stemness/aging-16-205717-s004.csv")
stem_res<-data.frame(tissue = unique(All_sample_edge_change$tissue),
                        cor=0,p=-1)
All_sample_edge_change<-merge(All_sample_edge_change,
                              stem,by.x=1,by.y=1)


for(i in 1:26){
  tmp_loss<-All_sample_edge_change[All_sample_edge_change$tissue == stem_res[i,]$tissue
                                   & All_sample_edge_change$type == "loss",]

  fit_loss<-cor.test(tmp_loss$counts,
                     tmp_loss$Stemness.Score,method = "spearman");
  stem_res[i,]$cor<-fit_loss$estimate;
  stem_res[i,]$p<-fit_loss$p.value
}
stem_res$padj<-p.adjust(stem_res$p,method = "BH")

stem_res$sig<-ifelse(test = stem_res$padj<0.05,yes = "Sig.",no = "N.S.")

library(ggrepel)
ggbarplot(stem_res,x = "tissue",y="cor",
          fill="sig",palette = c("grey","#FC4E07"),
          ylab = "Correlation with stemness score")+
  coord_flip()

