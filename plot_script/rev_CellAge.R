setwd("/home/owht/Data/KIZ/data/pCPM/Result/Summary/Loss_sig_metrics")
library(ggplot2)
library(stringr)
library(magrittr)
library(tidyverse)
library(ggpubr)

load("All_sample_edge_change.RData")
#All_sample_edge_change$geneage<-0
files<-list.files("./CellAge/")
geneage_res<-data.frame(tissue = str_split(files,
                                           pattern = "cell",
                                           simplify = T)[,1],
                        cor=0,p=-1)

for(i in 1:26){
  load(paste0("./CellAge/",files[i]))
  tmp_loss<-All_sample_edge_change[All_sample_edge_change$tissue == geneage_res[i,]$tissue
                                   & All_sample_edge_change$type == "loss",]
  Hsig_es<-t(Hsig_es) %>% as.data.frame
  Hsig_es$id<-rownames(Hsig_es)
  tmp_loss<-merge(tmp_loss,Hsig_es,by.x="ID", by.y="id")
  fit_loss<-cor.test(tmp_loss$counts,tmp_loss$GENEAGE_DATABASE,
                     method = "spearman");
  geneage_res[i,]$cor<-fit_loss$estimate;
  geneage_res[i,]$p<-fit_loss$p.value
}

geneage_res$padj<-p.adjust(geneage_res$p,method = "BH")
geneage_res$sig<-ifelse(test = geneage_res$padj<0.05,yes = "Sig.",no = "N.S.")

ggbarplot(geneage_res,x = "tissue",y="cor",
          fill="sig",palette = c("grey","#FC4E07"),
          ylab = "Correlation with aging gene expression (GenAge)")+
  coord_flip()
