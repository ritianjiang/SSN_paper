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

library(readr)
info<-read_tsv("../../../Dataset/GTEx/GTEx_2017-06-05_v8_sampleinfo.txt")
info<-info[,c(1,7)]
All_sample_edge_change<-merge(All_sample_edge_change,info,
                              by.x=1,by.y=1)
All_sample_edge_change$SMTSD<-str_replace_all(All_sample_edge_change$SMTSD,
                                              pattern = " - ",replacement = "_")
All_sample_edge_change$SMTSD<-str_replace_all(All_sample_edge_change$SMTSD,
                                              pattern = "\\)",replacement = "")
All_sample_edge_change$SMTSD<-str_replace_all(All_sample_edge_change$SMTSD,
                                              pattern = " \\(",replacement = "_")
All_sample_edge_change$SMTSD<-str_replace_all(All_sample_edge_change$SMTSD,
                                              pattern = " ",replacement = "_")
intersect(paste0(unique(All_sample_edge_change$SMTSD),"_covariates.csv"),
          list.files("../PEER/PEER/"))
All_sample_edge_change$PEER_file<-paste0(All_sample_edge_change$SMTSD,
                                         "_covariates.csv")

peer_file<-data.frame()
files<-intersect(paste0(unique(All_sample_edge_change$SMTSD),"_covariates.csv"),
                 list.files("../PEER/PEER/"))

age_cor_res<-data.frame(tissue = files, beta=0, p=-1)

for(i in 1:27){
  tmp_peer<-read.csv(paste0("../PEER/PEER/",files[i]))
  tmp_loss<-All_sample_edge_change[All_sample_edge_change$PEER_file == files[i] &
                                     All_sample_edge_change$type == "loss",]
  tmp<-merge(tmp_loss,tmp_peer,by.x="Subject.ID", by.y="SUBJID")
  tmp_fit<-lm(data = tmp,counts~Age+V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18)
  fit_sum<-summary(tmp_fit)
  fit_sum$coefficients[2,1]->age_cor_res[i,]$beta
  fit_sum$coefficients[2,4]->age_cor_res[i,]$p
}

age_cor_res$sig<-ifelse(test = age_cor_res$p<0.05,yes = "Sig.",no = "N.S.")

library(ggrepel)
age_cor_res$v<--log10(age_cor_res$p)
age_cor_res$tissue<-str_replace(age_cor_res$tissue,pattern = "_cova.*",replacement = "")
ggscatter(age_cor_res,x = "beta",y="v",
          color="sig",palette = c("grey","#FC4E07"),
          xlab = "Estimated coefficient beta with age",
          ylab = "-log10(P.value)",
          label = "tissue",repel=T,
          label.select = c("Esophagus_Mucosa","Whole_Blood","Lung"))+
  xlim(c(-240,240))
