setwd("/data/user_data/wanght/CPM_methods/xCell/")
library(magrittr)
library(e1071)
library(parallel)
library(preprocessCore)
library(data.table)
library(tidyverse)
library(ggplot2)
library(psych)
library(reshape)
source("CIBERSORT_1_04.R")
LM22<-read.table("LM22.txt",header = T,row.names = 1,sep="\t")
tissues<-list.dirs("../Exps/")
tissues<-str_split(tissues[2:27],pattern = "//",simplify = T)[,2]
load("../SumExp_GTExv8_tpm.RData")
se_CIBERSORT<-se_GTExv8_tpm[rownames(se_GTExv8_tpm) %in% rownames(LM22),]
rm(se_GTExv8_tpm)

if(!dir.exists("CIBERSORT_abs/")){
  dir.create("CIBERSORT_abs/")
}

####### STEP 1 CIBERSORT ###############
#LM22<-read.table("LM22.txt",header = T,row.names = 1,sep="\t")

for(i in 1:26){
  cat(tissues[i]);cat("\n")
  sample_id<-read.table(paste0("../Exps/",tissues[i],"/samples.txt"),header = F,
                        nrows = 1)
  sample_id<-t(sample_id)[2:ncol(sample_id)]
  
  exp_tissue<-assay(se_CIBERSORT[,sample_id])
  tmp<-cbind(data.frame(`Gene symbol`=rownames(exp_tissue)),exp_tissue)
  colnames(tmp)[1]<-"Gene symbol"
  # write.table(tmp,row.names = F,col.names = T,quote = F,sep="\t",
  #             file="tmp.txt")
  CIBERSORT_res<-CIBERSORT(sig_matrix = LM22,
                           mixture_file = exp_tissue,
                           absolute = TRUE,perm = 1,abs_method = "no.sumto1")
  write.table(CIBERSORT_res,file=paste0("./CIBERSORT_abs/",tissues[i]))
}

