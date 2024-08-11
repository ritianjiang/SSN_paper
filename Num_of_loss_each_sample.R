setwd("data/CPM_methods/Summary/")
load("Loss_sig_sample_nums.RData")
library(magrittr)
library(stringr)
setwd("/home/wanght/data/CPM_methods")
ind_info<-read.table("./Dataset/GTEx_2017-06-05_v8_indinfo.txt",
                     fill = T, # From GTEx Protal
                     header = T,row.names = 1)
ind_info$AGE<-ind_info$AGE %>%
  str_split(,pattern = "-",simplify = T) %>%
  .[,1] %>% as.numeric() %>% `+`(5)        #convert range into int, e.g. 20-29 ==> 25
ind_info[ind_info$SEX == 2,]$SEX<-"F"
ind_info[ind_info$SEX == "1",]$SEX<-"M"
ind_info$SEX<-factor(ind_info$SEX)

Dirs<-list.dirs("Result/")
library(ggplot2)
Summarys<-data.frame(Tissues = str_split(Dirs[2:27],pattern = "//",simplify = T)[,2],
                     gain_tau=0,gain_p=0,
                     loss_tau=0,loss_p=0)

All_sample_edge_change<-data.frame()

for(i in 1:26){
    load(paste0(Dirs[i+1],"/Result.RData"))
    fileS<-list.files(paste0("./Networks/",Summarys$Tissues[i]))
    sample_names<-str_split(fileS,pattern = "_|\\.",simplify = T)[,2]
    names(Edges_gain)<-sample_names;names(Edges_loss)<-sample_names
    tmp<-data.frame(ID=c(names(Edges_gain),names(Edges_loss)),
                    tissue = Summarys$Tissues[i],
                    type = c(rep("gain",length(Edges_gain)),
                             rep("loss",length(Edges_loss))),
                    counts = c(Edges_gain,Edges_loss))
    All_sample_edge_change<-rbind(All_sample_edge_change,tmp)
}
save(All_sample_edge_change,file="All_sample_edge_change.RData")
