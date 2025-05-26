library(stringr,quietly = T,verbose = F)
library(dplyr,verbose = F,quietly = T)
library(igraph)
library(magrittr)

# args<-commandArgs()
# TISSUE<-args[6]
# cat(paste0("Current Tissue is ",TISSUE,"\n"))

##FUNCTIONS
get_ind_files<-function(filepath){ #load sampleinfo files
  ind_info<-read.table(filepath,
                       fill = T, # From GTEx Protal
                       header = T,row.names = 1)
  ind_info$AGE<-ind_info$AGE %>%
    str_split(,pattern = "-",simplify = T) %>%
    .[,1] %>% as.numeric() %>% `+`(5)        #convert range into int, e.g. 20-29 ==> 25
  ind_info[ind_info$SEX == 2,]$SEX<-"F"
  ind_info[ind_info$SEX == "1",]$SEX<-"M"
  ind_info$SEX<-factor(ind_info$SEX)
  return(ind_info)
}

##PROCUDURES
setwd("/data/user_data/wanght/CPM_methods")

tissues<-list.dirs("./Exps/")
tissues<-str_split(tissues[2:27],pattern = "//",simplify = T)[,2]

PPIs<-read.table("9606.protein.links.v11.5.PPI700.symbols.txt",header = T)
PPIs$Int<-paste0(PPIs$Symbol1,"_",PPIs$Symbol2)
load("SumExp_GTExv8_tpm.filtered_0.3_100_PPI700.RData")

#Load TISSUE
for(TISSUE in tissues){
  cat(TISSUE);cat("\n")
  load(paste0("./Result/",TISSUE,"/Result.RData"))
  sample_id<-read.table(paste0("./Exps/",TISSUE,"/samples.txt"),header = F,
                        nrows = 1)
  sample_id<-t(sample_id)[2:ncol(sample_id)]
  
  #retain_PPIs<-PPIs[!PPIs$Int %in% tissue_results$edges,]
  tissue_results_sig<-tissue_results[tissue_results$p<0.001 & tissue_results$tau>0,]
  #Summarys[i,]$loss_sig_nums<-nrow(tissue_results_sig)
  #cat(i);cat(" ")
  tissue_results_sig$Gene1<-str_split(tissue_results_sig$edges,pattern = "_",simplify = T)[,1]
  tissue_results_sig$Gene2<-str_split(tissue_results_sig$edges,pattern = "_",simplify = T)[,2]
  Genes<-data.frame(c(tissue_results_sig$Gene1,tissue_results_sig$Gene2)) %>% distinct
  colnames(Genes)<-c("Gene")
  edges<-tissue_results_sig[seq(1,nrow(tissue_results_sig),2),c(4,5,2)] %>% distinct()
  net_pc<-graph_from_data_frame(
    d=edges,vertices=Genes,
    directed=F)
  deg<-degree(net_pc,mode="all")
  deg<-sort(deg,decreasing = T)
  
  hub_gene<-names(deg)[1]
  hub_neib<-neighborhood(graph = net_pc,order = 1,nodes = hub_gene)[[1]] %>% names
  hubs_exp<-se_GTExv8_tpm.filtered[hub_neib,sample_id]
  apply(assay(hubs_exp),1,function(x){cor(x,ages)})->hub_cor
  if(median(hub_cor)>0){
    MR_genes<-names(hub_cor[hub_cor>0])
  }else{MR_genes<-names(hub_cor[hub_cor<0])}
  cat("Pos:");cat(sum(hub_cor>0));cat(" Neg:");cat(sum(hub_cor<0));cat("\n")
  
  save(MR_genes,file=paste0("./MR/",TISSUE,".RData"))
}


