######################################## STEP 0 Functions           #############################
load_GTEx_tpm<-function(loc){   #loc: File position
  filename<-paste0(loc,"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
  tpm<-fread(filename,
             skip = 2)
  tpm<-as.data.frame(tpm)
  tpm$Name<-str_split(tpm$Name,pattern = "\\.",simplify = T)[,1]
  tpm<-tpm %>% distinct(Name, .keep_all = T)
  tpm<-tpm %>% distinct(Description, .keep_all = T)
  rowData<-DataFrame(Name=tpm$Name,SYMBOL=tpm$Description)
  tpm<-tpm[,-1];rownames(tpm)<-tpm[,1];tpm<-tpm[,-1]
  return(list(tpm,rowData))
}

load_GTEx_sampleinfo<-function(loc,tpm){
  filename<-paste0(loc,"GTEx_2017-06-05_v8_sampleinfo.txt")
  coda<-read.table(filename,header = T,sep="\t")
  coda<-coda[coda$SAMPID %in% colnames(tpm),]
  
  ind_info<-read.table(paste0(loc,"GTEx_2017-06-05_v8_indinfo.txt"),fill = T, # From GTEx Protal
                       header = T)
  ind_info$AGE<-ind_info$AGE %>%
    str_split(,pattern = "-",simplify = T) %>%
    .[,1] %>% as.numeric() %>% `+`(5)        #convert range into int, e.g. 20-29 ==> 25
  ind_info[ind_info$SEX == 2,]$SEX<-"F"
  ind_info[ind_info$SEX == "1",]$SEX<-"M"
  ind_info$SEX<-factor(ind_info$SEX)
  
  coda$individuals<-coda$SAMPID %>% 
    str_split(pattern = "-",simplify = T) %>% 
    .[,1:2] %>% apply(1,FUN=str_c,collapse="-")
  coda<-merge(coda,ind_info,by.x="individuals",by.y=1)
  return(coda)
}

SE_filter<-function(data,min_sample_per_tissue=100,
                    ppi_file="none",
                    perc_TPM_lessthan1=0.3){
  cat("Remove tissues less than:");cat(min_sample_per_tissue);cat("\n")
  tissues<-table(colData(data)$SMTS) %>% #select only tissues > 100 samples
    as.data.frame() %$% .[which(Freq>min_sample_per_tissue),1] %>%
    as.character()
  data<-data[,data$SMTS %in% tissues]
  
  if(ppi_file!="none"){
    cat("Retain genes in PPi network:");cat(ppi_file);cat("\n")
    ppi<-read.table(ppi_file,header = T)
    data<-data[rownames(data) %in% unique(c(ppi$Symbol1,ppi$Symbol2)),]
  }
  
  cat("Remove TPM<1 more than proportion:");cat(perc_TPM_lessthan1);cat("\n")
  genes_filtered<-apply(assay(data),1,function(x){
    sum(x<1)< perc_TPM_lessthan1*(dim(data)[2])
  })
  data<-data[genes_filtered,]
  
  return(data)
}



######################################## STEP 1 Setup the environment ###########################
setwd("/home/wanght/data/CPM_methods")

library(stringr)
library(SummarizedExperiment)
library(data.table)
library(magrittr)
library(dplyr)
library(tidyr)
library(UpSetR)
#load TPM and sampleinfo
GTEx_LOCATION="/home/wanght/data/Evo_longe/Exp_GTEx/"
tpm_result<-load_GTEx_tpm(GTEx_LOCATION)
tpm<-tpm_result[[1]];rowDataX<-tpm_result[[2]];
rm(tpm_result)
coda<-load_GTEx_sampleinfo(GTEx_LOCATION,tpm)

#make SummarizedExperiment object
se_GTExv8_tpm<-SummarizedExperiment(assays = tpm,colData = DataFrame(coda),rowData = rowDataX)
save(se_GTExv8_tpm,file="SumExp_GTExv8_tpm.RData")
rm(coda,rowDataX,tpm);gc()
table(se_GTExv8_tpm$SMTS)

#Filter SummarizedExperiment object
se_GTExv8_tpm.filtered<-SE_filter(se_GTExv8_tpm,min_sample_per_tissue = 100,
                                  ppi_file = "9606.protein.links.v11.5.PPI700.symbols.txt",
                                  perc_TPM_lessthan1 = 0.3)
save(se_GTExv8_tpm.filtered,file="SumExp_GTExv8_tpm.filtered_0.3_100_PPI700.RData")
se_GTExv8_tpm.filtered<-SE_filter(se_GTExv8_tpm,min_sample_per_tissue = 100,
                                  ppi_file = "9606.protein.links.v11.5.PPI700.symbols.txt",
                                  perc_TPM_lessthan1 = 0.3)

table(colData(se_GTExv8_tpm.filtered)$AGE,colData(se_GTExv8_tpm.filtered)$SMTS)

indids<-colData(se_GTExv8_tpm.filtered)[,c("individuals","SMTS")]
indids$fills<-T
as.data.frame(table(indids$individuals,indids$SMTS)) %>% spread(Var2,Freq) ->ind_mat


######################################## STEP 2 Scale (or not) and export files ###########################
load("SumExp_GTExv8_tpm.filtered_0.3_100_PPI700.RData")

tissues<-unique(colData(se_GTExv8_tpm.filtered)$SMTS)

for(tis in tissues){
  try(expr = {dir.create(paste0("Exps/",str_replace(tis,pattern = " ",replacement = "_")))},
      silent = T)
  refs<-assay(se_GTExv8_tpm.filtered[,se_GTExv8_tpm.filtered$SMTS == tis &
                                      se_GTExv8_tpm.filtered$AGE == 25]) %>% as.matrix()
  refs<-asinh(refs)
  refs<-cbind(data.frame(ID=rownames(refs)),refs)
  samples<-assay(se_GTExv8_tpm.filtered[,se_GTExv8_tpm.filtered$SMTS == tis &
                                          se_GTExv8_tpm.filtered$AGE != 25]) %>% as.matrix()
  samples<-asinh(samples)
  samples<-cbind(data.frame(ID=rownames(samples)),samples)
    
  write.table(refs,file=paste0("Exps/",
                               str_replace(tis,pattern = " ",replacement = "_"),
                               "/reference.txt"),quote=F,col.names = T,row.names = F,sep="\t")
  write.table(samples,file=paste0("Exps/",
                               str_replace(tis,pattern = " ",replacement = "_"),
                               "/samples.txt"),quote=F,col.names = T,row.names = F,sep="\t")
}

bg<-read.table("9606.protein.links.v11.5.PPI700.symbols.txt",header=T)
bg<-bg[(bg$Symbol1 %in% rownames(se_GTExv8_tpm.filtered)) & (bg$Symbol2 %in% rownames(se_GTExv8_tpm.filtered)),]
write.table(bg[,c(4,5)],col.names = F,
            file="STRING.v11.5.PPI700.background.symbols.txt",
            row.names = F,quote=F,sep="\t")
