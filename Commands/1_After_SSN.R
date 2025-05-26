library(stringr,quietly = T,verbose = F)
library(dplyr,verbose = F,quietly = T)
args<-commandArgs()
TISSUE<-args[6]
cat(paste0("Current Tissue is ",TISSUE,"\n"))

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

#Load individual files
ind_info<-get_ind_files("/home/wanght/data/CPM_methods/Dataset/GTEx_2017-06-05_v8_indinfo.txt")
cat("Loading individual sampleinfo files...... \n")

#Load network files
fileS<-list.files(paste0("./Networks/",TISSUE))
filed<-paste0("./Networks/",TISSUE)
refs<-read.table(paste0("Exps/",TISSUE,"/reference.txt"),header = T,row.names = 1)

inds<-fileS %>%
  str_split(pattern = "-|_",simplify = T) %>%
  .[,2:3] %>% apply(1,FUN=str_c,collapse="-")
Edges_gain<-rep(0,length(fileS))
Edges_loss<-rep(0,length(fileS))

cat("Loading SSN files...... \n")
Loss_list<-list()
cor_each<-function(tsx){cor(refs[tmp[tsx,]$Gene1,] %>% as.numeric(),
                            refs[tmp[tsx,]$Gene2,] %>% as.numeric())}
for(f in 1:length(fileS)){
  tmp<-read.table(paste0(filed,"/",fileS[f]),sep="\t",header = T)
  tmp$ref_cor<-0;tmp$trend<-"Gain"
  tmp$ref_cor<-sapply(1:nrow(tmp),cor_each,simplify="array")
  tmp[tmp$deltaPCC * tmp$ref_cor<0,]$trend<-"Loss"
  Edges_gain[f]<-nrow(tmp[tmp$trend == "Gain",])
  Edges_loss[f]<-nrow(tmp[tmp$trend == "Loss",])
  Loss_list<-c(Loss_list,list(paste0(tmp[tmp$trend == "Loss",1],"_",tmp[tmp$trend == "Loss",2])))
}

cat("Perform SSN-age examination......\n")
Loss_tissue<-unique(unlist(Loss_list))
Loss_tissue_mat<-matrix(0,nrow=length(Loss_tissue), #Row means interaction
                       ncol=length(Edges_gain)) #col mean sample
for(i in 1:nrow(Loss_tissue_mat)){
  tx<-Loss_tissue[i]
  Loss_tissue_mat[i,]<-as.integer(sapply(Loss_list,function(x){return(tx %in% x)},simplify="array"))
}
tissue_results<-data.frame(edges=Loss_tissue,tau=0,p=-1)
ages<-ind_info[inds,]$AGE
for(i in 1:nrow(tissue_results)){
  fit<-cor.test(Loss_tissue_mat[i,],
                ages,method = "kendall")
  fit$estimate->tissue_results[i,]$tau
  fit$p.value->tissue_results[i,]$p
}
names(ages)<-inds

cat(paste0("Saving results for ",TISSUE,"......\n"))

try(expr = {dir.create(paste0("./Result/",TISSUE))},
    silent = T)
ResD<-paste0("./Result/",TISSUE,"/Result.RData")
save(tissue_results,Edges_gain,Edges_loss,Loss_list,Loss_tissue_mat,ages,file=ResD)