setwd("/data/user_data/wanght/CPM_methods/xCell/")
library(magrittr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(psych)
library(reshape)
source("CIBERSORT_func.R")
LM22<-read.table("LM22.txt",header = T,row.names = 1,sep="\t")
tissues<-list.dirs("../Exps/")
tissues<-str_split(tissues[2:27],pattern = "//",simplify = T)[,2]
load("../SumExp_GTExv8_tpm.RData")
se_CIBERSORT<-se_GTExv8_tpm[rownames(se_GTExv8_tpm) %in% rownames(LM22),]
rm(se_GTExv8_tpm)

if(!dir.exists("CIBERSORT/")){
  dir.create("CIBERSORT/")
}

####### STEP 1 CIBERSORT ###############
#LM22<-read.table("LM22.txt",header = T,row.names = 1,sep="\t")

for(i in 5:26){
  cat(tissues[i]);cat("\n")
  sample_id<-read.table(paste0("../Exps/",tissues[i],"/samples.txt"),header = F,
                        nrows = 1)
  sample_id<-t(sample_id)[2:ncol(sample_id)]
  
  exp_tissue<-assay(se_CIBERSORT[,sample_id])
  exp_tissue<-as.matrix(exp_tissue)
  CIBERSORT_res<-CIBERSORT(sig_matrix = LM22,
                           mixture_file = exp_tissue,
                           filename = paste0("./CIBERSORT/",tissues[i]))
}

###### STEP 2 Corr#######################

CIBERSORT_loss_cor<-matrix(nrow=22,ncol=26,data=0)
CIBERSORT_loss_p<-matrix(nrow=22,ncol=26,data=0)
colnames(CIBERSORT_loss_cor)<-tissues
rownames(CIBERSORT_loss_cor)<-colnames(CIBERSORT_tissue)[1:22]
colnames(CIBERSORT_loss_p)<-tissues
rownames(CIBERSORT_loss_p)<-rownames(CIBERSORT_loss_cor)

#Load TISSUE
for(i in 1:26){
  cat(tissues[i]);cat("\n")
  load(paste0("../Result/",tissues[i],"/Result.RData"))
  sample_id<-read.table(paste0("../Exps/",tissues[i],"/samples.txt"),header = F,
                        nrows = 1)
  sample_id<-t(sample_id)[2:ncol(sample_id)]
  
  CIBERSORT_tissue<-read.table(paste0("./CIBERSORT/",tissues[i]),row.names = 1,header = T)
  CIBERSORT_tissue<-CIBERSORT_tissue[,1:22]
  
  apply(CIBERSORT_tissue,2,cor,y=as.numeric(Edges_loss),method="spearman")->CIBERSORT_loss_cor[,i]
  apply(CIBERSORT_tissue,2,function(x){return(cor.test(x,Edges_loss,method="spearman")$p.value)})->CIBERSORT_loss_p[,i]
  
}
CIBERSORT_loss_p[is.na(CIBERSORT_loss_p)]<-1
CIBERSORT_loss_cor[is.na(CIBERSORT_loss_cor)]<-0

myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  if(is.na(pval))
    stars = ""
  stars
}


heatmap <- melt(CIBERSORT_loss_cor) %>% rename(replace=c("X1"="Cell","X2"="Tissue",
                                                     "value"="cor")) %>%
  mutate(pvalue=melt(CIBERSORT_loss_p)[,3]) %>%
  mutate(signif = sapply(CIBERSORT_loss_p, function(x) myfun(x)))

ggplot(heatmap,aes(Tissue,Cell,col=cor))+
  geom_tile(color="grey90",fill="white",size=1)+
  geom_point(aes(size = abs(cor)),shape=20) + 
  geom_text(aes(label=signif),size=6,color="black",
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_gradient2(low = "blue",mid = "white",high="red")+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(family="Roboto"),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),axis.text.x.bottom = element_text(angle = 90),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid")) +
  scale_size(range=c(1,10),guide=NULL)+
  guides(color = guide_colorbar(direction = "vertical",
                                reverse = F,barwidth = unit(.5, "cm"),
                                barheight = unit(3, "cm")))
save(heatmap,file = "CIBERSORT_res.RData")
