setwd("/data/user_data/wanght/CPM_methods/xCell/")
library(magrittr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(psych)
library(reshape)

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

##PROCUDURES

tissues<-list.dirs("./Exps/")
tissues<-str_split(tissues[2:27],pattern = "//",simplify = T)[,2]
xCell<-fread("GTEx_Analysis_v8_xCell_scores_7_celltypes.txt.gz",)
xCell<-xCell[,2:17383]
xCell<-t(xCell)
colnames(xCell)<-c("Adipocytes","Epithelial_cells","Hepatocytes","Keratinocytes",
                   "Myocytes","Neurons","Neutrophils")

xCell_loss_cor<-matrix(nrow=7,ncol=26,data=0)
xCell_loss_p<-matrix(nrow=7,ncol=26,data=0)
colnames(xCell_loss_cor)<-tissues
rownames(xCell_loss_cor)<-colnames(xCell)
colnames(xCell_loss_p)<-tissues
rownames(xCell_loss_p)<-colnames(xCell)

#Load TISSUE
for(i in 1:26){
  cat(tissues[i]);cat("\n")
  load(paste0("../Result/",tissues[i],"/Result.RData"))
  sample_id<-read.table(paste0("../Exps/",tissues[i],"/samples.txt"),header = F,
                        nrows = 1)
  sample_id<-t(sample_id)[2:ncol(sample_id)]
  
  xCell_tissue<-xCell[sample_id,]

  apply(xCell_tissue,2,cor,y=as.numeric(Edges_loss),method="spearman")->xCell_loss_cor[,i]
  apply(xCell_tissue,2,function(x){return(cor.test(x,Edges_loss,method="spearman")$p.value)})->xCell_loss_p[,i]

}
xCell_loss_p[is.na(xCell_loss_p)]<-1
xCell_loss_cor[is.na(xCell_loss_cor)]<-0


heatmap <- melt(xCell_loss_cor) %>% rename(replace=c("X1"="Cell","X2"="Tissue",
                                          "value"="cor")) %>%
  mutate(pvalue=melt(xCell_loss_p)[,3]) %>%
  mutate(signif = sapply(xCell_loss_p, function(x) myfun(x)))

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
