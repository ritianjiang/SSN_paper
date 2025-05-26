load("/home/owht/Data/KIZ/data/pCPM/Result/Summary/Loss_sig_metrics/All_sample_edge_change.RData")
All_sample_edge_change<-All_sample_edge_change[All_sample_edge_change$type == "loss",]
indinfo<-read.table("/home/owht/Data/KIZ/data/pCPM/Dataset/GTEx/GTEx_2017-06-05_v8_indinfo.txt",sep="\t",header=T)

library(stringr)
library(ggunchained)
str_split(All_sample_edge_change$ID,pattern = "-",simplify = T)[,1:2]->tmp

All_sample_edge_change$ind<-paste0(tmp[,1],"-",tmp[,2])
All_sample_edge_change<-merge(All_sample_edge_change,indinfo,by.x=5,by.y=1)
All_sample_edge_change$disease = ifelse(All_sample_edge_change$DTHHRDY %in% c(3,4),
                                        yes = "Disease",no = "healthy")

tmp<-All_sample_edge_change[All_sample_edge_change$tissue == "Brain",]
tmp<-na.omit(tmp)
ggplot(tmp,aes(x=AGE,y=counts,fill=disease,color=disease))+
  geom_split_violin(draw_quantiles = c(0.5))+theme_bw()

library(ggpubr)

wilcox.test(tmp[tmp$AGE=="60-69" & tmp$disease == "Disease",]$counts,
            tmp[tmp$AGE=="60-69" & tmp$disease != "Disease",]$counts)

Res<-expand.grid(unique(All_sample_edge_change$tissue),
                 unique(All_sample_edge_change$AGE))
colnames(Res)<-c("tissue","age")
Res$p<- -1; Res$bias<-0
Res$Disease_n<-0; Res$healthy_n<-0

ts<-table(All_sample_edge_change$disease,All_sample_edge_change$tissue) %>% t
ts<-ts[ts[,1]>50 & ts[,2]>50,]

Res<-Res[Res$tissue %in% rownames(ts),]
for(i in 1:nrow(Res)){
  tmp<-All_sample_edge_change[All_sample_edge_change$tissue == Res[i,]$tissue &
                                All_sample_edge_change$AGE == Res[i,]$age,]
  Res[i,]$p<-wilcox.test(tmp[tmp$disease == "Disease",]$counts,
                         tmp[tmp$disease != "Disease",]$counts)$p.value
  Res[i,]$bias<-mean(tmp[tmp$disease == "Disease",]$counts) - 
                 mean(tmp[tmp$disease != "Disease",]$counts)
  Res[i,]$Disease_n<-sum(tmp$disease == "Disease")
  Res[i,]$healthy_n<-sum(tmp$disease != "Disease")
}

library(qvalue)
#qvalue(Res$p,lambda=seq(0.1,0.7,0.05))$qvalues->Res$padj
p.adjust(Res$p,method = "BH")->Res$padj


library(reshape2)

Res2<-cast(Res,formula=bias~tissue+age)

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

#colnames(Res)[3]<-"pval"
Res$signif<-sapply(Res$padj,function(x) myfun(x))
Res$age<-factor(Res$age,levels = c("30-39","40-49","50-59","60-69","70-79"))
ggplot(Res,aes(age,tissue,col=bias))+
  geom_tile(color="grey90",fill="white",size=1)+
  geom_point(aes(size = abs(bias)),shape=20) + 
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
