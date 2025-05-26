library(ggunchained)
library(ggforce)
library(data.table)
library(reshape2)
load("/home/owht/Data/KIZ/data/pCPM2/Boxplot_35_65_75.RData")
ggplot(Boxplot_35_65_75,
       aes(fill=ages>40,
           y=log(Edges_gain+1),
           x=tissue))+
  geom_boxplot(draw_quantiles = c(0.5))+theme_classic()
Boxplot_35_65_75$loss_ratio<-Boxplot_35_65_75$Edges_loss/
  (Boxplot_35_65_75$Edges_gain+Boxplot_35_65_75$Edges_loss)

Boxplot_35_65_75<-as.data.table(Boxplot_35_65_75)
Boxplot_35_65_75$all<-Boxplot_35_65_75$Edges_gain + Boxplot_35_65_75$Edges_loss
res<-data.frame(tissue = unique(Boxplot_35_65_75$tissue),
                p=0)
for(i in 1:nrow(res)){
  fit<-wilcox.test(Boxplot_35_65_75[tissue == res[i,]$tissue & ages>40,]$all,
                   Boxplot_35_65_75[tissue == res[i,]$tissue & ages<40,]$all)
  res[i,]$p<-fit$p.value
}
res$padj<-p.adjust(res$p,method = "BH")
