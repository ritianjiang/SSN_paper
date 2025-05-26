library(ggplot2)
library(stringr)
load("/home/owht/Data/KIZ/data/pCPM2/CIBERSORT_res.RData")
heatmap$padj<-p.adjust(heatmap$pvalue,method = "BH")
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
heatmap$padj<-p.adjust(heatmap$pvalue,method = "BH")
heatmap$signif<-sapply(heatmap$padj, function(x) myfun(x))
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
