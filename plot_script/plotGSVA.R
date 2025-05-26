setwd("/home/owht/Data/KIZ/data/pCPM/Result/Summary/GSVA")
library(ggplot2)

#PLOT TNF-NFKB vs Heart
load("./Hearthsig.es.RData")
load("../../Result/Heart/Result.RData")
tmp<-data.frame(NFKB_TNF=Hsig_es[1,],Loss=Edges_loss)
ggplot(tmp,aes(y=log(Edges_loss),x=NFKB_TNF))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[1,],log(Edges_loss),method = "spearman")
#0.6428 p<2.2e-16

#PLOT BILE vs Pituitary
load("./Pituitaryhsig.es.RData")
load("../../Result/Pituitary/Result.RData")
tmp<-data.frame(NFKB_TNF=Hsig_es[44,],Loss=Edges_loss)
ggplot(tmp,aes(y=log(Edges_loss),x=NFKB_TNF))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[44,],log(Edges_loss),method = "spearman")
#-0.4227 p=2.639e-13

#PLOT ROS vs OVA
load("./Ovaryhsig.es.RData")
load("../../Result/Ovary/Result.RData")
tmp<-data.frame(ROS=Hsig_es[36,],Loss=Edges_loss)
ggplot(tmp,aes(y=log(Edges_loss),x=ROS))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[36,],Edges_loss,method = "spearman")
#0.586 p=5.845e-16

#PLOT PEROXIOME vs BLOOD
load("./Bloodhsig.es.RData")
load("../../Result/Blood/Result.RData")
tmp<-data.frame(ROS=Hsig_es[45,],Loss=Edges_loss)
ggplot(tmp,aes(y=log(Edges_loss),x=ROS))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[45,],Edges_loss,method = "spearman")$p.value
#-0.62523 p<2.2e-16

#PLOT OXI vs LUNG
load("./Lunghsig.es.RData")
load("../../Result/Lung/Result.RData")
tmp<-data.frame(ROS=Hsig_es[34,],Loss=Edges_loss)
ggplot(tmp,aes(y=log(Edges_loss),x=ROS))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[34,],Edges_loss,method = "spearman")
#-0.599 p<2.2e-16
