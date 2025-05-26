setwd("/home/owht/Data/KIZ/data/pCPM/Result/Summary/GSVA")
library(ggplot2)
library(data.table)

#PLOT TNF-NFKB vs Heart
load("./Hearthsig.es.RData")
load("../../Result/Heart/Result.RData")
tmp<-data.frame(NFKB_TNF=Hsig_es[1,],Loss=Edges_loss)
tmp$sample<-rownames(tmp)
samples<-fread("/home/owht/Data/KIZ/data/pCPM/Dataset/GTEx/GTEx_2017-06-05_v8_sampleinfo.txt")
samples<-as.data.frame(samples)
rownames(samples)<-samples$SAMPID
tmp<-cbind(tmp,samples[rownames(tmp),])

ggplot(tmp,aes(y=log(Edges_loss),x=NFKB_TNF))+
  geom_point(aes(color=))+
  geom_smooth(method = "lm")+theme_classic()

pcs<-read.table("~/Data/KIZ/data/pCPM2/PC_top5/heart.txt")
tmp<-cbind(tmp[rownames(pcs),],pcs)
fit<-lm(data = tmp,
        Loss~NFKB_TNF+SMTSD)

cor.test(Hsig_es[1,],log(Edges_loss),method = "spearman")
#0.6428 p<2.2e-16

#PLOT BILE vs Pituitary
load("./Pituitaryhsig.es.RData")
load("../../Result/Pituitary/Result.RData")
tmp<-data.frame(NFKB_TNF=Hsig_es[44,],Loss=Edges_loss)
pcs<-read.table("~/Data/KIZ/data/pCPM2/PC_top5/heart.txt")
tmp<-cbind(tmp[rownames(pcs),],pcs)
fit<-lm(data = tmp,
        Loss~NFKB_TNF+PC1+PC2+PC3+PC4+PC5)
ggplot(tmp,aes(y=log(Edges_loss),x=NFKB_TNF))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[44,],log(Edges_loss),method = "spearman")
#-0.4227 p=2.639e-13

#PLOT ROS vs OVA
load("./Ovaryhsig.es.RData")
load("../../Result/Ovary/Result.RData")
tmp<-data.frame(ROS=Hsig_es[36,],Loss=Edges_loss)
pcs<-read.table("~/Data/KIZ/data/pCPM2/PC_top5/ovary.txt")
tmp<-cbind(tmp[rownames(pcs),],pcs)
fit<-lm(data = tmp,
        Loss~ROS+PC1+PC2+PC3+PC4+PC5)

ggplot(tmp,aes(y=log(Edges_loss),x=ROS))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[36,],Edges_loss,method = "spearman")
#0.586 p=5.845e-16

#PLOT PEROXIOME vs BLOOD
load("./Bloodhsig.es.RData")
load("../../Result/Blood/Result.RData")
tmp<-data.frame(ROS=Hsig_es[45,],Loss=Edges_loss)
pcs<-read.table("~/Data/KIZ/data/pCPM2/PC_top5/blood.txt")
tmp<-cbind(tmp[rownames(pcs),],pcs)
fit<-lm(data = tmp,
        Loss~ROS+PC1+PC2+PC3+PC4+PC5)
ggplot(tmp,aes(y=log(Edges_loss),x=ROS))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[45,],Edges_loss,method = "spearman")$p.value
#-0.62523 p<2.2e-16

#PLOT OXI vs LUNG
load("./Lunghsig.es.RData")
load("../../Result/Lung/Result.RData")
tmp<-data.frame(ROS=Hsig_es[34,],Loss=Edges_loss)
pcs<-read.table("~/Data/KIZ/data/pCPM2/PC_top5/lung.txt")
tmp<-cbind(tmp[rownames(pcs),],pcs)
fit<-lm(data = tmp,
        Loss~ROS+PC1+PC2+PC3+PC4+PC5)

ggplot(tmp,aes(y=log(Edges_loss),x=ROS))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()
cor.test(Hsig_es[34,],Edges_loss,method = "spearman")
#-0.599 p<2.2e-16

####PC sample
heart<-read.table("/home/owht/Data/KIZ/data/pCPM2/PC_top5/heart.txt")
sampleinfo<-read.csv("/home/owht/Data/KIZ/data/pCPM/Result/Summary/Stemness/aging-16-205717-s004.csv")
sample_heart<-sampleinfo[sampleinfo$Sample.ID %in% rownames(heart),]
sample_heart<-cbind(sample_heart,heart[sample_heart$Sample.ID,])
ggplot(sample_heart,aes(x=PC1,y=PC2,color=Tissue.Site.Detail))+
  geom_point()+theme_classic()

tmp<-cbind(sample_heart,tmp[sample_heart$Sample.ID,])
ggplot(tmp,aes(y=log(Edges_loss),x=NFKB_TNF))+
  geom_point(aes(color=Age))+
  geom_smooth(method = "lm")+theme_classic()+
  facet_wrap(~SMTSD)

ggplot(tmp,aes(y=log(Edges_loss),x=NFKB_TNF))+
  geom_point(aes(color=Age))+
  geom_smooth(method = "lm")+theme_classic()+
  facet_wrap(~SMTSD)