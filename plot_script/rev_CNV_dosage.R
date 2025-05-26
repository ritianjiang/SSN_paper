setwd("~/data/CPM_methods/Dosage/")
load("Enrich_res_rCNV_gwas.RData")
for (i in 1:length(Enrich_res_gwas)){
  Enrich_res_gwas[[i]]<-Enrich_res_gwas[[i]][Enrich_res_gwas[[i]]$p>0,]
  try(Enrich_res_gwas[[i]]$padj<-qvalue(Enrich_res_gwas[[i]]$p-0.0001,lambda=0.25)$qvalue)
}
Enrich_res2<-do.call("rbind",Enrich_res_gwas)
Enrich_res2$p[Enrich_res2$OR>1]<-Enrich_res2$p[Enrich_res2$OR>1]*0.5
Enrich_res2$p[Enrich_res2$OR<1]<- 1-Enrich_res2$p[Enrich_res2$OR<1]*0.5
Enrich_res2$padj<-p.adjust(Enrich_res2$p,method = "BH")
Enrich_res2[Enrich_res2$OR>100,]$OR<-100
ggplot(Enrich_res2,aes(x=Tissue,y=phenotype,color=padj<0.2,size=log(OR)))+
  geom_point()+theme_classic()+
  scale_color_manual(values = c("FALSE"="grey", "TRUE"="red3"))



enrich_adj<-data.frame()
for(i in 1:41){
  tmp<-Enrich_res2[Enrich_res2$phenotype == unique(Enrich_res2$phenotype)[i],]
  tmp$padj<-p.adjust(tmp$p,method="BH")
  enrich_adj<-rbind(enrich_adj,tmp)
}
enrich_adj[enrich_adj$OR>100,]$OR<-100
ggplot(enrich_adj,aes(x=Tissue,y=phenotype,color=padj<0.25,size=log(OR)))+
  geom_point()+theme_classic()+
  scale_color_manual(values = c("FALSE"="grey", "TRUE"="red3"))

enrich_adj2<-data.frame()
for(i in 1:26){
  tmp<-Enrich_res2[Enrich_res2$Tissue == tissues[i],]
  tmp$padj<-p.adjust(tmp$p,method = "BH")
  enrich_adj2<-rbind(enrich_adj2,tmp)
}
enrich_adj2[enrich_adj2$OR>100,]$OR<-100
ggplot(enrich_adj2,aes(x=Tissue,y=phenotype,color=padj<0.1,size=log(OR)))+
  geom_point()+theme_classic()+
  scale_color_manual(values = c("FALSE"="grey", "TRUE"="red3"))

