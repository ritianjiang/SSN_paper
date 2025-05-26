setwd("~/Data/KIZ/data/pCPM2/")
library(ggplot2)
load("Enrich_adj_pheno_fdr.Rdata")

#HP0001250 Seizure
#HP0002086 Abnormality of the respiratory system

#HP0011442 Abnormal central motor function
#HP0100753 Schizophrenia
ggplot(enrich_adj[enrich_adj$phenotype == "HP:0100753",],
       aes(x=Tissue,y=(OR),fill=padj<0.25))+
  geom_bar(stat = "identity") + 
  coord_flip()+theme_classic()+
  labs(title = "Schizophrenia")+
  scale_fill_manual(values = c("FALSE"="grey", "TRUE"="red3"))
