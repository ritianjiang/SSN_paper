#PI4KA PI4K2B CDIPT in small intestine
setwd("/home/wanght/data/CPM_methods/Dosage")
si<-fread("../Exps/Small_Intestine/samples.txt",header=T)
si<-si[si$ID %in% c("PI4KA","PI4K2B","CDIPT"),]

sampleinfo<-read.table("../Dataset/GTEx_2017-06-05_v8_sampleinfo.txt",sep="\t")

load("../SumExp_GTExv8_tpm.filtered_0.3_100_PPI700.RData")
si<-se_GTExv8_tpm.filtered[c("PI4KA","PI4K2B","CDIPT"),
                           colData(se_GTExv8_tpm.filtered)$SMTS == "Small Intestine"]
si_noref<-si[,si$AGE>=30]
si_ref<-si[,si$AGE<=30]

library(corrplot)
corrplot(cor(t(assay(si_ref))),col = colorRampPalette(c("blue","white","red"))(20),
         type="upper")

corrplot(cor(t(assay(si_noref))),col = colorRampPalette(c("blue","white","red"))(20),
         type="upper")
mat<-t(assay(si_noref[,order(si_noref$AGE)]))
library(pheatmap)
pheatmap(mat,show_rownames = F,scale = "column")
