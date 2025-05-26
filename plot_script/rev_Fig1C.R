library(magrittr)
library(stringr)
setwd("/home/wanght/data/CPM_methods")
ind_info<-read.table("./Dataset/GTEx_2017-06-05_v8_indinfo.txt",
                     fill = T, # From GTEx Protal
                     header = T,row.names = 1)
ind_info$AGE<-ind_info$AGE %>%
  str_split(,pattern = "-",simplify = T) %>%
  .[,1] %>% as.numeric() %>% `+`(5)        #convert range into int, e.g. 20-29 ==> 25
ind_info[ind_info$SEX == 2,]$SEX<-"F"
ind_info[ind_info$SEX == "1",]$SEX<-"M"
ind_info$SEX<-factor(ind_info$SEX)

Dirs<-list.dirs("Result/")
library(ggplot2)
Summarys<-data.frame(Tissues = str_split(Dirs[2:27],pattern = "//",simplify = T)[,2],
                     gain_tau=0,gain_p=0,
                     loss_tau=0,loss_p=0)
gain_totals<-data.frame()
loss_totals<-data.frame()

for(i in 1:26){
  load(paste0(Dirs[i+1],"/Result.RData"))
  tmp_gain<-data.frame(tissue = Summarys[i,]$Tissues,age=ages,gain=Edges_gain,
                       DTHHRDY = ind_info[names(ages),]$DTHHRDY)
  tmp_loss<-data.frame(tissue = Summarys[i,]$Tissues,age=ages,loss=Edges_loss,
                       DTHHRDY = ind_info[names(ages),]$DTHHRDY)
  gain_totals<-rbind(gain_totals,tmp_gain)
  loss_totals<-rbind(loss_totals,tmp_loss)
  cat(i);cat(" ")
}

ggplot(loss_totals,aes(x=age,y=loss))+geom_point()+geom_smooth(method = "lm")+facet_wrap(~tissue)
ggplot(gain_totals,aes(x=age,y=gain))+geom_point()+geom_smooth(method = "lm")+facet_wrap(~tissue)

ggplot(loss_totals[loss_totals$DTHHRDY!=0,],aes(x=DTHHRDY,y=loss))+geom_point(position = position_jitter())+
  geom_smooth(method = "lm")+facet_wrap(~tissue)

for(i in 1:26){
  tmp_gain<-gain_totals[gain_totals$tissue == Summarys[i,]$Tissues,]
  tmp_loss<-loss_totals[loss_totals$tissue == Summarys[i,]$Tissues,]
  fit_gain<-cor.test(tmp_gain$age,tmp_gain$gain);Summarys[i,]$gain_tau<-fit_gain$estimate;Summarys[i,]$gain_p<-fit_gain$p.value
  fit_loss<-cor.test(tmp_loss$age,tmp_loss$loss);Summarys[i,]$loss_tau<-fit_loss$estimate;Summarys[i,]$loss_p<-fit_loss$p.value
}
Summarys$Tissues<-factor(Summarys$Tissues,
                         level=Summarys$Tissues[order(Summarys$loss_tau)])
Summarys$loss_padj<-p.adjust(Summarys$loss_p,method = "BH")
Summarys$gain_padj<-p.adjust(Summarys$gain_p,method = "BH")

p1<-ggplot(Summarys,aes(x=Tissues,y=gain_tau,fill=gain_padj<0.01))+geom_bar(stat = "identity")+ theme_classic()+
  theme(axis.text.x.bottom = element_text(angle = 45))
p2<-ggplot(Summarys,aes(x=Tissues,y=loss_tau,fill=loss_padj<0.01))+geom_bar(stat = "identity")+theme_classic()+
  theme(axis.text.x.bottom = element_text(size = 0),
        axis.title.x = element_text(size=0))

library(gridExtra)
grid.arrange(p2,p1,nrow=2)
save(loss_totals,gain_totals,file="Summary.RData")
# 1   Adipose_Tissue -0.009770869  0.07164191 2.473614e-02 0.775198462
# 2    Adrenal_Gland  0.018048467  0.16903452 1.483036e-02 0.782236168
# 3            Blood  0.103695336  0.17942442 1.336228e-06 0.026393373
# 4     Blood_Vessel  0.065251799  0.10547287 8.068113e-04 0.064643173
# 5            Brain  0.045972662  0.05528558 9.243535e-03 0.064096039
# 6           Breast  0.064734650  0.13696611 9.243535e-03 0.321480547
# 7            Colon  0.028441458  0.19770664 1.336228e-06 0.619174540
# 8        Esophagus  0.123465824  0.16730757 3.538664e-08 0.000215549
# 9            Heart  0.088198978  0.09698062 9.243535e-03 0.042019456
# 10           Liver  0.037598687 -0.05168763 4.466341e-01 0.718064264
# 11            Lung  0.126705075  0.18222017 1.214311e-04 0.026393373
# 12          Muscle  0.080737846  0.14282544 5.253204e-04 0.070517794
# 13           Nerve  0.020706283  0.09219066 3.624695e-02 0.735198832
# 14           Ovary  0.213956873  0.24455579 5.026932e-03 0.032121218
# 15        Pancreas -0.019669709  0.09063263 1.276763e-01 0.775198462
# 16       Pituitary -0.022685772  0.07191221 2.448609e-01 0.775198462
# 17        Prostate  0.070130072  0.20841358 5.026932e-03 0.497308631
# 18  Salivary_Gland  0.073213909  0.19156142 2.693524e-02 0.575858538
# 19            Skin -0.068132507 -0.03943606 1.260015e-01 0.032121218
# 20 Small_Intestine  0.167201474  0.12671668 1.260015e-01 0.076173071
# 21          Spleen  0.045803951  0.23740350 1.281105e-03 0.650147314
# 22         Stomach  0.107518519  0.15788172 9.243535e-03 0.113246637
# 23          Testis  0.075716468  0.10632904 6.871761e-02 0.321480547
# 24         Thyroid  0.032910651  0.15302246 6.755624e-04 0.604762795
# 25          Uterus  0.242305167  0.28272949 4.845693e-03 0.032121218
# 26          Vagina  0.183673138  0.20692778 2.166346e-02 0.070517794
