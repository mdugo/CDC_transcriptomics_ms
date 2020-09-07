
library(Biobase)
library(survminer)
library(cowplot)
library(singscore)
library(survival)

dir.create("~/CDC3/results/supp_figure2")

gs<-scan("~/CDC3/results/figure3/INT_CDC_signature_gene_list.txt",what="character")


# KIRC

load("~/CDC3/data/kidney_datasets/tcga_kirc/fpkm_tcga_kirc_collSum.RData")
kirc<-kirc[,kirc$definition=="Primary solid Tumor"]
exprs(kirc)<-log2(exprs(kirc)+1)
clinical<-read.table("~/CDC3/data/TCGA-CDR-SupplementalTableS1.txt",header=T,sep="\t",as.is=T,comment.char="")
cs<-intersect(clinical$bcr_patient_barcode,kirc$patient)
clinical<-clinical[match(cs,clinical$bcr_patient_barcode),]
kirc<-kirc[,match(cs,kirc$patient)]
identical(clinical$bcr_patient_barcode,kirc$patient)

rankData<-rankGenes(kirc)
score<-simpleScore(rankData,upSet=gs)

time<-as.numeric(clinical$OS.time)/365.25
event<-as.numeric(clinical$OS)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
hrtext<-paste("HR for 'High' group = ",round(res$conf.int[1,1],2)," (95% CI, ",round(res$conf.int[1,3],2),"-",round(res$conf.int[1,4],2),")",sep="")
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit1 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR for 'High' group = ",round(res$conf.int[1,1],2)," (95% CI, ",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",format(pval,digits=3,scientific=T),sep="")

g1<-ggsurvplot(fit1, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "KIRC - OS") +
  xlab("Time (years)")
g1$plot<-g1$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0)+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14),title = element_text(size = 16))

time<-as.numeric(clinical$PFI.time)/365.25
event<-as.numeric(clinical$PFI)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit2 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR for 'High' group = ",round(res$conf.int[1,1],2)," (95% CI, ",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",format(pval,digits=3,scientific=T),sep="")

g2<-ggsurvplot(fit2, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "KIRC - PFI")
g2$plot<-g2$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0)+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14),title = element_text(size = 16))

# KIRP

load("~/CDC3/data/kidney_datasets/tcga_kirp/fpkm_tcga_kirp_collSum.RData")
kirp<-kirp[,kirp$definition=="Primary solid Tumor"]
exprs(kirp)<-log2(exprs(kirp)+1)
clinical<-read.table("~/CDC3/data/TCGA-CDR-SupplementalTableS1.txt",header=T,sep="\t",as.is=T,comment.char="")
cs<-intersect(clinical$bcr_patient_barcode,kirp$patient)
clinical<-clinical[match(cs,clinical$bcr_patient_barcode),]
kirp<-kirp[,match(cs,kirp$patient)]
identical(clinical$bcr_patient_barcode,kirp$patient)

rankData<-rankGenes(kirp)
score<-simpleScore(rankData,upSet=gs)

time<-as.numeric(clinical$OS.time)/365.25
event<-as.numeric(clinical$OS)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit1 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR for 'High' group = ",round(res$conf.int[1,1],2)," (95% CI, ",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",round(pval,3),sep="")

g3<-ggsurvplot(fit1, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "KIRP - OS")
g3$plot<-g3$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0)+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14),title = element_text(size = 16))

time<-as.numeric(clinical$PFI.time)/365.25
event<-as.numeric(clinical$PFI)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit2 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR for 'High' group = ",round(res$conf.int[1,1],2)," (95% CI, ",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",round(pval,3),sep="")

g4<-ggsurvplot(fit2, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "KIRP - PFI")
g4$plot<-g4$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0)+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14),title = element_text(size = 16))

# KICH

load("~/CDC3/data/kidney_datasets/tcga_kich/fpkm_tcga_kich_collSum.RData")
kich<-kich[,kich$definition=="Primary solid Tumor"]
exprs(kich)<-log2(exprs(kich)+1)
clinical<-read.table("~/CDC3/data/TCGA-CDR-SupplementalTableS1.txt",header=T,sep="\t",as.is=T,comment.char="")
cs<-intersect(clinical$bcr_patient_barcode,kich$patient)
clinical<-clinical[match(cs,clinical$bcr_patient_barcode),]
kich<-kich[,match(cs,kich$patient)]
identical(clinical$bcr_patient_barcode,kich$patient)

rankData<-rankGenes(kich)
score<-simpleScore(rankData,upSet=gs)

time<-as.numeric(clinical$OS.time)/365.25
event<-as.numeric(clinical$OS)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit1 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR for 'High' group = ",round(res$conf.int[1,1],2)," (95% CI, ",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",round(pval,3),sep="")

g5<-ggsurvplot(fit1, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "KICH - OS")
g5$plot<-g5$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0)+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14),title = element_text(size = 16))

time<-as.numeric(clinical$PFI.time)/365.25
event<-as.numeric(clinical$PFI)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit2 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR for 'High' group = ",round(res$conf.int[1,1],2)," (95% CI, ",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",round(pval,3),sep="")

g6<-ggsurvplot(fit2, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "KICH - PFI")
g6$plot<-g6$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0)+
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14),title = element_text(size = 16))

glist<-list(g1,g3,g5,g2,g4,g6)
res<-arrange_ggsurvplots(glist,ncol=2,nrow=3,print = F)
ggsave("~/CDC3/results/supp_figure2/supplementary_figure2.tiff",plot=res,width=8.5,height = 13,units = "in",dpi=600,compression="lzw")


