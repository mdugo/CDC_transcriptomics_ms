

library(ggplot2)
library(ggsci)
library(reshape2)
library(genefilter)
library(singscore)
library(GSEABase)
library(Hmisc)
library(limma)
library(ggpubr)

dir.create("~/CDC3/results/figure3",recursive=T)

#*************************
# Figure 3A
#*************************

load("~/CDC3/data/kidney_datasets/int/INT_processed.RData")
patient<-factor(dataset$Patient_ID)
f<-factor(dataset$histology)
design <- model.matrix(~ 0+f)
colnames(design)<-levels(f)
corfit <- duplicateCorrelation(dataset,design,block=patient)
corfit$consensus
fit <- lmFit(dataset,design,block=patient,correlation=corfit$consensus)
contrast.matrix <- makeContrasts(CDC-Normal,
                                 CDC-ccRCC,
                                 ccRCC-Normal,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res<-topTable(fit2, number=nrow(dataset), adjust="BH")
res<-res[,-c(2,3)]
x<-decideTests(fit2,method="separate",adjust.method="BH",p.value=0.25,lfc=1)@.Data
x<-x[rownames(res),] 
identical(rownames(res),rownames(x)) 
res<-cbind(res, x)
colnames(res)[4:6]<-c("logFC_CDC_vs_Normal","logFC_CDC_vs_ccRCC","logFC_ccRCC_vs_Normal")
colnames(res)<-gsub(" - ","_vs_",colnames(res))

mygenes<-res$SYMBOL[res$logFC_CDC_vs_Normal>=1 & res$logFC_CDC_vs_ccRCC>=1 & res$logFC_ccRCC_vs_Normal<log2(1.5) & res$adj.P.Val<0.25]
write.table(mygenes,file="~/CDC3/results/figure3/INT_CDC_signature_gene_list.txt",sep="\t",quote=F,col.names=F,row.names=F)
subset<-dataset[fData(dataset)$SYMBOL%in%mygenes,]
rownames(subset)<-fData(subset)$SYMBOL
medie<-sort(rowMeans(exprs(subset)[,subset$histology=="CDC"]),decreasing=T)
subset<-subset[names(medie),]

x<-reshape2::melt(t(exprs(subset)))
pdata<-data.frame(Var1=colnames(subset),Histology=subset$histology,stringsAsFactors = F)
x<-merge(x,pdata,by="Var1",all.x=T)
x<-x[,-1]
x$Histology<-factor(x$Histology,levels=c("CDC","ccRCC","Normal"))

g1<-ggplot(data=x) +
  geom_boxplot(aes(x=Var2, y=value, fill=Histology),coef=NULL) + 
  theme_pubr() +
  labs(x="",y="log2(expression)") +
  scale_fill_manual(values=c("#BC3C29FF","#0072B5FF","#FFDC91FF")) +
  theme(text=element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1,color="black"),
        axis.text.y=element_text(color="black"),
        legend.position = "top",
        plot.margin=margin(t=0,r=0.1,b=1.5,l=0.1,unit="cm"))

#*************************
# Figure 3B
#*************************

gs<-scan("~/CDC3/results/figure3/INT_CDC_signature_gene_list.txt",what="character")

load("~/CDC3/data/kidney_datasets/all_kidney_histologies_expression_dataset.RData")
unscaled<-unscaled[rowSds(exprs(unscaled))>0,unscaled$histology!="normal"]

rankData<-rankGenes(unscaled)
score<-simpleScore(rankData,upSet=gs)

score<-data.frame(sample=colnames(unscaled),Histology=unscaled$histology,score=score$TotalScore,dataset=unscaled$dataset,stringsAsFactors = F)
score$dataset2<-score$dataset
score$dataset2[score$dataset2=="GSE11151" & score$Histology=="CDC"]<-"GSE11151-CDC"
score<-score[-which(score$dataset2=="INT" & score$Histology=="ccRCC"),]
score$dataset2[score$dataset2=="UTUC-CORNELL-BAYLOR-MDACC"]<-"CBM-UTUC"
score$dataset2[score$dataset2=="RMC-PRJNA605003"]<-"PRJNA605003-RMC"
score$dataset2[score$dataset2=="ICGC"]<-"ICGC-RECA-EU"

score$dataset2<-factor(score$dataset2,levels=levels(factor(score$dataset2))[c(25,23,30,3,26,9,1,27,28,29,24,2,4:8,10:22)])
score$Histology[score$Histology!="ccRCC"]<-capitalize(score$Histology[score$Histology!="ccRCC"])
score$Histology<-factor(score$Histology,levels=c("CDC","ccRCC","Papillary","Chromophobe","Oncocytoma","RMC","UTUC"))


g2<-ggplot(data=score, aes(x=dataset2,y=score,color=Histology)) + 
  geom_jitter(position=position_jitter(0.15),cex=0.8) +
  theme_pubr() +
  theme(text=element_text(size=12),
        axis.text.x = element_text(angle=45,hjust=1),
        plot.margin=margin(t=0,r=0.1,b=0,l=0.1,unit="cm")) +
  scale_color_manual(values=c("#BC3C29","#0072B5","#E18727","#20854E","#EE4C97","#7876B1","#7E6148")) +
  stat_summary(fun.data=mean_sd,geom="pointrange",color="black",size=0.3) + 
  labs(x="", y="CDC signature score") +
  guides(color=guide_legend(override.aes = list(size=2.5)))

plotlist<-list(g1, g2)
tiff("~/CDC3/results/figure3/figure3.tiff",width=9,height=9,units = "in", res = 600, compression = "lzw")
cowplot::plot_grid(plotlist=plotlist,ncol=1, labels=c("A","B"))
dev.off()
