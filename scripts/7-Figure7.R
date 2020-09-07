
library(Biobase)
library(edgeR)
library(limma)
library(RColorBrewer)
library(cola)
library(NMF)
library(preprocessCore)
library(ComplexHeatmap)
library(fgsea)
library(data.table)
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(singscore)
library(GSEABase)
library(survminer)
library(ArrayTools)
library(CePa)
library(msigdbr)
library(circlize)

dir.create("~/CDC3/results/figure7")

load("~/CDC3/data/kidney_datasets/metadataset_CDC_zscore_quantile.RData")
exprs(dataset.zscore)<-exprs(dataset.zscore)+abs(min(exprs(dataset.zscore)))+1
colnames(pData(dataset.zscore))[2]<-"Batch"

mat<-exprs(dataset.zscore)

register_NMF()

set.seed(123)
rl = run_all_consensus_partition_methods(
  mat,
  mc.cores = 14,
  max_k = 4,
  scale_rows = F,
  sample_by = "row",
  anno = pData(dataset.zscore)[, c("Batch"), drop = FALSE],
  anno_col = c("INT" = "#E41A1C", "GSE11151" = "#984EA3", "GSE89122" = "#377EB8", "WACH" = "#4DAF4A")
)
save(rl, file = "~/CDC3/results/figure7/cola_subgroups_CDC.RData")
cola_report(rl, output_dir = "~/CDC3/results/figure7/cola_subgroups_CDC_report", mc.cores = 14)

#**************************************
# Figure 7A
#**************************************

load("~/CDC3/results/figure7/cola_subgroups_CDC.RData")
x<-collect_classes(rl, k = 2)
load("~/CDC3/data/kidney_datasets/metadataset_CDC_zscore_quantile.RData")
dataset.zscore$class<-ifelse(get_classes(rl,k=2)$class==1,"C2","C1")

f<-factor(dataset.zscore$class)
design<-model.matrix(~0+f)
fit <- lmFit(dataset.zscore,design)
contrast.matrix <- makeContrasts(fC1-fC2,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,trend=T)
res<-topTable(fit2, number=nrow(dataset.zscore), adjust="BH")
write.table(res,file="~/CDC3/results/figure7/limma_C1_vs_C2.txt",sep="\t",row.names=F,quote=F)

rnk<-res$t
names(rnk)<-res$Symbol
gs = as.data.frame(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
setlist<-vector("list",length = length(unique(gs$gs_name)))
names(setlist)<-unique(gs$gs_name)
for(i in 1:length(setlist)){
  setlist[[i]]<-gs$human_gene_symbol[gs$gs_name==names(setlist)[i]]
}
setlist<-sapply(setlist,function(x){x[x%in%rownames(dataset.zscore)]})

set.seed(1234)
fgseares<-fgsea(pathways = setlist,
                stats = rnk,
                nperm = 10000,
                minSize = 15,
                maxSize = 500)
fgseares<-data.frame(fgseares[order(fgseares$pval),])
fgseares<-fgseares[fgseares$padj<0.05,]
fgseares$Enrichment<-ifelse(fgseares$NES<0,"Negative","Positive")

# attach reactome hierarchy

fgseares$pathway2<-tolower(gsub("\\:","",gsub("_"," ",gsub("REACTOME_","",fgseares$pathway))))
fgseares$pathway2<-gsub(" ","",fgseares$pathway2)
fgseares$reactome_hierarchy<-""
df<-read.table("~/CDC3/data/reactome_hierarchy.txt",sep="\t",header=T,as.is=T,quote="")
for(i in 1:nrow(fgseares)){
  temp<-df[grep(paste0(";",fgseares$pathway2[i],"$"),df$Hierarchies_GSEA2),]
  fgseares$reactome_hierarchy[i]<-paste(temp$Hierarchy,collapse=" // ")
}
fgseares$reactome_hierarchy[fgseares$reactome_hierarchy==""]<-capitalize(tolower(gsub("\\:","",gsub("_"," ",gsub("REACTOME_","",fgseares$pathway[fgseares$reactome_hierarchy==""])))))
fgseares<-fgseares[,-10]
fgseares<-fgseares[order(fgseares$Enrichment,fgseares$reactome_hierarchy),]
write.table(fgseares,file="~/CDC3/results/figure7/fgsea_C1_vs_C2_c2cp_reactome.txt",sep="\t",row.names=F,quote=F)

# figure gene sets

cl<-read.table("~/CDC3/results/figure7/fgsea_C1_vs_C2_c2cp_reactome.txt",sep="\t",header=T,as.is=T)
cl<-cl[cl$padj<0.05,]
freq<-table(cl$Category)
cl<-cl[cl$Category%in%names(freq)[freq>=2],]
cl$Category<-gsub("rna","RNA",gsub("dna","DNA",gsub("gpcr/rtk","GPCR/RTK",cl$Category)))
cl$Category<-capitalize(cl$Category)
cl$Category[cl$Category=="Diseases of signal transduction by growth factor receptors and second messengers"]<-"Diseases of signal transduction"
medie<-sort(tapply(cl$NES,cl$Category,mean),decreasing=F)
medie<-data.frame(Category=names(medie),Media=medie)
cl<-merge(cl,medie,by="Category")
cl$Category<-factor(cl$Category,levels=rownames(medie))

tiff("~/CDC3/results/figure7/figure7A.tiff",width=9.5,height=10,units = "in",res=600,compression="lzw")
ggplot(data=cl, aes(y=Category,x=NES,fill=Enrichment)) +
  geom_jitter(position=position_jitter(width=0.18,height=0.12),cex=2.3,shape=21) +
  theme_pubr(legend="top") +
  scale_fill_manual(values=c("royalblue","red")) +
  xlab("NES") +
  ylab("") +
  theme(axis.text.x = element_text(angle=0,hjust=1,size=15),
        axis.text.y = element_text(angle=0,hjust=1,size=15),
        axis.title.x = element_text(size=15),
        plot.margin=margin(t=0,r=0.1,b=0,l=2,unit="cm"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) +
  geom_hline(yintercept=1:length(cl$Category), linetype="dashed", color = "grey70") +
  guides(fill=guide_legend(ncol=1))
dev.off()

#**************************************
# Figure 7B
#**************************************

# input files for SubMap

dir.create("~/CDC3/results/figure7/submap")

load("~/CDC3/data/kidney_datasets/tcga_kirc/fpkm_tcga_kirc_collSum.RData")
kirc<-kirc[,kirc$definition=="Primary solid Tumor"]
kirc<-kirc[,!is.na(kirc$subtype_mRNA_cluster)]
exprs(kirc)<-log2(exprs(kirc)+1)

output.gct(kirc,filename = "~/CDC3/results/figure7/submap/tcga_kirc")

r1<-paste(ncol(kirc),"4 1")
r2<-c("#",paste(as.vector(kirc$subtype_mRNA_cluster),collapse=" "))
r2<-paste(r2[1],r2[2],collapse=" ")
r3<-paste(as.numeric(kirc$subtype_mRNA_cluster),collapse=" ")
cls<-matrix(c(r1,r2,r3),ncol=1)
write.table(cls,file="~/CDC3/results/figure7/submap/tcga_kirc_subtypes.cls",sep="\t",row.names=F,col.names=F,quote=F)

# TCGA-KIRP

load("~/CDC3/data/kidney_datasets/tcga_kirp/fpkm_tcga_kirp_collSum.RData")
kirp<-kirp[,kirp$definition=="Primary solid Tumor"]
kirp<-kirp[,!is.na(kirp$subtype_mRNA.clusters..3.group.NMF..Rathmell.group.)]
exprs(kirc)<-log2(exprs(kirc)+1)

output.gct(kirp,filename = "~/CDC3/results/figure7/submap/tcga_kirp")

r1<-paste(ncol(kirp),"3 1")
r2<-c("#",paste(as.vector(kirp$subtype_mRNA.clusters..3.group.NMF..Rathmell.group.),collapse=" "))
r2<-paste(r2[1],r2[2],collapse=" ")
r3<-paste(as.numeric(kirp$subtype_mRNA.clusters..3.group.NMF..Rathmell.group.),collapse=" ")
cls<-matrix(c(r1,r2,r3),ncol=1)
write.table(cls,file="~/CDC3/results/figure7/submap/tcga_kirp_subtypes.cls",sep="\t",row.names=F,col.names=F,quote=F)

# CDC

load("~/CDC3/results/figure7/cola_subgroups_CDC.RData")
x<-collect_classes(rl, k = 2)
load("~/CDC3/data/kidney_datasets/metadataset_CDC_zscore_quantile.RData")
dataset.zscore<-dataset.zscore[,dataset.zscore$histology=="CDC"]
dataset.zscore$class<-ifelse(get_classes(rl,k=2)$class==1,"2","1")
output.gct(dataset.zscore,filename = "~/CDC3/results/figure7/submap/CDC_zscore_quantile")

r1<-paste(ncol(dataset.zscore),"2 1")
r2<-c("#",paste(as.vector(dataset.zscore$class),collapse=" "))
r2<-paste(r2[1],r2[2],collapse=" ")
r3<-paste(as.numeric(dataset.zscore$class),collapse=" ")
cls<-matrix(c(r1,r2,r3),ncol=1)
write.table(cls,file="~/CDC3/results/figure7/submap/CDC_subtypes.cls",sep="\t",row.names=F,col.names=F,quote=F)

#### HEATMAPS OF RESULTS

pval<-read.gct("~/CDC3/results/figure7/submap/SubMap_CDC_TCGA_KIRC_Bonferroni_SAmatrix.gct")
colnames(pval)<-c("m1","m2","m3","m4")
rownames(pval)<-c("CDC-S1","CDC-S2")
pval<-round(pval,2)
col_fun = colorRamp2(c(0, 0.5, 1), c("red","orange","cornsilk"))
h1<-Heatmap(pval,
            name = "Bonferroni p-value",
            show_heatmap_legend = F,
            show_row_dend = F,
            show_column_dend = F,
            cluster_rows = F,
            cluster_columns = F,
            height=unit(1.2,"cm"),
            heatmap_legend_param = list(direction = "horizontal"),
            row_names_gp = gpar(fontsize = 10),
            row_names_side = "left",
            column_title = "TCGA-KIRC",
            row_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 0,
            column_names_centered = T,
            col = col_fun,
            rect_gp = gpar(col = "black", lwd = 0.5),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(round(pval[i, j],3), x, y, gp = gpar(fontsize = 9))
            })
tiff("~/CDC3/results/figure7/figure7B.kirc.tiff",height=5,width=5.5,units = "cm",res=600,compression="lzw")
draw(h1, heatmap_legend_side = "bottom")
dev.off()

pval<-read.gct("~/CDC3/results/figure7/submap/SubMap_CDC_TCGA_KIRP_Bonferroni_SAmatrix.gct")
rownames(pval)<-c("m1","m2","m3")
colnames(pval)<-c("CDC-S1","CDC-S2")
col_fun = colorRamp2(c(0, 0.5, 1), c("red","orange","cornsilk"))
h2<-Heatmap(pval,
            show_heatmap_legend = F,
            show_row_dend = F,
            show_column_dend = F,
            cluster_rows = F,
            cluster_columns = F,
            height=unit(1.2,"cm"),
            row_names_gp = gpar(fontsize = 10),
            row_names_side = "left",
            column_title = "TCGA-KIRP",
            row_title_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 0,
            column_names_centered = T,
            col = col_fun,
            rect_gp = gpar(col = "black", lwd = 0.5),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(round(pval[i, j],3), x, y, gp = gpar(fontsize = 9))
            })

tiff("~/CDC3/results/figure7/figure7B.kirp.tiff",height=5,width=5.5,units = "cm",res=600,compression="lzw")
h2
dev.off()


#**************************************
# Figure 7C
#**************************************

res<-read.table("~/CDC3/results/figure7/limma_C1_vs_C2.txt",sep="\t",header=T,as.is=T)
res<-res[order(res$t,decreasing=T),]
up<-res$Symbol[1:150]
res<-res[order(res$t,decreasing=F),]
dn<-res$Symbol[1:150]

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
score<-simpleScore(rankData,upSet=dn)

time<-as.numeric(clinical$OS.time)/365.25
event<-as.numeric(clinical$OS)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit1 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR = ",round(res$conf.int[1,1],2)," (",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",format(pval,digits=3,scientific=T),sep="")

g1<-ggsurvplot(fit1, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "TCGA-KIRC, Genes UP in CDC-S2") +
  xlab("Time (years)") +
  ylab("Survival probability (OS)")
g1$plot<-g1$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0,size=6)+
  theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),axis.title = element_text(size=18),title = element_text(size=20), legend.text = element_text(size = 16), legend.title = element_text(size = 16))

time<-as.numeric(clinical$PFI.time)/365.25
event<-as.numeric(clinical$PFI)
sc<-factor(ifelse(score$TotalScore>=median(score$TotalScore),"high","low"),levels=c("low","high"))
res<-summary(coxph(Surv(time,event)~sc))
cox.zph(coxph(Surv(time,event)~sc))
df<-data.frame(time=time,event=event,class=sc)
fit2 <- survfit(Surv(time, event) ~ class,data=df)
sdf <- survdiff(Surv(time, event) ~ class,data=df)
pval <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
hrtext<-paste("HR = ",round(res$conf.int[1,1],2)," (",round(res$conf.int[1,3],2)," - ",round(res$conf.int[1,4],2),")\n","Log-rank p-value = ",format(pval,digits=3,scientific=T),sep="")

g2<-ggsurvplot(fit2, data=df,pval = F,pval.method = F,legend.labs=c("Low","High")) +
  labs(title = "TCGA-KIRC, Genes UP in CDC-S2") +
  xlab("Time (years)") +
  ylab("Survival probability (PFI)")
g2$plot<-g2$plot+geom_text(x=0, y=0.20, label=hrtext, hjust=0,size=6)+
  theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18),axis.title = element_text(size=18),title = element_text(size=20), legend.text = element_text(size = 16), legend.title = element_text(size = 16))

glist<-list(g1,g2)
res<-arrange_ggsurvplots(glist,ncol=1,nrow=2,print = F)
ggsave("~/CDC3/results/figure7/figure7C.tiff",plot=res,height=9.5,width = 4.75,units = "in",dpi=600,compression="lzw")


