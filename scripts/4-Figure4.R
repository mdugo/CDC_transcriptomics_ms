

library(ggplot2)
library(ggsci)
library(Biobase)
library(genefilter)
library(singscore)
library(GSEABase)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)

dir.create("~/CDC/results/figure4")

#********************************************************
# correlation between INT-CDC scores and AUC CTRP values
#********************************************************

load("~/CDC/data/ccle_ctrp/CCLE_RNAseq_genes_rpkm_20180929_collSum.RData")
dataset<-dataset[rowSums(exprs(dataset)>0)>=1,]

gs<-scan("~/CDC/results/figure3/INT_CDC_signature_gene_list.txt", what="character")

rankData<-rankGenes(dataset)
score<-simpleScore(rankData,upSet=gs)

ctrpcells<-read.table("~/CDC/data/ccle_ctrp/ctrp_cell_lines.txt",header=T,sep="\t",as.is=T)
ctrpcomp<-read.table("~/CDC/data/ccle_ctrp/ctrp_compounds.txt",header=T,sep="\t",as.is=T,quote="")[,c(1,2,4,5,6,7)]
ctrpdrugs<-read.table("~/CDC/data/ccle_ctrp/ctrp_drugs_auc.txt",header=T,sep="\t",as.is=T)

ctrpcells<-merge(ctrpcells,ctrpdrugs,by="index_ccl",all.x=T)
ctrpcells<-merge(ctrpcells,ctrpcomp,by="index_cpd",all.x=T)

info<-pData(dataset)[,c(1,3,4,5,29,32,33)]
info$cell_line_name<-gsub("_.*","",info$CCLE_ID)
info$CDC_score_up<-score$TotalScore

df<-merge(ctrpcells,info,by="cell_line_name")

cpd<-data.frame(compound=unique(df$compound_name),correl=0,p=0)
for(i in 1:nrow(cpd)){
  temp<-df[df$compound_name==cpd$compound[i],]
  cpd$correl[i]<-cor(temp$CDC_score_up,temp$area_under_sensitivity_curve,method="spearman")
  cpd$p[i]<-cor.test(temp$CDC_score_up,temp$area_under_sensitivity_curve,method="spearman")$p.value
}
cpd$fdr<-p.adjust(cpd$p,method="BH")
cpd<-cpd[order(cpd$p),]
cpd[cpd$correl<0 & cpd$fdr<0.05,]

annot<-unique(ctrpcells[,c(10,13,14)])
annot<-annot[match(cpd$compound,annot$compound_name),]
identical(as.vector(cpd$compound),annot$compound_name)
cpd<-cbind(annot,cpd)
cpd<-cpd[,-4]
write.table(cpd,file="~/CDC/results/figure4/correlation_CTRP_AUC_INT_signature.txt",sep="\t",row.names=F,quote=F)


#********************************************************
# correlation between INT-CDC scores and AUC GDSC values
#********************************************************

drug<-read.table("~/CDC/data/gdsc/TableS4B_cell_lines_drug_auc.txt",header=T,sep="\t",as.is=T,row.names=1)
drug<-t(drug[,-1])
rownames(drug)<-gsub("X","",rownames(drug))

drugannot<-read.table("~/CDC/data/gdsc/Screened_Compounds.txt",sep="\t",header=T,as.is=T)
rownames(drugannot)<-drugannot$DRUG_ID
drugannot<-drugannot[rownames(drug),]
identical(rownames(drug),rownames(drugannot))

load("~/CDC/data/gdsc/E_MTAB_3610_collMaxVar.RData")
rownames(dataset.coll)<-fData(dataset.coll)$Symbol
rankData<-rankGenes(dataset.coll)
score<-simpleScore(rankData,upSet=gs)

info<-pData(dataset.coll)[,c(1,20,26,27)]
info<-unique(info)
CDC_score_up<-tapply(score$TotalScore,dataset.coll$COSMIC_ID,mean)
CDC_score_up<-data.frame(COSMIC_ID=names(CDC_score_up),CDC_score_up=CDC_score_up)
info<-merge(info,CDC_score_up,by="COSMIC_ID")
rownames(info)<-info$COSMIC_ID

cl<-intersect(colnames(drug),rownames(info))
drug<-drug[,cl]
info<-info[cl,]
identical(colnames(drug),rownames(info))

cpd<-data.frame(drugannot,correl=0,p=0)
for(i in 1:nrow(cpd)){
  cpd$correl[i]<-cor(as.vector(as.matrix(drug[rownames(drug)==cpd$DRUG_ID[i],])),info$CDC_score_up,method="spearman",use = "pairwise.complete.obs")
  cpd$p[i]<-cor.test(as.vector(as.matrix(drug[rownames(drug)==cpd$DRUG_ID[i],])),info$CDC_score_up,method="spearman")$p.value
}
cpd$fdr<-p.adjust(cpd$p,method="BH")
cpd<-cpd[order(cpd$p),]
cpd[cpd$correl<0 & cpd$fdr<0.05,]
write.table(cpd,file="~/CDC/results/figure4/correlation_GDSC_AUC_INT_signature.txt",sep="\t",row.names=F,quote=F)


#*******************************************************************
# retrieve targeted pathway from GDSC data and attach to CTRP data
#*******************************************************************

ctrp<-read.table("~/CDC/results/figure4/correlation_CTRP_AUC_INT_signature.txt",sep="\t",header=T,as.is=T,quote="")
ctrp$drug2<-tolower(gsub(" ","",gsub("-","",ctrp$compound)))
gdsc<-read.table("~/CDC/results/figure4/correlation_GDSC_AUC_INT_signature.txt",sep="\t",header=T,as.is=T)
gdsc$drug2<-tolower(gsub(" ","",gsub("-","",gdsc$DRUG_NAME)))
gdsc$synonym2<-tolower(gsub(" ","",gsub("-","",gdsc$SYNONYMS)))

lista<-vector("list",length=nrow(gdsc))
for(i in 1:length(lista)){
  lista[[i]]<-c(gdsc$drug2[i],unlist(strsplit(gdsc$synonym2[i],",")))
}

ctrp$gdsc_name<-""
for(i in 1:nrow(ctrp)){
  index<-sapply(lista,function(z)ctrp$drug2[i]%in%z)
  if(sum(index)>0){
    ctrp$gdsc_name[i]<-paste(unlist(lista[index]),collapse=";")
  }
}

ctrp$TARGET_PATHWAY<-""
for(i in 1:nrow(ctrp)){
  pathway<-unique(gdsc$TARGET_PATHWAY[gdsc$drug2==ctrp$drug2[i]])
  if(length(pathway)>0){
    ctrp$TARGET_PATHWAY[i]<-pathway
  }
}
ctrp$TARGET_PATHWAY[ctrp$TARGET_PATHWAY==""]<-ctrp$target_or_activity_of_compound[ctrp$TARGET_PATHWAY==""]
write.table(ctrp,file="~/CDC/results/figure4/correlation_CTRP_AUC_INT_signature.txt",sep="\t",row.names=F,quote=F)

#*************************
# Figure 4A
#*************************

ctrp<-read.table("~/CDC/results/figure4/correlation_CTRP_AUC_INT_signature.txt",sep="\t",header=T,as.is=T,quote="")
ctrp$sensitivity<-""
ctrp$sensitivity[ctrp$correl<0&ctrp$fdr<0.05]<-"Active compound"
ctrp$sensitivity[ctrp$correl>0&ctrp$fdr<0.05]<-"Inactive compound"
ctrp$sensitivity<-factor(ctrp$sensitivity,levels=c("Active compound","","Inactive compound"))
ctrp$drug2<-tolower(gsub(" ","",gsub("-","",ctrp$compound)))
ctrp$dataset<-"CTRP"
gdsc<-read.table("~/CDC/results/figure4/correlation_GDSC_AUC_INT_signature.txt",sep="\t",header=T,as.is=T)
gdsc$sensitivity<-""
gdsc$sensitivity[gdsc$correl<0&gdsc$fdr<0.05]<-"Active compound"
gdsc$sensitivity[gdsc$correl>0&gdsc$fdr<0.05]<-"Inactive compound"
gdsc$sensitivity<-factor(gdsc$sensitivity,levels=c("Active compound","","Inactive compound"))
gdsc$drug2<-tolower(gsub(" ","",gsub("-","",gdsc$DRUG_NAME)))
gdsc$dataset<-"GDSC"

x1<-ctrp$drug2[ctrp$correl<0 & ctrp$fdr<0.05]
x2<-gdsc$drug2[gdsc$correl<0 & gdsc$fdr<0.05]
common<-intersect(x1,x2)

gdsc[gdsc$drug2%in%common,]
ctrp[ctrp$drug2%in%common,]

ctrp2<-ctrp[,c(1,7,9,4:6,10,11)]
gdsc2<-gdsc[,c(2,10,5:8,9,11)]
colnames(gdsc2)<-colnames(ctrp2)
df<-rbind(ctrp2,gdsc2)
g1<-ggplot(data=df,aes(x=correl,y=-log10(fdr),fill=sensitivity)) +
  geom_point(shape=21,size=3) +
  theme_pubr(border = T) +
  scale_fill_manual(breaks=c("Active compound","Inactive compound"),values=c("#9970AB","grey70","#E08214")) +
  facet_wrap(~dataset,ncol=2) +
  xlab("Spearman correlation") +
  ylab(bquote(~-Log[10]~ "false discovery rate")) +
  labs(fill="") +
  guides(fill=guide_legend(ncol=2)) +
  theme(strip.text.x = element_text(size=12,face="bold"),legend.text = element_text(size=12), legend.title = element_text(size=12))

#*************************
# Figure 4B
#*************************

cmpd<-intersect(ctrp2$drug2,gdsc2$drug2)

ctrp3<-ctrp2[ctrp2$drug2%in%cmpd,c(4,7,2)]
gdsc3<-gdsc2[gdsc2$drug2%in%cmpd,c(4,7,2)]
gdsc3<-aggregate(gdsc3[,1],by=list(drug2=gdsc3$drug2),mean)
gdsc3$sensitivity<-""
for(i in 1:nrow(gdsc3)){
  temp<-gdsc2[gdsc2$drug2==gdsc3$drug2[i],]
  if(length(temp)>1){
    gdsc3$sensitivity[i]<-as.vector(temp$sensitivity[temp$fdr==min(temp$fdr)])
  } else {
    gdsc3$sensitivity[i]<-as.vector(temp$sensitivity)
  }
}
colnames(gdsc3)[2]<-"correl"
ctrp3<-ctrp3[order(ctrp3$drug2),]
gdsc3<-gdsc3[order(gdsc3$drug2),]
df<-merge(ctrp3,gdsc3,by="drug2",suffixes=c("_ctrp","_gdsc"))
df$sensitivity<-""
df$sensitivity[df$sensitivity_ctrp=="Active compound" & df$sensitivity_gdsc=="Active compound"]<-"CTRP & GDSC"
df$show_name<-""
df$show_name[df$sensitivity=="CTRP & GDSC"]<-capitalize(df$drug2[df$sensitivity=="CTRP & GDSC"])
df$show_name<-gsub("Bms","BMS-",gsub("Tgx","TGX-",gsub("Azd","AZD",df$show_name)))

correl<-cor(df$correl_ctrp,df$correl_gdsc,method="spearman")
pval<-cor.test(df$correl_ctrp,df$correl_gdsc,method="spearman")$p.value

g2<-ggplot(data=df,aes(x=correl_ctrp,y=correl_gdsc,fill=sensitivity,label=show_name)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype="dashed", color = "grey60") +
  geom_smooth(inherit.aes = F,aes(x=correl_ctrp,y=correl_gdsc),method="lm",se=F,color="blue",size=0.5) +
  geom_point(shape=21,size=3) +
  theme_pubr(legend="none") +
  theme(legend.text = element_text(size=12)) +
  scale_fill_manual(breaks=c("CTRP & GDSC"),values=c("grey70","forestgreen","turquoise3")) +
  labs(fill="") +
  xlab("Spearman correlation (CTRP)") +
  ylab("Spearman correlation (GDSC)") +
  geom_text_repel(data=df[df$sensitivity=="CTRP & GDSC",],segment.size = 0.5,segment.color = "grey30",size=4,fontface=2,force=20,box.padding = 0.5,seed = 222) +
  geom_text(x=-0.35,y=0.2,label=paste("r =",round(correl,2),"\np =",format(pval,digits=2,scientific=T)),size=5,hjust=0)

plotlist<-list(g1, g2)
tiff("~/CDC/results/figure4/figure4A-B.tiff",width=12,height=5,units = "in", res = 600, compression = "lzw")
cowplot::plot_grid(plotlist=plotlist, ncol=2, nrow=1, rel_widths = c(1,1), labels=c("A","B"))
dev.off()

#*************************
# Figure 4C
#*************************

# calculate differential expression between CDC and normal

load("~/CDC/data/metadataset_CDC_normal_zscore_quantile.RData")
f<-factor(dataset.zscore$histology)
design<-model.matrix(~0+f)
fit<-lmFit(dataset.zscore,design)
contrast.matrix <- makeContrasts(fCDC-fNormal,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,trend=T)
res<-topTable(fit2, number=nrow(dataset.zscore), adjust="BH")

targets.ctrp<-ctrp[ctrp$drug2%in%df$drug2[df$sensitivity=="CTRP & GDSC"],]
targets.ctrp$gene_name_of_protein_target[targets.ctrp$drug2=="dasatinib"]<-paste("ABL1",targets.ctrp$gene_name_of_protein_target[targets.ctrp$drug2=="dasatinib"],sep=";")
targets.gdsc<-gdsc[gdsc$drug2%in%df$drug2[df$sensitivity=="CTRP & GDSC"],]

t1<-unique(unlist(strsplit(targets.ctrp$gene_name_of_protein_target,";")))
t2<-unique(unlist(strsplit(targets.gdsc$TARGET,", ")))
targets<-union(t1,t2)
targets<-targets[1:13]

td<-matrix(0,nrow=length(targets),ncol=length(targets.ctrp$drug2))
rownames(td)<-targets
colnames(td)<-targets.ctrp$compound_name
td<-td[,colnames(td)!="bleomycin"]
for(i in 1:ncol(td)){
  td[,i]<-ifelse(rownames(td)%in%unlist(strsplit(targets.ctrp$gene_name_of_protein_target[targets.ctrp$compound_name==colnames(td)[i]],";")),1,0)
}
td[td==0]<-NA
td[!is.na(td)]<-"0"
fdr<-td

res<-res[res$Symbol%in%rownames(td),]
res<-res[match(rownames(td),res$Symbol),]
res$direction<-"NS"
res$direction[res$adj.P.Val<0.05 & res$logFC>0]<-"UP"
res$direction[res$adj.P.Val<0.05 & res$logFC<0]<-"DOWN"

for(i in 1:nrow(td)){
  td[i,][td[i,]%in%"0"]<-res$logFC[res$Symbol==rownames(td)[i]]
}
for(i in 1:nrow(fdr)){
  fdr[i,][fdr[i,]%in%"0"]<-res$adj.P.Val[res$Symbol==rownames(fdr)[i]]
}
fdr2<-apply(fdr,2,as.numeric)
rownames(fdr2)<-rownames(fdr)
colnames(fdr2)<-capitalize(colnames(fdr))
fdr2[!is.na(fdr2) & fdr2 < 0.05]<-"*"
fdr2[is.na(fdr2) | fdr2!="*"]<-""

td2<-apply(td,2,as.numeric)
rownames(td2)<-rownames(td)
colnames(td2)<-capitalize(colnames(td2))
td2<-td2[order(td2[,1],td2[,2],td2[,3],td2[,4],td2[,5],td2[,6],td2[,7],td2[,8],td2[,9],td2[,10],decreasing = T),order(td2[1,],td2[2,],td2[3,],td2[4,],td2[5,],td2[6,],td2[7,],td2[8,],td2[9,],td2[10,],td2[11,],td2[12,],td2[13,],decreasing=T)]
fdr2<-fdr2[rownames(td2),colnames(td2)]

lgd_list<-list(Legend(labels = c("Down in CDC", "Not significant", "Up in CDC"), title = "Target modulation", 
                      legend_gp = gpar(fill = c("royalblue","grey60","red2")),ncol = 1,title_position = "topcenter",gap = unit(5, "mm")))
classdrug<-ctrp[ctrp$drug2%in%tolower(gsub("-","",colnames(td2))),c(1,9)]
rownames(classdrug)<-capitalize(classdrug$compound_name)
classdrug<-classdrug[colnames(td2),]
classdrug$TARGET_PATHWAY<-factor(classdrug$TARGET_PATHWAY,levels=unique(classdrug$TARGET_PATHWAY))
topha<-HeatmapAnnotation("Drug class"=classdrug$TARGET_PATHWAY,
                         col=list("Drug class"=c("ABL signaling"="#FF7F0E","EGFR signaling"="#E377C2","ERK MAPK signaling"="#BCBD22","IGFR signaling"="#17BECF")),
                         na_col="white",
                         gp = gpar(col = "white"),
                         annotation_height = 0.5,
                         annotation_name_side="right",
                         annotation_legend_param = list(nrow = 2))

h1<-Heatmap(td2,
            name="log2(FC)",
            show_heatmap_legend = T,
            heatmap_legend_param = list(direction = "horizontal"),
            cluster_rows = F,
            cluster_columns = F,
            top_annotation = topha,
            col=c("royalblue","grey60","red2"),
            rect_gp = gpar(col = "black", lwd = 1),
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 11),
            row_names_gp = gpar(fontsize = 10),
            width=unit(7,"cm"),
            height=unit(4.5,"cm"),
            na_col="white",
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(fdr2[i, j], x, y, gp = gpar(fontsize = 12))
            })
tiff("~/CDC/results/figure4/figure4C.tiff",width=5,height=5,units = "in", res = 600, compression = "lzw")
h2<-draw(h1,heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()


#*************************
# Figure 4D
#*************************

valdrugs<-c("dasatinib","afatinib","lapatinib","gefitinib","trametinib","selumetinib","erlotinib","saracatinib","bosutinib","BMS-536924")
df<-read.table("~/CDC/data/clue_INT_CDC_31_genes.gct",sep="\t",skip=2,header=T,as.is=T,quote="",comment.char="")
colnames(df)[10:ncol(df)]<-df[1,10:ncol(df)]
df<-df[-c(1:5),]
for(i in 10:ncol(df)){
  df[,i]<-as.numeric(df[,i])
}
df<-df[order(df$summary),]


cp<-df[df$name%in%valdrugs,]
cp<-cp[-7,]

scores<-as.matrix(cp[,10:ncol(cp)])
colnames(scores)<-capitalize(colnames(scores))
nomi<-capitalize(cp$name)

col_fun = colorRamp2(c(-100, 0, 100), c("#9970AB","grey95","#E08214"))
h2<-Heatmap(scores,
            name="Score",
            col=col_fun,
            row_labels = nomi,
            cluster_rows = F,
            cluster_columns = F,
            column_split = c(rep("Cell lines",9),"Summary"),
            column_title = c("",""),
            rect_gp = gpar(col = "black", lwd = 1),
            column_names_rot = 45,
            column_names_gp = gpar(fontsize = 12),
            row_names_gp = gpar(fontsize = 13),
            width=unit(9.4,"cm"),
            height=unit(6.3,"cm"),
            na_col="white",
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(round(scores[i, j],1), x, y, gp = gpar(fontsize = 9))
            },
            heatmap_legend_param = list(
              legend_direction = "horizontal", 
              legend_width = unit(4, "cm"),
              title_position = "topcenter")
)
tiff("~/CDC/results/figure4/figure4D.tiff",width=14,height=10.5,units = "cm",res=600,compression="lzw")
draw(h2,heatmap_legend_side = "bottom")
dev.off()
