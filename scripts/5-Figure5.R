
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggridges)
library(cowplot)
library(AUCell)
library(gridExtra)
library(Biobase)
library(genefilter)
library(ggpubr)
library(Hmisc)
library(GSEABase)

dir.create("~/CDC3/results/figure5")

#************************************
# Figure 5A
#************************************

load("~/CDC3/data/kidney_datasets/all_kidney_histologies_expression_dataset.RData")
unscaled<-unscaled[rowSds(exprs(unscaled))>0,]
unscaled.s<-unscaled
exprs(unscaled.s)<-t(scale(t(exprs(unscaled)),scale=T))

load("~/CDC3/data/kidney_datasets/GSE131685/GSE131685_processed.RData")
kid.subset<-subset(kid, idents = c(1,2,3,6,7,8,10))
kid.subset<-FindVariableFeatures(kid.subset, selection.method = "vst", nfeatures = 5000)
varfeat<-VariableFeatures(kid.subset)
exp.kid<-as.matrix(GetAssayData(kid.subset,slot="scale.data")[varfeat,])
exp.kid.s<-t(scale(t(exp.kid),scale=T))

cg<-intersect(rownames(unscaled),rownames(exp.kid))
unscaled.s<-unscaled.s[cg,]
exp.kid.s<-exp.kid.s[cg,]
identical(rownames(unscaled.s),rownames(exp.kid.s))

correl.s<-cor(exprs(unscaled.s),exp.kid.s,method="pearson")
correl.s.mean<-t(apply(correl.s,1,function(x)tapply(x,Idents(kid.subset),median)))
mat<-t(correl.s.mean)

histology<-unscaled.s$histology
names(histology)<-colnames(unscaled.s)
histology<-sort(histology)
mat<-mat[,names(histology)]

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
df<-data.frame(Spearman=correl.s.mean[,1],cluster=colnames(correl.s.mean)[1],histology=unscaled.s$histology,stringsAsFactors = F)
for(i in 1:ncol(correl.s.mean)){
  temp<-data.frame(Spearman=correl.s.mean[,i],cluster=colnames(correl.s.mean)[i],histology=unscaled.s$histology,stringsAsFactors = F)
  df<-rbind(df,temp)
}
df$cluster[df$cluster==1]<-"Proximal convoluted tubule cells"
df$cluster[df$cluster=="2"]<-"Proximal tubule cells"
df$cluster[df$cluster=="3"]<-"Proximal straight tubule cells"
df$cluster[df$cluster=="6"]<-"Glomerular parietal epithelial cells"
df$cluster[df$cluster=="7"]<-"Distal tubule cells"
df$cluster[df$cluster=="8"]<-"Collecting duct principal cells"
df$cluster[df$cluster=="10"]<-"Collecting duct intercalated cells"
df$cluster<-factor(df$cluster,levels=c("Proximal convoluted tubule cells","Proximal tubule cells","Proximal straight tubule cells",
                                       "Glomerular parietal epithelial cells","Distal tubule cells","Collecting duct intercalated cells",
                                       "Collecting duct principal cells"
))
df<-df[df$histology!="normal",]
df$histology[df$histology!="ccRCC"]<-capitalize(df$histology[df$histology!="ccRCC"])
df$histology<-factor(df$histology,levels=c("CDC","ccRCC","Papillary","Chromophobe","Oncocytoma","RMC","UTUC"))

tiff("~/CDC3/results/figure5/figure5A.tiff",width=8,height=5,units = "in",res=600,compression="lzw")
ggplot(data=df, aes(x=cluster,y=Spearman, fill=cluster)) +
  geom_boxplot(coef=NULL) +
  scale_fill_manual(values=cbPalette) +
  theme_pubr(border=T,x.text.angle=45,legend="top") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size=12,face="bold"),
        legend.text = element_text(size=9),
        legend.title = element_blank()) +
  labs(fill="Cell type") +
  ylab("Spearman's rho") +
  facet_wrap(~histology,nrow=2) +
  guides(fill=guide_legend(ncol=3))
dev.off()

#************************************
# Figure 5B
#************************************

gs<-scan("~/CDC3/results/figure3/INT_CDC_signature_gene_list.txt",what="character")
setlist<-GeneSet(setName="INT-CDC",geneIds=gs)

exp<-GetAssayData(kid,slot="scale.data")
cells_rankings <- AUCell_buildRankings(exp,nCores = 12)
cells_AUC <- AUCell_calcAUC(setlist, cells_rankings, nCores = 12)
kid@meta.data$CDC_AUCell<-cells_AUC@assays@data@listData$AUC[1,]

kid.subset<-subset(kid, idents = c(1,2,3,6,7,8,10))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
df<-data.frame(CDC_score=kid.subset$CDC_AUCell,cluster=Idents(kid.subset))
df$cluster<-as.vector(df$cluster)
df$cluster[df$cluster==1]<-"Proximal convoluted tubule cells"
df$cluster[df$cluster=="2"]<-"Proximal tubule cells"
df$cluster[df$cluster=="3"]<-"Proximal straight tubule cells"
df$cluster[df$cluster=="6"]<-"Glomerular parietal epithelial cells"
df$cluster[df$cluster=="7"]<-"Distal tubule cells"
df$cluster[df$cluster=="8"]<-"Collecting duct principal cells"
df$cluster[df$cluster=="10"]<-"Collecting duct intercalated cells"
df$cluster<-factor(df$cluster,levels=c("Proximal convoluted tubule cells","Proximal tubule cells","Proximal straight tubule cells",
                                       "Glomerular parietal epithelial cells","Distal tubule cells","Collecting duct intercalated cells",
                                       "Collecting duct principal cells"
))
p1<-format(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Collecting duct intercalated cells"])$p.value,digits=3,scientific=T)
p2<-format(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Distal tubule cells"])$p.value,digits=3,scientific=T)
p3<-format(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Glomerular parietal epithelial cells"])$p.value,digits=3,scientific=T)
p4<-format(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Proximal straight tubule cells"])$p.value,digits=3,scientific=T)
p5<-format(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Proximal tubule cells"])$p.value,digits=3,scientific=T)
p6<-format(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Proximal convoluted tubule cells"])$p.value,digits=3,scientific=T)

plot1<-ggplot(df, aes(x=CDC_score,y=cluster,fill=cluster)) + 
  geom_density_ridges(quantile_lines=T, quantiles=2,scale=1.5) +
  scale_fill_manual(values = cbPalette) +
  labs(x="AUCell score",y="") +
  theme_pubr(legend="none") +
  theme(axis.text.y = element_text(size=10.5),axis.text.x=element_text(size=10),
        axis.title.x = element_text(size=10),
        plot.margin = margin(t=0.05, r=0, b=0.05, l=0,unit = "inches")) +
  geom_text(x=0.23,y=1.22,label=paste("p =",p6),fontface="italic",size=3) +
  geom_text(x=0.23,y=2.22,label=paste("p =",p5),fontface="italic",size=3) +
  geom_text(x=0.23,y=3.22,label=paste("p =",p4),fontface="italic",size=3) +
  geom_text(x=0.23,y=4.22,label=paste("p =",p3),fontface="italic",size=3) +
  geom_text(x=0.23,y=5.22,label=paste("p =",p2),fontface="italic",size=3) +
  geom_text(x=0.23,y=6.22,label=paste("p =",p1),fontface="italic",size=3)

# random permutations

gs<-gs[gs%in%rownames(exp)]
randomList<-vector("list",length=1000)
names(randomList)<-paste0("randomGS",1:1000)
set.seed(111)
for(i in 1:length(randomList)){
  randomList[[i]]<-sample(rownames(exp),size=27,replace = F)
}
for(i in 1:length(randomList)){
  randomList[[i]]<-GeneSet(setName=names(randomList)[[i]],geneIds=randomList[[i]])
}
gscolc<-GeneSetCollection(randomList)
cells_AUC <- AUCell_calcAUC(gscolc, cells_rankings, nCores = 4)
AUCscore<-cells_AUC@assays@data@listData$AUC
mean.AUCscore<-t(apply(AUCscore,1,function(x)tapply(x,Idents(kid),median)))
mean.AUCscore<-mean.AUCscore[,-c(4,5,9)]

df<-data.frame(CDC_score=as.vector(mean.AUCscore),cluster=rep(colnames(mean.AUCscore),each=nrow(mean.AUCscore)),stringsAsFactors = F)
df$cluster[df$cluster==1]<-"Proximal convoluted tubule cells"
df$cluster[df$cluster=="2"]<-"Proximal tubule cells"
df$cluster[df$cluster=="3"]<-"Proximal straight tubule cells"
df$cluster[df$cluster=="6"]<-"Glomerular parietal epithelial cells"
df$cluster[df$cluster=="7"]<-"Distal tubule cells"
df$cluster[df$cluster=="8"]<-"Collecting duct principal cells"
df$cluster[df$cluster=="10"]<-"Collecting duct intercalated cells"
df$cluster<-factor(df$cluster,levels=c("Proximal convoluted tubule cells","Proximal tubule cells","Proximal straight tubule cells",
                                       "Glomerular parietal epithelial cells","Distal tubule cells","Collecting duct intercalated cells",
                                       "Collecting duct principal cells"
))
p1<-round(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Collecting duct intercalated cells"])$p.value,3)
p2<-round(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Distal tubule cells"])$p.value,3)
p3<-round(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Glomerular parietal epithelial cells"])$p.value,3)
p4<-round(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Proximal straight tubule cells"])$p.value,3)
p5<-round(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Proximal tubule cells"])$p.value,3)
p6<-round(wilcox.test(df$CDC_score[df$cluster=="Collecting duct principal cells"],df$CDC_score[df$cluster=="Proximal convoluted tubule cells"])$p.value,3)

plot2<-ggplot(df, aes(x=CDC_score,y=cluster,fill=cluster)) + 
  geom_density_ridges(quantile_lines=T, quantiles=2,scale=0.5) +
  scale_fill_manual(values = cbPalette) +
  labs(x="AUCell score",y="") +
  theme_pubr(legend="none") +
  theme(axis.text.y = element_blank(),axis.text.x=element_text(size=10),
        axis.title.x = element_text(size=10),
        plot.margin = margin(t=0.05, r=2.5, b=0.05, l=0,unit="inches")) +
  geom_text(x=0.065,y=1.22,label=paste("p =",p6),fontface="italic",size=3) +
  geom_text(x=0.065,y=2.22,label=paste("p =",p5),fontface="italic",size=3) +
  geom_text(x=0.065,y=3.22,label=paste("p =",p4),fontface="italic",size=3) +
  geom_text(x=0.065,y=4.22,label=paste("p =",p3),fontface="italic",size=3) +
  geom_text(x=0.065,y=5.22,label=paste("p =",p2),fontface="italic",size=3) +
  geom_text(x=0.065,y=6.22,label=paste("p =",p1),fontface="italic",size=3)

tiff("~/CDC3/results/figure5/figure5B.tiff",width=10,height=4.5,units = "in",res=600,compression="lzw")
plot_grid(plot1,plot2)
dev.off()
