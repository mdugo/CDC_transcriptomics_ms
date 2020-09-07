
library(Biobase)
library(limma)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(edgeR)
library(RColorBrewer)
library(VennDiagram)
library(gprofiler2)
library(data.table)
library(psych)
library(igraph)
library(randomcoloR)

dir.create("~/CDC3/results/figure2")

#*************************
# Figure 2A
#*************************

degsList<-vector("list",length=4)
names(degsList)<-c("INT","GSE89122","WACH","GSE11151")

# INT - CDC vs NORMAL

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
degsList$INT<-topTable(fit2, coef=1, number=nrow(dataset), adjust="BH")

# GSE89122 - CDC vs NORMAL

load("~/CDC3/data/kidney_datasets/GSE89122/GSE89122_raw_counts.RData")
dataset<-dataset[,dataset$patient != "CDC5"]
dge<-DGEList(counts=exprs(dataset), genes=fData(dataset))
isexpr<-rowSums(exprs(dataset)>10) >= 1
dge<-dge[isexpr,]
dge <- calcNormFactors(dge,method="TMM")
f<-factor(dataset$histology)
patient<-factor(dataset$patient)
design<-model.matrix(~0+f+patient)
v<-voom(dge,design=design,plot=F)
fit <- lmFit(v,design)
contrast.matrix <- makeContrasts(fCDC-fNormal, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
degsList$GSE89122<-topTable(fit2, number=nrow(dataset), adjust="BH")

# WACH - CDC vs NORMAL

load("~/CDC3/data/kidney_datasets/wach/wach_raw_counts.RData")
dge<-DGEList(counts=exprs(dataset), genes=fData(dataset))
isexpr<-rowSums(exprs(dataset)>10) >= 1
dge<-dge[isexpr,]
dge <- calcNormFactors(dge,method="TMM")
f<-factor(dataset$histology)
design<-model.matrix(~0+f)
v<-voom(dge,design=design,plot=F)
fit <- lmFit(v,design)
contrast.matrix <- makeContrasts(fCDC-fNormal,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
degsList$WACH<-topTable(fit2, number=nrow(dataset), adjust="BH")

# GSE11151 - CDC vs NORMAL

load("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151/GSE11151_series_matrix_collMaxMean.RData")
dataset.coll<-dataset.coll[,dataset.coll$tumor_type%in%c("adult normal kidney","collecting duct carcinoma")]

f<-factor(ifelse(dataset.coll$tumor_type=="collecting duct carcinoma","CDC","Normal"))
design <- model.matrix(~ 0+f)
fit <- lmFit(dataset.coll,design)
contrast.matrix <- makeContrasts(fCDC-fNormal,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
degsList$GSE11151<-topTable(fit2, number=nrow(dataset.coll), adjust="BH")


up1<-degsList$INT$SYMBOL[degsList$INT$logFC>=1 & degsList$INT$adj.P.Val<0.25]
dn1<-degsList$INT$SYMBOL[degsList$INT$logFC<=-1 & degsList$INT$adj.P.Val<0.25]

up2<-degsList$GSE89122$Symbol[degsList$GSE89122$logFC>=1 & degsList$GSE89122$adj.P.Val<0.05]
dn2<-degsList$GSE89122$Symbol[degsList$GSE89122$logFC<=-1 & degsList$GSE89122$adj.P.Val<0.05]

up3<-degsList$WACH$Symbol[degsList$WACH$logFC>=1 & degsList$WACH$adj.P.Val<0.05]
dn3<-degsList$WACH$Symbol[degsList$WACH$logFC<=-1 & degsList$WACH$adj.P.Val<0.05]

up4<-degsList$GSE11151$SYMBOL[degsList$GSE11151$logFC>=1 & degsList$GSE11151$adj.P.Val<0.05]
dn4<-degsList$GSE11151$SYMBOL[degsList$GSE11151$logFC<=-1 & degsList$GSE11151$adj.P.Val<0.05]

venn.plot.up <- draw.quad.venn(
  area1 = length(up1),
  area2 = length(up2),
  area3 = length(up3),
  area4 = length(up4),
  n12 = length(intersect(up1,up2)),
  n13 = length(intersect(up1,up3)),
  n14 = length(intersect(up1,up4)),
  n23 = length(intersect(up2,up3)),
  n24 = length(intersect(up2,up4)),
  n34 = length(intersect(up3,up4)),
  n123 = length(intersect(intersect(up1,up2),up3)),
  n124 = length(intersect(intersect(up1,up2),up4)),
  n134 = length(intersect(intersect(up1,up3),up4)),
  n234 = length(intersect(intersect(up2,up3),up4)),
  n1234 = length(intersect(intersect(up1,up2),intersect(up3,up4))),
  category = c("INT", "GSE89122", "WACH", "GSE11151"),
  fill = brewer.pal(4,"Set1"),
  cat.cex = 1.1,
  cat.col = brewer.pal(4,"Set1"),
  cat.fontfamily = rep("sans", 4),
  #cat.col = c("black", "black", "black", "black"),
  lwd = rep(1, 4),
  margin=c(0.07),
  fontfamily = rep("sans",15),
  label.col = c("black","black","black","black","firebrick","firebrick","firebrick","black","black","black","firebrick","firebrick","black","black","black"),
  cex=c(2,2,2,2,3,3,3,2,2,2,3,3,2,2,2)/2,
  ind=F
)

venn.plot.dn <- draw.quad.venn(
  area1 = length(dn1),
  area2 = length(dn2),
  area3 = length(dn3),
  area4 = length(dn4),
  n12 = length(intersect(dn1,dn2)),
  n13 = length(intersect(dn1,dn3)),
  n14 = length(intersect(dn1,dn4)),
  n23 = length(intersect(dn2,dn3)),
  n24 = length(intersect(dn2,dn4)),
  n34 = length(intersect(dn3,dn4)),
  n123 = length(intersect(intersect(dn1,dn2),dn3)),
  n124 = length(intersect(intersect(dn1,dn2),dn4)),
  n134 = length(intersect(intersect(dn1,dn3),dn4)),
  n234 = length(intersect(intersect(dn2,dn3),dn4)),
  n1234 = length(intersect(intersect(dn1,dn2),intersect(dn3,dn4))),
  category = c("INT", "GSE89122", "WACH", "GSE11151"),
  fill = brewer.pal(4,"Set1"),
  cat.cex = 1.1,
  cat.col = brewer.pal(4,"Set1"),
  cat.fontfamily = rep("sans", 4),
  #cat.col = c("black", "black", "black", "black"),
  lwd = rep(1, 4),
  margin=c(0.07),
  fontfamily = rep("sans",15),
  label.col = c("black","black","black","black","firebrick","firebrick","firebrick","black","black","black","firebrick","firebrick","black","black","black"),
  cex=c(2,2,2,2,3,3,3,2,2,2,3,3,2,2,2)/2,
  ind=F
)

tiff("~/CDC3/results/figure2/figure2A.tiff",width=4,height=8,res=600,units = "in",compression="lzw")
gridExtra::grid.arrange(grobTree(venn.plot.up),grobTree(venn.plot.dn))
dev.off()

# coun matrix 

up.tot<-unique(c(up1,up2,up3,up4))
matrix.up<-matrix(0,ncol=4,nrow=length(up.tot))
rownames(matrix.up)<-up.tot
colnames(matrix.up)<-c("INT", "WANG", "WACH", "GSE11151")
matrix.up[,1]<-ifelse(rownames(matrix.up)%in%up1,1,0)
matrix.up[,2]<-ifelse(rownames(matrix.up)%in%up2,1,0)
matrix.up[,3]<-ifelse(rownames(matrix.up)%in%up3,1,0)
matrix.up[,4]<-ifelse(rownames(matrix.up)%in%up4,1,0)
matrix.up<-data.frame(Symbol=rownames(matrix.up),matrix.up)
matrix.up$Count<-rowSums(matrix.up[,2:5])
matrix.up$Direction<-1
matrix.up$Direction[matrix.up$Count==4]<-2

dn.tot<-unique(c(dn1,dn2,dn3,dn4))
matrix.dn<-matrix(0,ncol=4,nrow=length(dn.tot))
rownames(matrix.dn)<-dn.tot
colnames(matrix.dn)<-c("INT", "WANG", "WACH", "GSE11151")
matrix.dn[,1]<-ifelse(rownames(matrix.dn)%in%dn1,1,0)
matrix.dn[,2]<-ifelse(rownames(matrix.dn)%in%dn2,1,0)
matrix.dn[,3]<-ifelse(rownames(matrix.dn)%in%dn3,1,0)
matrix.dn[,4]<-ifelse(rownames(matrix.dn)%in%dn4,1,0)
matrix.dn<-data.frame(Symbol=rownames(matrix.dn),matrix.dn)
matrix.dn$Count<-rowSums(matrix.dn[,2:5])
matrix.dn$Direction<--1
matrix.dn$Direction[matrix.dn$Count==4]<--2

matrix.full<-rbind(matrix.up,matrix.dn)
write.table(matrix.full,file="~/CDC3/results/figure2/validation_degs_cdc_vs_normal.txt",sep="\t",row.names=F,quote=F)

#**********************************
# Figure 2B - gProfiler enrichment
#**********************************

up<-matrix.full$Symbol[matrix.full$Count>=3 & matrix.full$Direction>0]
dn<-matrix.full$Symbol[matrix.full$Count>=3 & matrix.full$Direction<0]

gostup<-gost(query=up,evcodes=T,sources = c("GO:BP","REAC"),correction_method = "g_SCS")$result
gostdn<-gost(query=dn,evcodes=T,sources = c("GO:BP","REAC"),correction_method = "g_SCS")$result

fwrite(gostup,file="~/CDC3/results/figure2/gprofiler_validated_genes_up.txt",sep="\t",row.names=F,sep2=c("",";",""))
fwrite(gostdn,file="~/CDC3/results/figure2/gprofiler_validated_genes_dn.txt",sep="\t",row.names=F,sep2=c("",";",""))


#### calculate kappa matrix for gene set overlap ####

# up

intersection<-strsplit(gostup$intersection,",")
names(intersection)<-paste(gostup$term_id,gostup$term_name,sep="//")
flat.up<-matrix(0,nrow=length(intersection),ncol=length(unique(unlist(intersection))))
rownames(flat.up)<-names(intersection)
colnames(flat.up)<-unique(unlist(intersection))
for(i in 1:nrow(flat.up)){
  flat.up[i,]<-ifelse(colnames(flat.up)%in%intersection[[i]],1,0)
}

kappa.up<-matrix(0,ncol=nrow(flat.up),nrow=nrow(flat.up))
rownames(kappa.up)<-rownames(flat.up)
colnames(kappa.up)<-rownames(flat.up)
for(i in 1:nrow(kappa.up)){
  print(paste0(i,"/",nrow(kappa.up)))
  for(j in 1:nrow(kappa.up)){
    kappa.up[i,j]<-cohen.kappa(table(flat.up[i,],flat.up[j,]))$kappa
  }
}

kappa.up[upper.tri(kappa.up)]<-NA
kappa.up<-reshape2::melt(kappa.up, na.rm=T)

# dn

intersection<-strsplit(gostdn$intersection,",")
names(intersection)<-paste(gostdn$term_id,gostdn$term_name,sep="//")
flat.dn<-matrix(0,nrow=length(intersection),ncol=length(unique(unlist(intersection))))
rownames(flat.dn)<-names(intersection)
colnames(flat.dn)<-unique(unlist(intersection))
for(i in 1:nrow(flat.dn)){
  flat.dn[i,]<-ifelse(colnames(flat.dn)%in%intersection[[i]],1,0)
}

kappa.dn<-matrix(0,ncol=nrow(flat.dn),nrow=nrow(flat.dn))
rownames(kappa.dn)<-rownames(flat.dn)
colnames(kappa.dn)<-rownames(flat.dn)
for(i in 1:nrow(kappa.dn)){
  print(paste0(i,"/",nrow(kappa.dn)))
  for(j in 1:nrow(kappa.dn)){
    kappa.dn[i,j]<-cohen.kappa(table(flat.dn[i,],flat.dn[j,]))$kappa
  }
}

kappa.dn[upper.tri(kappa.dn)]<-NA
kappa.dn<-reshape2::melt(kappa.dn, na.rm=T)


# Cytoscape input files

kappa.up<-kappa.up[kappa.up$value>=0.5,]
kappa.up<-kappa.up[kappa.up$Var1!=kappa.up$Var2,] # remove self-loop
kappa.dn<-kappa.dn[kappa.dn$value>=0.5,]
kappa.dn<-kappa.dn[kappa.dn$Var1!=kappa.dn$Var2,] # remove self-loop
kmat<-rbind(kappa.up,kappa.dn)
write.table(kmat,file="~/CDC3/results/figure2/input_edge_list_cytoscape_fig2B.txt",sep="\t",row.names=F,quote=F)
nodeatt<-data.frame(Node=c(union(kappa.up$Var1,kappa.up$Var2),union(kappa.dn$Var1,kappa.dn$Var2)),Direction=c(rep("Up",length(union(kappa.up$Var1,kappa.up$Var2))),rep("Down",length(union(kappa.dn$Var1,kappa.dn$Var2)))))
write.table(nodeatt,file="~/CDC3/results/figure2/input_node_attributes_cytoscape_fig2B.txt",sep="\t",row.names=F,quote=F)

# Network clustering is performed using Cytoscape and Community Cluster (gLay) implemented in ClusterMaker2 App.
# Input node attributes include gene set name and direction of modulation (Up or Down).
# From Cytoscape, run Community Cluster then export the node attributes table with the glayCluster column 
# Clusters are then manually annotated in functional categories. An annotated version of the file is
# provided in folder ~/CDC3/data

tiff("~/CDC3/results/figure2/figure2B.tiff",width=16,height=8,res=600,units = "in",compression="lzw")

net<-read.table("~/CDC3/results/figure2/input_edge_list_cytoscape_fig2B.txt",header=T,sep="\t",as.is=T)
vertattall<-read.csv("~/CDC3/data/input_edge_list_cytoscape_fig2B_annotated.csv",as.is=T)
vertattall<-vertattall[,c(4,2,3,7)]
colnames(vertattall)[2]<-"glayCluster"
vert.up<-vertattall[vertattall$Direction=="Up",]
net.up<-net[net$Var1%in%vert.up$name | net$Var2%in%vert.up$name,]
vert.up<-vert.up[order(vert.up$Category),]

x<-graph_from_data_frame(net.up, directed = F, vertices = vert.up)
identical(names(V(x)),vert.up$name)

markg<-list()
for(i in 1:length(unique(vert.up$Category))){
  markg[[i]]<-which(vert.up$Category%in%unique(vert.up$Category)[i])
}

labels<-unique(data.frame(glayCluster=V(x)$glayCluster,Category=V(x)$Category,stringsAsFactors = F))
labels<-aggregate(labels$glayCluster,by=list(Category=labels$Category),paste,collapse="/")
colnames(labels)[2]<-"glayCluster"
set.seed(113)
labels$clustercol<-distinctColorPalette(length(unique(vert.up$Category)))
rownames(labels)<-1:nrow(labels)

par(mfrow=c(1,2))
set.seed(149)
plot(x,
     xlim=c(-1.2,1.2),
     ylim=c(-2,1),
     layout=layout_with_fr(x),
     vertex.label=NA,
     vertex.size=5,
      mark.groups = markg,
     mark.shape=0.9,
     mark.expand=8,
     mark.border = "black",
     mark.col = labels$clustercol,
     vertex.color=ifelse(V(x)$Direction=="Up","red","blue"),
     edge.width=net$value*2,
     edge.color="grey30",
     axes=F)

legend(x=-0.9,
       y=-1.1,
       legend=paste(labels$glayCluster,labels$Category,sep=" - "),
       fill=labels$clustercol,
       ncol=1,
       bty="n",
       cex=0.8,
       xjust=0,
       yjust=1,
       x.intersp = 0.3,
       y.intersp = 0.9,
       xpd=T)
text(0.30,0.16,labels$glayCluster[1],cex=0.8,font=2,adj=0)
text(0.54,-0.33,labels$glayCluster[2],cex=0.8,font=2,adj=0)
text(-0.94,0.62,labels$glayCluster[3],cex=0.8,font=2,adj=0)
text(0.10,0.87,labels$glayCluster[4],cex=0.8,font=2,adj=0)
text(1.07,0.05,labels$glayCluster[5],cex=0.8,font=2,adj=0)
text(-0.75,0.07,labels$glayCluster[6],cex=0.8,font=2,adj=0)
text(-0.43,1.11,labels$glayCluster[13],cex=0.8,font=2,adj=0)
text(0.89,0.59,labels$glayCluster[7],cex=0.8,font=2,adj=0)
text(-0.97,-0.85,labels$glayCluster[8],cex=0.8,font=2,adj=0)
text(0.97,-0.73,labels$glayCluster[9],cex=0.8,font=2,adj=0)
text(-0.27,-0.87,labels$glayCluster[10],cex=0.8,font=2,adj=0)
text(0.63,0.85,labels$glayCluster[11],cex=0.8,font=2,adj=0)
text(0.91,0.34,labels$glayCluster[12],cex=0.8,font=2,adj=0)

###########

vert.dn<-vertattall[vertattall$Direction=="Down",]
net.dn<-net[net$Var1%in%vert.dn$name | net$Var2%in%vert.dn$name,]
vert.dn<-vert.dn[order(vert.dn$Category),]

x<-graph_from_data_frame(net.dn, directed = F, vertices = vert.dn)
identical(names(V(x)),vert.dn$name)

markg<-list()
for(i in 1:length(unique(vert.dn$Category))){
  markg[[i]]<-which(vert.dn$Category%in%unique(vert.dn$Category)[i])
}

labels<-unique(data.frame(glayCluster=V(x)$glayCluster,Category=V(x)$Category,stringsAsFactors = F))
set.seed(222)
labels$clustercol<-distinctColorPalette(length(unique(vert.dn$Category)))
rownames(labels)<-1:nrow(labels)

set.seed(134)
plot(x,
     xlim=c(-1.2,1.2),
     ylim=c(-2,1),
     layout=layout_with_fr(x),
     vertex.label=NA,
     vertex.size=5,mark.groups = markg,
     mark.shape=0.9,
     mark.expand=8,
     mark.border = "black",
     mark.col = labels$clustercol,
     vertex.color=ifelse(V(x)$Direction=="Up","red","blue"),
     edge.width=net$value*2,
     edge.color="grey30",
     axes=F)

legend(x=-0.9,
       y=-1.1,
       legend=paste(labels$glayCluster,labels$Category,sep=" - "),
       fill=labels$clustercol,
       ncol=1,
       bty="n",
       cex=0.8,
       xjust=0,
       yjust=1,
       x.intersp = 0.3,
       y.intersp = 0.9,
       xpd=T)

text(-0.82,0.16,labels$glayCluster[1],cex=0.8,font=2,adj=0)
text(1.0,0.24,labels$glayCluster[2],cex=0.8,font=2,adj=0)
text(1.07,-0.14,labels$glayCluster[3],cex=0.8,font=2,adj=0)
text(0.25,1.10,labels$glayCluster[4],cex=0.8,font=2,adj=0)
text(-0.12,-0.9,labels$glayCluster[5],cex=0.8,font=2,adj=0)
text(0.78,0.62,labels$glayCluster[6],cex=0.8,font=2,adj=0)
text(-0.3,0.19,labels$glayCluster[7],cex=0.8,font=2,adj=0)
text(-0.15,1.12,labels$glayCluster[8],cex=0.8,font=2,adj=0)
text(0.06,-0.4,labels$glayCluster[9],cex=0.8,font=2,adj=0)
text(0.13,0.38,labels$glayCluster[10],cex=0.8,font=2,adj=0)
text(-0.43,-0.3,labels$glayCluster[11],cex=0.8,font=2,adj=0)
dev.off()