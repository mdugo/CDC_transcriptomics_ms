
library(Biobase)
library(genefilter)
library(singscore)
library(GSEABase)
library(randomcoloR)
library(ComplexHeatmap)
library(limma)
library(msigdbr)

dir.create("~/CDC3/results/figure6")

load("~/CDC3/data/kidney_datasets/all_kidney_histologies_expression_dataset.RData")
unscaled<-unscaled[rowSds(exprs(unscaled))>0,]

gs = as.data.frame(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
setlist<-vector("list",length = length(unique(gs$gs_name)))
names(setlist)<-unique(gs$gs_name)
for(i in 1:length(setlist)){
  setlist[[i]]<-gs$human_gene_symbol[gs$gs_name==names(setlist)[i]]
}
setlist<-sapply(setlist,function(x){x[x%in%rownames(unscaled)]})
setlist<-setlist[sapply(setlist,length)>=15 & sapply(setlist,length)<=500]
for(i in 1:length(setlist)){
  setlist[[i]]<-GeneSet(setName=names(setlist)[i],geneIds=setlist[[i]])
}
setlist<-GeneSetCollection(setlist)

rankData<-rankGenes(exprs(unscaled))
score<-multiScore(rankData,upSetColc=setlist,knownDirection = F)
score<-score$Scores
save(score,file="~/CDC3/results/figure6/c2cp_reactome_singscore.RData")

## differential expression

load("~/CDC3/results/figure6/c2cp_reactome_singscore.RData")

histo<-factor(unscaled$histology,levels=c("CDC","ccRCC","chromophobe","normal","oncocytoma","papillary","RMC","UTUC"))
cc<-apply(combn(levels(histo),m=2),2,paste,collapse="-")
f<-factor(unscaled$histology)
design<-model.matrix(~0+f)
colnames(design)<-levels(f)
fit<-lmFit(score,design)
contrast.matrix <- makeContrasts(contrasts=cc,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res<-topTable(fit2, number=nrow(score), adjust="BH")

dt<-decideTests(fit2,method="separate",p.value = 0.05)@.Data
dt<-dt[rownames(res),]
identical(rownames(res),rownames(dt))
res<-cbind(res,dt)
colnames(res)[1:28]<-paste0("logFC_",gsub("\\.","_",colnames(res)[1:28]))
res<-data.frame(GeneSet=rownames(res),res)

up<-rownames(dt)[rowSums(dt[,1:7])==7]
dn<-rownames(dt)[rowSums(dt[,1:7])==-7]
res$Direction_CDC<-""
res$Direction_CDC[rownames(res)%in%up]<-"Up"
res$Direction_CDC[rownames(res)%in%dn]<-"Down"
write.table(res,file="~/CDC3/results/figure6/limma_cdc_vs_other_c2cp_reactome.txt",sep="\t",row.names=F,quote=F)


########### attach reactome hierarchy ##############

res<-read.table("~/CDC3/results/figure6/limma_cdc_vs_other_c2cp_reactome.txt",sep="\t",header=T,as.is=T)
res<-res[res$Direction_CDC!="",]
colnames(res)[1]<-"pathway"

res$pathway2<-tolower(gsub("\\:","",gsub("_"," ",gsub("REACTOME_","",res$pathway))))
res$pathway2<-gsub(" ","",res$pathway2)
res$reactome_hierarchy<-""

df<-read.table("~/CDC3/data/reactome_hierarchy.txt",sep="\t",header=T,as.is=T,quote="")

for(i in 1:nrow(res)){
  temp<-df[grep(paste0(";",res$pathway2[i],"$"),df$Hierarchies_GSEA2),]
  res$reactome_hierarchy[i]<-paste(temp$Hierarchy,collapse=" // ")
}
res$reactome_hierarchy[res$reactome_hierarchy==""]<-"TOPLEVEL"
res<-res[,c(1,62,64)]
write.table(res,file="~/CDC3/results/figure6/annotation_table_c2cp_reactome_cdc_vs_other.txt",sep="\t",row.names=F,quote=F)

# heatmap paper

load("~/CDC3/data/kidney_datasets/all_kidney_histologies_expression_dataset.RData")
load("~/CDC3/results/figure6/c2cp_reactome_singscore.RData")
identical(colnames(score),colnames(unscaled))

x<-read.table("~/CDC3/results/figure6/annotation_table_c2cp_reactome_cdc_vs_other.txt",sep="\t",header=T,as.is=T)
x<-x[order(x$Direction,x$reactome_hierarchy),]
x$Categories<-c("Adherens junctions interactions","Detoxification of Reactive Oxygen Species","Biological oxidations","Biological oxidations",
              "Biological oxidations","Metabolism of amino acids and derivatives","Metabolism of amino acids and derivatives",
              "Signaling by GPCR","Signaling by GPCR","Interleukin-37 signaling", rep("Innate immune system",6),"Intraflagellar transport",
              "Signaling by Rho GTPases","Signaling by Rho GTPases","DNA repair")
x$Categories<-factor(x$Categories,levels=unique(x$Categories))
score<-score[x$pathway,]
identical(rownames(score),x$pathway)
medie<-t(apply(score,1,function(x)tapply(x,unscaled$histology,mean)))
medie<-t(scale(t(medie),scale=T))
rownames(medie)<-gsub("REACTOME_","",rownames(medie))

set.seed(124)
catcol<-distinctColorPalette(k=10)
names(catcol)<-unique(x$Categories)
catcol[table(x$Categories)<2]<-"white"

ra<-rowAnnotation(Category=x$Categories,
                  col=list(Category=catcol),
                  na_col="white",
                  #gp = gpar(col = "white"),
                  annotation_height = 0.5,
                  show_legend=F)

rowtitle<-table(x$Categories)
names(rowtitle)[rowtitle<2]<-""
rowtitle<-names(rowtitle)

hm<-Heatmap(medie,
            cluster_rows = F,
            column_labels = c("ccRCC","CDC","Chromophobe","Normal kidney","Oncocytoma","Papillary","RMC","UTUC"),
            name="Average singscore",
            heatmap_legend_param = list(direction = "horizontal"),
            row_names_gp = gpar(fontsize=8),
            column_names_gp = gpar(fontsize=10),
            column_names_rot=90,
            show_row_names = T,
            row_split=x$Categories,
        row_title=rowtitle,
            row_title_rot = 0,
        row_title_gp = gpar(fontsize=9),
            left_annotation = ra,
            width=unit(1.8,"in"))
tiff("~/CDC3/results/figure6/figure6.tiff",width=11,height=6,units="in",res=600,compression="lzw")
draw(hm, heatmap_legend_side = "bottom")
dev.off()