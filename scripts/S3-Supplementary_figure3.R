
library(Biobase)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(genefilter)
library(Hmisc)

dir.create("~/CDC3/results/supp_figure3")


load("~/CDC3/data/kidney_datasets/all_kidney_histologies_expression_dataset.RData")
unscaled<-unscaled[rowSds(exprs(unscaled))>0,]
unscaled.s<-unscaled
exprs(unscaled.s)<-t(scale(t(exprs(unscaled)),scale=T))

exp<-read.table("~/CDC3/data/cheval_et_al_nephron_dataset.txt",sep="\t",header=T,as.is=T)
exp<-aggregate(as.matrix(exp[,2:9]),by=list(Symbol=exp$HUGO),sum)
rownames(exp)<-exp$Symbol
exp<-as.matrix(exp[,grep("Hs\\.",colnames(exp))])
exp<-exp[rowSds(exp)>0,]
exp.kid.s<-t(scale(t(exp),scale=T))

cg<-intersect(rownames(unscaled.s),rownames(exp.kid.s))
unscaled.s<-unscaled.s[cg,]

exp.kid.s<-exp.kid.s[cg,]
identical(rownames(unscaled.s),rownames(exp.kid.s))

correl.s<-cor(exprs(unscaled.s),exp.kid.s,method="pearson")
mat<-t(correl.s)

histology<-unscaled.s$histology
names(histology)<-colnames(unscaled.s)
histology<-sort(histology)
mat<-mat[,names(histology)]

df<-data.frame(Spearman=correl.s[,1],cluster=colnames(correl.s)[1],histology=unscaled.s$histology,stringsAsFactors = F)
for(i in 1:ncol(correl.s)){
  temp<-data.frame(Spearman=correl.s[,i],cluster=colnames(correl.s)[i],histology=unscaled.s$histology,stringsAsFactors = F)
  df<-rbind(df,temp)
}
df<-df[df$histology!="normal",]
df$histology[df$histology!="ccRCC"]<-capitalize(df$histology[df$histology!="ccRCC"])
df$histology<-factor(df$histology,levels=c("CDC","ccRCC","Papillary","Chromophobe","Oncocytoma","RMC","UTUC"))
df$cluster<-gsub("Hs.","",df$cluster)
df$cluster<-factor(df$cluster,levels=c("Glom","S1","S3","mTAL","cTAL","DCT","CCD","OMCD"))

ggplot(data=df, aes(x=cluster,y=Spearman, fill=cluster)) +
  geom_boxplot(coef=NULL) +
  scale_fill_simpsons() +
  theme_pubr(border=T,x.text.angle=90,legend="top") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size=12,face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(fill="Cell type") +
  ylab("Spearman's rho") +
  facet_wrap(~histology,ncol=4) +
  guides(fill=guide_legend(ncol=4))
ggsave("~/CDC3/results/supp_figure3/supplementary_figure3.tiff",width=8,height=5,dpi=600,compression="lzw")

