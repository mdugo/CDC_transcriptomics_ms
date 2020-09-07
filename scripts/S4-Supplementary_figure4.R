
library(Biobase)
library(ggplot2)
library(ggpubr)
library(cola)
library(singscore)
library(RColorBrewer)

dir.create("~/CDC3/results/supp_figure4")

#**************************************
# Supplementary Figure 4A
#**************************************

load("~/CDC3/results/figure7/cola_subgroups_CDC.RData")
tiff("~/CDC3/results/supp_figure4/supplementary_figure4A.tiff",width=8,height=7,units="in",res=600,compression="lzw")
classes<-collect_classes(rl, k = 2)
dev.off()

#**************************************
# Supplementary Figure 4B
#**************************************

load("~/CDC3/results/figure7/cola_subgroups_CDC.RData")
x<-collect_classes(rl, k = 2)
load("~/CDC3/data/kidney_datasets/metadataset_CDC_zscore_quantile.RData")
dataset.zscore$class<-ifelse(get_classes(rl,k=2)$class==1,"C2","C1")

gs<-scan("~/CDC3/results/figure3/INT_CDC_signature_gene_list.txt",what="character")

rankData<-rankGenes(dataset.zscore)
score<-simpleScore(rankData,upSet=gs)

df<-data.frame(Score=score$TotalScore,Cluster=ifelse(dataset.zscore$class=="C1","CDC-S1","CDC-S2"))
ggplot(data=df, aes(x=Cluster,y=Score, fill=Cluster)) +
  geom_boxplot(coef=NULL) +
  theme_pubr(legend="none") +
  stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y = -0.05) +
  ylab("INT-CDC signature score") +
  scale_fill_manual(values=brewer.pal(3,"Set2")[1:2]) +
  theme(text = element_text(size=10))
ggsave("~/CDC3/results/supp_figure4/supplementary_figure4B.tiff",width=7.5,height=7.5,units = "cm",dpi=600,compression="lzw")
