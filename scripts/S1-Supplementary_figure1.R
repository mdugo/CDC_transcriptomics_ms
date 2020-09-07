
library(gprofiler2)
library(ComplexHeatmap)

dir.create("~/CDC3/results/supp_figure1")

gs<-scan("~/CDC3/results/figure3/INT_CDC_signature_gene_list.txt",what="character")

res<-gost(query=gs,
              evcodes=T,
              correction_method = "g_SCS",
              sources=c("GO:BP","REAC"))$result

res$term_name<-paste(res$source,capitalize(res$term_name))

genes<-unique(unlist(strsplit(res$intersection,",")))
mat<-matrix(0,ncol=nrow(res),nrow=length(genes))
colnames(mat)<-res$term_name
rownames(mat)<-genes
for(i in 1:ncol(mat)){
  mat[,i]<-ifelse(rownames(mat)%in%unlist(strsplit(res$intersection[i],",")),1,0)
}
mat<-t(mat)
rownames(mat)[8]<-"GO:BP Cell diff. involved in embryonic placenta development"

ha<-rowAnnotation("-log10(FDR)"=anno_barplot(-log10(res$p_value),width = unit(2, "cm")),show_annotation_name = T,annotation_name_gp = gpar(fontsize=9))

tiff("~/CDC3/results/supp_figure1/supplementary_figure1.tiff",width=10,height=6, units = "in", res=600, compression="lzw")
Heatmap(mat,
        col=c("white","black"),
        border=T,
        rect_gp = gpar(col = "grey70", lwd = 1),
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        show_column_dend = F,
        show_row_dend = F,
        show_heatmap_legend = F,
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 90,
        right_annotation = ha,
        heatmap_width = unit(6, "inches"),
        heatmap_height = unit(5, "inches"),
        row_names_max_width = unit(6, "cm"))
dev.off()
