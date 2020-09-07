

library(Biobase)
library(limma)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(edgeR)
library(fgsea)
library(ggplotify)
library(cowplot)

dir.create("~/CDC3/results/figure1", recursive = T)
degsList <- vector("list", length = 4)
names(degsList) <- c("INT", "GSE89122", "WACH", "GSE11151")

#*************************
# Figure 1A
#*************************

load("~/CDC3/data/kidney_datasets/int/INT_processed.RData")
patient <- factor(dataset$Patient_ID)
f <- factor(dataset$histology)
design  <-  model.matrix(~ 0+f)
colnames(design) <- levels(f)
corfit  <-  duplicateCorrelation(dataset, design, block = patient)
corfit$consensus
fit  <-  lmFit(dataset, design, block = patient, correlation = corfit$consensus)
contrast.matrix  <-  makeContrasts(CDC-Normal, 
                                 CDC-ccRCC, 
                                 ccRCC-Normal, 
                                 levels = design)
fit2  <-  contrasts.fit(fit, contrast.matrix)
fit2  <-  eBayes(fit2)
degsList$INT <- topTable(fit2, coef = 1, number = nrow(dataset), adjust = "BH")

degsList$INT$Direction <- "N.S."
degsList$INT$Direction[degsList$INT$logFC> = 1 & degsList$INT$adj.P.Val<0.25] <- "Up-regulated"
degsList$INT$Direction[degsList$INT$logFC< = -1 & degsList$INT$adj.P.Val<0.25] <- "Down-regulated"
degsList$INT$Direction <- factor(degsList$INT$Direction, levels = c("N.S.", "Up-regulated", "Down-regulated"))
degsList$INT$Label <- ""
degsList$INT$Label[degsList$INT$logFC> = 1 & degsList$INT$adj.P.Val<0.25][1:10] <- degsList$INT$SYMBOL[degsList$INT$logFC> = 1 & degsList$INT$adj.P.Val<0.25][1:10]
degsList$INT$Label[degsList$INT$logFC< = -1 & degsList$INT$adj.P.Val<0.25][1:10] <- degsList$INT$SYMBOL[degsList$INT$logFC< = -1 & degsList$INT$adj.P.Val<0.25][1:10]

set.seed(85)
pC0 <- ggplot(data = degsList$INT, aes(x = logFC, y = -log10(degsList$INT$adj.P.Val), color = Direction, label = Label)) +
  geom_point() +
  theme_pubr(base_size = 15) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey40", "red2", "royalblue4")) +
  xlab(bquote(~Log[2]~ "fold change")) +
  ylab(bquote(~-Log[10]~ "false discovery rate")) +
  xlim(c(-4, 4)) +
  ylim(c(0, 4)) +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_text_repel(force = 20, size = 5, segment.alpha = 0.6)

#*************************
# Figure 1B
#*************************
source("~/CDC3/auxiliary_functions/enrichmentPlot_fgsea.R")

up.genes <- degsList$INT$SYMBOL[degsList$INT$logFC> = 1 & degsList$INT$adj.P.Val<0.25]
dn.genes <- degsList$INT$SYMBOL[degsList$INT$logFC< = -1 & degsList$INT$adj.P.Val<0.25]

setList <- list(up.genes, dn.genes)
names(setList) <- c("INT CDC vs normal UP", "INT CDC vs normal DN")

# GSE89122 - CDC vs NORMAL

load("~/CDC3/data/kidney_datasets/GSE89122/GSE89122_raw_counts.RData")
dataset <- dataset[, dataset$patient ! =  "CDC5"]
dge <- DGEList(counts = exprs(dataset), genes = fData(dataset))
isexpr <- rowSums(exprs(dataset)>10) > =  1
dge <- dge[isexpr, ]
dge  <-  calcNormFactors(dge, method = "TMM")
f <- factor(dataset$histology)
patient <- factor(dataset$patient)
design <- model.matrix(~0+f+patient)
v <- voom(dge, design = design, plot = F)
fit  <-  lmFit(v, design)
contrast.matrix  <-  makeContrasts(fCDC-fNormal, levels = design)
fit2  <-  contrasts.fit(fit, contrast.matrix)
fit2  <-  eBayes(fit2)
degsList$GSE89122 <- topTable(fit2, number = nrow(dataset), adjust = "BH")
rnk <- degsList$GSE89122$t
names(rnk) <- degsList$GSE89122$Symbol
set.seed(1234)
fgseaRes <- fgsea(pathways = setList, stats = rnk, nperm = 10000)
pC1 <- as.ggplot(~plotEnrichment.full(pathway = c("INT CDC vs normal UP", "INT CDC vs normal DN"), 
                                  contrast.name = "GSE89122", 
                                  collection = setList, 
                                  stats = rnk, 
                                  enrichment.resObj = fgseaRes, 
                                  cex.legend = 1.8, 
                                  cex.yaxis = 1.2, 
                                  cex.xaxis = 1.8, 
                                  cex.lab = 1.5, 
                                  cex.geneset.name = 2, 
                                  positive.name = "Up", 
                                  negative.name = "Down", 
                                  plotRankMetric = F, 
                                  matrix.layout = c(1:8), 
                                  matrix.layout.ncol = 1, 
                                  heights = c(2.0, 0.41, 0.35, 2), 
                                  inset.legend = c(-0.28, 0), 
                                  change.hit.col = F))

# WACH - CDC vs NORMAL

load("~/CDC3/data/kidney_datasets/wach/wach_raw_counts.RData")
dge <- DGEList(counts = exprs(dataset), genes = fData(dataset))
isexpr <- rowSums(exprs(dataset)>10) > =  1
dge <- dge[isexpr, ]
dge  <-  calcNormFactors(dge, method = "TMM")
f <- factor(dataset$histology)
design <- model.matrix(~0+f)
v <- voom(dge, design = design, plot = F)
fit  <-  lmFit(v, design)
contrast.matrix  <-  makeContrasts(fCDC-fNormal, levels = design)
fit2  <-  contrasts.fit(fit, contrast.matrix)
fit2  <-  eBayes(fit2)
degsList$WACH <- topTable(fit2, number = nrow(dataset), adjust = "BH")
rnk <- degsList$WACH$t
names(rnk) <- degsList$WACH$Symbol
set.seed(1234)
fgseaRes <- fgsea(pathways = setList, stats = rnk, nperm = 10000)
pC2 <- as.ggplot(~plotEnrichment.full(pathway = c("INT CDC vs normal UP", "INT CDC vs normal DN"), 
                                  contrast.name = "WACH", 
                                  collection = setList, 
                                  stats = rnk, 
                                  enrichment.resObj = fgseaRes, 
                                  cex.legend = 1.8, 
                                  cex.yaxis = 1.2, 
                                  cex.xaxis = 1.8, 
                                  cex.lab = 1.5, 
                                  cex.geneset.name = 2, 
                                  positive.name = "Up", 
                                  negative.name = "Down", 
                                  plotRankMetric = F, 
                                  matrix.layout = c(1:8), 
                                  matrix.layout.ncol = 1, 
                                  heights = c(2.0, 0.41, 0.35, 2), 
                                  inset.legend = c(-0.28, 0), 
                                  change.hit.col = F))

# GSE11151 - CDC vs NORMAL

load("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151/GSE11151_series_matrix_collMaxMean.RData")
dataset.coll <- dataset.coll[, dataset.coll$tumor_type%in%c("adult normal kidney", "collecting duct carcinoma")]

f <- factor(ifelse(dataset.coll$tumor_type  ==  "collecting duct carcinoma", "CDC", "Normal"))
design  <-  model.matrix(~ 0+f)
fit  <-  lmFit(dataset.coll, design)
contrast.matrix  <-  makeContrasts(fCDC-fNormal, levels = design)
fit2  <-  contrasts.fit(fit, contrast.matrix)
fit2  <-  eBayes(fit2)
degsList$GSE11151 <- topTable(fit2, number = nrow(dataset.coll), adjust = "BH")
rnk <- degsList$GSE11151$t
names(rnk) <- degsList$GSE11151$SYMBOL
set.seed(1234)
fgseaRes <- fgsea(pathways = setList, stats = rnk, nperm = 10000)

pC3 <- as.ggplot(~plotEnrichment.full(pathway = c("INT CDC vs normal UP", "INT CDC vs normal DN"), 
                                  contrast.name = "GSE11151", 
                                  collection = setList, 
                                  stats = rnk, 
                                  enrichment.resObj = fgseaRes, 
                                  cex.legend = 1.8, 
                                  cex.yaxis = 1.2, 
                                  cex.xaxis = 1.8, 
                                  cex.lab = 1.5, 
                                  cex.geneset.name = 2, 
                                  positive.name = "Up", 
                                  negative.name = "Down", 
                                  plotRankMetric = F, 
                                  matrix.layout = c(1:8), 
                                  matrix.layout.ncol = 1, 
                                  heights = c(2.0, 0.41, 0.35, 2), 
                                  inset.legend = c(-0.28, 0), 
                                  change.hit.col = F))

plotlist <- list(pC0, pC1, pC2, pC3)
tiff("~/CDC3/results/figure1/figure1.tiff", width = 16, height = 6, units = "in", res = 600, compression = "lzw")
cowplot::plot_grid(plotlist = plotlist, ncol = 4, rel_widths = c(1, 0.7, 0.7, 0.7), labels = c("A", "B", "", ""))
dev.off()



