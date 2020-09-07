library(GEOquery)
library(WGCNA)
library(DESeq2)
library(readxl)
library(utils)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(Seurat)
library(harmony)
library(filesstrings)
library(affy)
library(frma)
library(affycoretools)
library(hgu133plus2.db)
library(hgu133a.db)
library(reshape2)
library(CePa)
library(ArrayExpress)
library(hgu219.db)
library(goseq)
library(GenomicFeatures)
library(illuminaHumanv4.db)
library(oligo)
library(lumi)
library(edgeR)

dir.create("~/CDC3/data/kidney_datasets", recursive=T)

######## download annotation for TCGA datasets
download.file("https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f",destfile="~/CDC3/data/tcga_gencode.v22.annotation.gtf.gz")
gunzip(path.expand("~/CDC3/data/tcga_gencode.v22.annotation.gtf.gz"))

########  download GENCODE gene annotation 
download.file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz",destfile="~/CDC3/data/gencode.v32lift37.annotation.gtf.gz")
gunzip(path.expand("~/CDC3/data/gencode.v32lift37.annotation.gtf.gz"))

annot<-read.table("~/CDC3/data/gencode.v32lift37.annotation.gtf",sep="\t",header=F,as.is=T,skip=5)
annot<-annot[annot$V3=="gene",]
annot2<-strsplit(annot$V9,"; ")
annot$gene_id<-gsub(".* ","",sapply(annot2,function(x)x[grep("gene_id",x)]))
annot$gene_type<-gsub(".* ","",sapply(annot2,function(x)x[grep("gene_type",x)]))
annot$gene_name<-gsub(".* ","",sapply(annot2,function(x)x[grep("gene_name",x)]))
annot<-annot[,-c(1:9)]
annot$gene_length<- getlength(annot$gene_name, genome = "hg19",id = "geneSymbol")
write.table(annot,file="~/CDC3/data/gene_annotation_gencode32_hg19.txt",sep="\t",row.names=F,quote=F)
file.remove("~/CDC3/data/gencode.v32lift37.annotation.gtf")

########  generate REACTOME hierarchy for gene sets categorization
# The REACTOME file describing relations between pathways was downloaded from 
download.file("https://reactome.org/download/current/ReactomePathways.txt",destfile = "~/CDC3/data/ReactomePathways.txt")
download.file("https://reactome.org/download/current/ReactomePathwaysRelation.txt", destfile="~/CDC3/data/ReactomePathwaysRelation.txt")

pathway<-read.table("~/CDC3/data/ReactomePathwaysRelation.txt",sep="\t",as.is=T)
colnames(pathway)<-c("Parent","Son")
pathway<-pathway[grep("HSA-",pathway$Parent),]

annot<-read.table("~/CDC3/data/ReactomePathways.txt",sep="\t",as.is=T,quote="")
annot<-annot[annot$V3=="Homo sapiens",1:2]
colnames(annot)<-c("Parent","Name")
pathway<-merge(pathway,annot,by="Parent")
annot$Son<-annot$Parent
annot<-annot[,colnames(annot)!="Parent"]
pathway<-merge(pathway,annot,by="Son",suffixes=c("_Parent","_Son"))
pathway<-pathway[,c(2,1,3,4)]
pathway<-pathway[order(pathway$Parent,pathway$Son),]
pathway$Name_Parent<-gsub(";","",pathway$Name_Parent)
pathway$Name_Son<-gsub(";","",pathway$Name_Son)

pathway$hierarchy<-paste(pathway$Name_Parent,pathway$Name_Son,sep=";")
pathway$hierarchy2<-paste(pathway$Name_Parent,pathway$Name_Son,sep=";")
pathway$hierarchy3<-paste(pathway$Name_Parent,pathway$Name_Son,sep=";")
pathway$hierarchy4<-paste(pathway$Name_Parent,pathway$Name_Son,sep=";")
pathway$hierarchy5<-paste(pathway$Name_Parent,pathway$Name_Son,sep=";")
pathway$hierarchy6<-paste(pathway$Name_Parent,pathway$Name_Son,sep=";")

for(i in 1:nrow(pathway)){
  print(i)
  while(strsplit(pathway$hierarchy[i],";")[[1]][1]%in%pathway$Name_Son){
    parents<-pathway$Name_Parent[pathway$Name_Son==strsplit(pathway$hierarchy[i],";")[[1]][1]][1]
    pathway$hierarchy[i]<-paste(parents,pathway$hierarchy[i],sep=";")
  }
}

for(i in 1:nrow(pathway)){
  print(i)
  while(strsplit(pathway$hierarchy2[i],";")[[1]][1]%in%pathway$Name_Son){
    parents<-pathway$Name_Parent[pathway$Name_Son==strsplit(pathway$hierarchy2[i],";")[[1]][1]]
    if(length(parents)>1){
      parents<-parents[2]
    }
    pathway$hierarchy2[i]<-paste(parents,pathway$hierarchy2[i],sep=";")
  }
}

for(i in 1:nrow(pathway)){
  print(i)
  while(strsplit(pathway$hierarchy3[i],";")[[1]][1]%in%pathway$Name_Son){
    parents<-pathway$Name_Parent[pathway$Name_Son==strsplit(pathway$hierarchy3[i],";")[[1]][1]]
    if(length(parents)>1){
      parents<-parents[3]
    }
    pathway$hierarchy3[i]<-paste(parents,pathway$hierarchy3[i],sep=";")
  }
}

for(i in 1:nrow(pathway)){
  print(i)
  while(strsplit(pathway$hierarchy4[i],";")[[1]][1]%in%pathway$Name_Son){
    parents<-pathway$Name_Parent[pathway$Name_Son==strsplit(pathway$hierarchy4[i],";")[[1]][1]]
    if(length(parents)>1){
      parents<-parents[4]
    }
    pathway$hierarchy4[i]<-paste(parents,pathway$hierarchy4[i],sep=";")
  }
}

for(i in 1:nrow(pathway)){
  print(i)
  while(strsplit(pathway$hierarchy5[i],";")[[1]][1]%in%pathway$Name_Son){
    parents<-pathway$Name_Parent[pathway$Name_Son==strsplit(pathway$hierarchy5[i],";")[[1]][1]]
    if(length(parents)>1){
      parents<-parents[5]
    }
    pathway$hierarchy5[i]<-paste(parents,pathway$hierarchy5[i],sep=";")
  }
}

for(i in 1:nrow(pathway)){
  print(i)
  while(strsplit(pathway$hierarchy6[i],";")[[1]][1]%in%pathway$Name_Son){
    parents<-pathway$Name_Parent[pathway$Name_Son==strsplit(pathway$hierarchy6[i],";")[[1]][1]]
    if(length(parents)>1){
      parents<-parents[6]
    }
    pathway$hierarchy6[i]<-paste(parents,pathway$hierarchy6[i],sep=";")
  }
}

pathway$hierarchy2[pathway$hierarchy2==pathway$hierarchy]<-""
pathway$hierarchy3[pathway$hierarchy3==pathway$hierarchy]<-""
pathway$hierarchy4[pathway$hierarchy4==pathway$hierarchy]<-""
pathway$hierarchy5[pathway$hierarchy5==pathway$hierarchy]<-""
pathway$hierarchy6[pathway$hierarchy6==pathway$hierarchy]<-""

pathway$hierarchy2[grep("^NA;",pathway$hierarchy2)]<-""
pathway$hierarchy3[grep("^NA;",pathway$hierarchy3)]<-""
pathway$hierarchy4[grep("^NA;",pathway$hierarchy4)]<-""
pathway$hierarchy5[grep("^NA;",pathway$hierarchy5)]<-""
pathway$hierarchy6[grep("^NA;",pathway$hierarchy6)]<-""

hierarchies<-c(pathway$hierarchy,pathway$hierarchy2,pathway$hierarchy3,pathway$hierarchy4,pathway$hierarchy5,pathway$hierarchy6)
hierarchies<-hierarchies[hierarchies!=""]
hierarchies<-unique(hierarchies)

df<-data.frame(Hierarchy=hierarchies,Hierarchies_GSEA=tolower(hierarchies))
df$Hierarchies_GSEA<-gsub("\\.","",gsub("\\'","",gsub("\\+","",gsub("\\&","",gsub("-"," ",gsub("\\)","",gsub("\\(","",gsub("\\/"," ",gsub(",","",gsub("\\:"," ",df$Hierarchies_GSEA))))))))))
df$Hierarchies_GSEA2<-gsub(" ","",df$Hierarchies_GSEA)
write.table(df,file="~/CDC3/data/reactome_hierarchy.txt",sep="\t",row.names=F,quote=F)

file.remove("~/CDC3/data/ReactomePathways.txt")
file.remove("~/CDC3/data/ReactomePathwaysRelation.txt")

#***************************
# INT dataset
#***************************

dir.create("~/CDC3/data/kidney_datasets/int")

gse153965<-getGEO("GSE153965")
x<-gse153965$GSE153965_series_matrix.txt.gz
pdata<-pData(gse153965)
pdata<-pdata[,c(1,2,8,10)]
pdata$patient<-gsub(".* ","",pdata$title)
pdata$histology<-"CDC"
pdata$histology[grep("Normal",pdata$title)]<-"Normal"
pdata$tumor_grade<-gsub("sample","NA",gsub(".* ","",pdata$characteristics_ch1))

# download cel files
sapply(colnames(gse153965),getGEOSuppFiles, makeDirectory = F,baseDir = "~/CDC3/data/kidney_datasets/int")
celFiles<-list.celfiles("~/CDC3/data/kidney_datasets/int", full.names=TRUE)
pdata<-new("AnnotatedDataFrame",data=pdata)
rawData<-read.celfiles(celFiles,phenoData=pdata)
dataset<-rma(rawData,normalize=TRUE,target="core")


## NetAffx annotation
unzip(zipfile="~/CDC3/data/kidney_datasets/int/HTA-2_0.r3.na36.hg19.a1.transcript.csv.zip",overwrite = T)

annot<-read.csv("~/CDC3/data/kidney_datasets/int/HTA-2_0.r3.na36.hg19.a1.transcript.csv",sep=",",dec=".",header=T,as.is=T,comment.char="#")
annot<-annot[,-c(10:16,19:21)]
rownames(annot)<-annot$transcript_cluster_id
annot<-annot[rownames(dataset),]
identical(rownames(annot),rownames(dataset))

gene_list<-strsplit(annot$gene_assignment," /// ")
annot$affy_mrna<-sapply(gene_list,function(y)paste(unique(sapply(sapply(y,function(x)strsplit(x," // ")),"[",1)),collapse=";"))
annot$affy_symbol<-sapply(gene_list,function(y)paste(unique(sapply(sapply(y,function(x)strsplit(x," // ")),"[",2)),collapse=";"))
annot$affy_entrez<-sapply(gene_list,function(y)paste(unique(sapply(sapply(y,function(x)strsplit(x," // ")),"[",5)),collapse=";"))
annot$SYMBOL<-sapply(strsplit(annot$affy_symbol,";"),"[",1)
annot$ENTREZ<-sapply(strsplit(annot$affy_entrez,";"),"[",1)

annot<-annot[rownames(dataset),]
identical(rownames(annot),rownames(dataset))

fData(dataset)<-annot[,c(1,10,11,15,16)]

####filtering

dataset<-dataset[fData(dataset)$SYMBOL!="NA",]

###collapse
x<-collapseRows(exprs(dataset),
                rowID=rownames(dataset),
                rowGroup = fData(dataset)$SYMBOL,
                method="maxRowVariance"
)
dataset<-dataset[x$selectedRow,]
save(dataset,file="~/CDC3/data/kidney_datasets/int/INT_processed.RData")

#***************************
# GSE89122
#***************************

dir.create("~/CDC3/data/kidney_datasets/GSE89122")

# get series matrix to obtain sample names and phenodata

gse89122 <- getGEO(GEO = "GSE89122")
gse89122<-gse89122$GSE89122_series_matrix.txt.gz
pdata<-pData(gse89122)
pdata<-pdata[,c(1,2,8,10)]
pdata$patient<-gsub(".* ","",pdata$title)
pdata$histology<-"CDC"
pdata$histology[grep("Normal",pdata$title)]<-"Normal"
pdata$tumor_grade<-gsub("sample","NA",gsub(".* ","",pdata$characteristics_ch1))

# download raw counts
sapply(colnames(gse89122),getGEOSuppFiles, makeDirectory = F,baseDir = "~/CDC3/data/kidney_datasets/GSE89122")
files<-list.files("~/CDC3/data/kidney_datasets/GSE89122",full.names=T,"count.txt.gz")
sapply(files,gunzip,overwrite=T)

files<-list.files("~/CDC3/data/kidney_datasets/GSE89122",full.names=T,"count.txt")
expdata<-data.frame(NA)

for(i in files){
  temp<-read.table(i,sep="\t",header=T,quote="",row.names = 1)
  expdata<-cbind(expdata,temp)
}
expdata<-expdata[,-1]
colnames(expdata)<-gsub("_.*","",gsub(".*\\/","",files))
fdata<-data.frame(Symbol=rownames(expdata),stringsAsFactors = F)

annot<-read.table("~/CDC3/data/gene_annotation_gencode32_hg19.txt",header=T,sep="\t",as.is=T)
annot<-annot[,c(3:4)]
annot<-unique(annot)
annot<-aggregate(annot$gene_length,by=list(Symbol=annot$gene_name),max)
colnames(annot)[2]<-"gene_length"
fdata<-merge(fdata,annot,by="Symbol",all.x=T)
rownames(fdata)<-fdata$Symbol
fdata<-fdata[rownames(expdata),]
identical(rownames(fdata),rownames(expdata))

dataset<-ExpressionSet(assayData=as.matrix(expdata),
                       phenoData=new("AnnotatedDataFrame",pdata),
                       featureData=new("AnnotatedDataFrame",fdata)
)
dataset<-dataset[-grep("^__",rownames(dataset)),]
save(dataset,file="~/CDC3/data/kidney_datasets/GSE89122/GSE89122_raw_counts.RData")
file.remove(list.files("~/CDC3/data/kidney_datasets/GSE89122",full.names=T,"count.txt"))

#***************************
# WACH
#***************************

dir.create("~/CDC3/data/kidney_datasets/wach")

download.file("https://www.mdpi.com/2072-6694/12/1/64/s1", destfile="~/CDC3/data/kidney_datasets/wach/cancers-12-00064-s001.zip")
unzip(path.expand("~/CDC3/data/kidney_datasets/wach/cancers-12-00064-s001.zip"), overwrite=T, unzip="unzip", exdir = path.expand("~/CDC3/data/kidney_datasets/wach"))
file.remove("~/CDC3/data/kidney_datasets/wach/cancers-12-00064-s001.zip")

exp<-as.data.frame(read_xlsx("~/CDC3/data/kidney_datasets/wach/cancers-667778-suppl-XML/cancers-667778-suppl-Tables.xlsx", sheet=11))
exp<-exp[-1,]
colnames(exp)[1]<-"gene_id"
for(i in 2:ncol(exp)){
  exp[,i]<-as.numeric(exp[,i])
}

annot<-read.table("~/CDC3/data/gene_annotation_gencode32_hg19.txt",header=T,sep="\t",as.is=T)
annot2<-annot[,c(1,3)]
annot2$gene_id<-gsub("\\..*","",annot2$gene_id)
annot2<-unique(annot2)
exp<-merge(annot2,exp,by="gene_id",all.y=T)
exp<-aggregate(as.matrix(exp[,-c(1:2)]),by=list(Symbol=exp$gene_name),sum)

fdata<-data.frame(Symbol=exp$Symbol, stringsAsFactors = F)
annot<-annot[,c(3:4)]
annot<-aggregate(annot$gene_length,by=list(Symbol=annot$gene_name),max)
colnames(annot)[2]<-"gene_length"
fdata<-merge(fdata,annot,by="Symbol",all.x=T)
rownames(fdata)<-fdata$Symbol
rownames(exp)<-exp$Symbol
exp<-as.matrix(exp[,-1])
fdata<-fdata[rownames(exp),]
identical(rownames(fdata),rownames(exp))

pdata<-data.frame(sample_ID=colnames(exp),histology=c(rep("Normal",8),rep("CDC",2)))
rownames(pdata)<-colnames(exp)

dataset<-ExpressionSet(assayData = as.matrix(exp),
                       phenoData = new("AnnotatedDataFrame",pdata),
                       featureData = new("AnnotatedDataFrame",fdata))
save(dataset,file="~/CDC3/data/kidney_datasets/wach/wach_raw_counts.RData")

#***************************
# GEO dataset
#***************************

dir.create("~/CDC3/data/kidney_datasets/ncbi_geo/",recursive=T)

mydatasets<-c("GSE2748","GSE11151","GSE11904","GSE12090","GSE14994","GSE17895","GSE19949","GSE22541",
              "GSE23629","GSE26574","GSE29609","GSE36895","GSE40435","GSE41137","GSE46699","GSE47032",
              "GSE53757","GSE65615","GSE68417")
for(i in mydatasets){
  dir.create(paste0("~/CDC3/data/kidney_datasets/ncbi_geo/",i))
}

for(i in mydatasets){
  print(paste0("Downloading dataset ",i))
  x<-getGEO(i,getGPL = T)
  if(class(x)=="list" & length(x)==1){
    x<-x[[1]]
  }
  save(x,file=paste0("~/CDC3/data/kidney_datasets/ncbi_geo/",i,"/",i,"_series_matrix.RData"))
}

files<-list.files("~/CDC3/data/kidney_datasets/ncbi_geo/","_series_matrix.RData",recursive=T,full.names=T)

#-------------------------------------------------------------------------------------------
load(files[1])
pData(x)<-pData(x)[,c(1,2,10)]
x$tissue<-ifelse(x$characteristics_ch1=="renal tumor","tumor","normal")
x$tumor_type<-gsub(";.*","",x$title)
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\..*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[1]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[2])
pData(x)<-pData(x)[,c(1,2,40:47)]
x$tissue<-"tumor"
x$tumor_type<-x$`Disease State:ch1`
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE11904")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11904")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11904"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\..*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133a2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[2]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11904",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11904","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[3])
pData(x)<-pData(x)[,c(1,2,10)]
x$tissue<-"tumor"
x$tumor_type<-gsub("\\:.*","",x$title)
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE12090")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE12090")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE12090"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\..*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[3]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE12090",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE12090","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[4])
x<-x[[2]]
pData(x)<-pData(x)[,c(1,2,33:35)]
x$tissue<-"tumor"
x$tissue[grep("normal",x$`disease state:ch1`)]<-"normal"
x$tumor_type<-x$`disease state:ch1`
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE14994")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE14994")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE14994"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\..*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hthgu133a.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[4]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE14994",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE14994","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[5])
pData(x)<-pData(x)[,c(1,2,52)]
x$tissue<-ifelse(x$`subtype:ch1`=="normal tissue","normal","tumor")
x$tumor_type<-x$`subtype:ch1`
colnames(fData(x))[3]<-"Symbol"
fData(x)$Symbol<-gsub(" .*","",fData(x)$Symbol)
sum(is.na(fData(x)$Symbol) | fData(x)$Symbol=="")
x<-x[!is.na(fData(x)$Symbol) & fData(x)$Symbol!="",]
collrow<-collapseRows(exprs(x),
                      rowID=rownames(x),
                      rowGroup=fData(x)$Symbol,
                      method="MaxMean")
dataset.coll<-x[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[5]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE17895","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[6])
x<-x[[1]]
x<-x[,is.na(x$`cell line:ch1`)]
x<-x[,-grep("mixed|metastasis",x$`tissue type:ch1`)]
pData(x)<-pData(x)[,c(1,2,42:50)]
x$tissue<-"tumor"
x$tumor_type<-x$`icd-o 3 diagnosis text:ch1`
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE19949")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE19949")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE19949"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\_.*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hthgu133a.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[6]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE19949",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE19949","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[7])
pData(x)<-pData(x)[,c(1,2,8,10:17)]
x$tissue<-"tumor"
x$tumor_type<-as.vector(x$source_name_ch1)
x<-x[,-grep("metastasis",x$tumor_type)]
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE22541")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE22541")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE22541"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\..*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[7]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE22541",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE22541","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[8])
pData(x)<-pData(x)[,c(1,2,36,37)]
x$tissue<-x$`tissue:ch1`
x$tumor_type<-"clear cell renal cell carcinoma"
x<-x[,-grep("Metast",x$tissue)]
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE23629")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE23629")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE23629"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\_.*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[8]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE23629",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE23629","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[9])
pData(x)<-pData(x)[,c(1,2,30)]
x$tissue<-"tumor"
x$tissue[grep("normal",x$`disease state:ch1`)]<-"normal"
x$tumor_type<-x$`disease state:ch1`
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE26574")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE26574")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE26574"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\..*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[9]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE26574",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE26574","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[10])
pData(x)<-pData(x)[,c(1,2,27)]
x$tissue<-"tumor"
x$tumor_type<-"papillary renal cell carcinoma"

x<-annotateEset(x,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
x<-x[!is.na(fData(x)$SYMBOL),]

collrow<-collapseRows(exprs(x),
                      rowID=rownames(x),
                      rowGroup=fData(x)$SYMBOL,
                      method="MaxMean")
dataset.coll<-x[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[10]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE2748","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[11])
pData(x)<-pData(x)[,c(1,2,58:76)]
x$tissue<-"tumor"
x$tumor_type<-"clear cell renal cell carcinoma"
fData(x)<-fData(x)[,c(5,9,10)]
colnames(fData(x))[3]<-"Symbol"
sum(is.na(fData(x)$Symbol) | fData(x)$Symbol=="")
x<-x[!is.na(fData(x)$Symbol) & fData(x)$Symbol!="",]
collrow<-collapseRows(exprs(x),
                      rowID=rownames(x),
                      rowGroup=fData(x)$Symbol,
                      method="MaxMean")
dataset.coll<-x[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[11]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE29609","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[12])
pData(x)<-pData(x)[,c(1,2,44:54)]
x<-x[,-grep("Mouse tumorgraft",x$`tissue:ch1`)]
x$tissue<-ifelse(x$`tissue:ch1`=="Normal cortex","normal","tumor")
x$tumor_type<-ifelse(x$`tissue:ch1`=="Normal cortex","normal","clear cell renal cell carcinoma")
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE36895")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE36895")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE36895"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\_.*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[12]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE36895",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE36895","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[13])
pData(x)<-pData(x)[,c(1,2,22,35:39)]
x$tissue<-ifelse(x$`tissue type:ch1`=="adjacent non-tumour renal tissue","normal","tumor")
x$tumor_type<-ifelse(x$`tissue type:ch1`=="adjacent non-tumour renal tissue","normal","clear cell renal cell carcinoma")
pdata<-pData(x)

getGEOSuppFiles("GSE40435",makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE40435")
gunzip("~/CDC3/data/kidney_datasets/ncbi_geo/GSE40435/GSE40435_non_normalized.txt.gz")
raw<-read.table("~/CDC3/data/kidney_datasets/ncbi_geo/GSE40435/GSE40435_non_normalized.txt",header=T,as.is=T,sep="\t",row.names=1)
exp.raw<-raw[,grep("Signal",colnames(raw))]
colnames(exp.raw)<-gsub("\\..*","",gsub("X","",colnames(exp.raw)))
det.raw<-raw[,grep("Detection",colnames(raw))]
colnames(det.raw)<-gsub("\\..*","",gsub("X","",colnames(det.raw)))
identical(colnames(exp.raw),as.vector(pdata$description))
colnames(exp.raw)<-rownames(pdata)
colnames(det.raw)<-rownames(pdata)

fdata<-data.frame(ILMN_ID=rownames(exp.raw),stringsAsFactors=F)
fdata$Symbol<-unlist(mget(fdata$ILMN_ID,illuminaHumanv4SYMBOL))
rownames(fdata)<-fdata$ILMN_ID

dataset<-ExpressionSet(assayData = as.matrix(exp.raw),
                       phenoData=new("AnnotatedDataFrame",pdata),
                       featureData=new("AnnotatedDataFrame",fdata))
detection(dataset)<-det.raw
dataset<-lumiN(lumiT(dataset,method="log2"),method="rsn")

dataset<-dataset[rowSums(detection(dataset)<0.01)>0,]
dataset<-dataset[!is.na(fData(dataset)$Symbol),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$Symbol,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[13]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE40435","non_normalized.txt",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE40435","RAW.tar",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE40435","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[14])
x<-x[,x$`treatment:ch1`=="Snap Frozen"]
pData(x)<-pData(x)[,c(1,2,36:38)]
x$tissue<-"tumor"
x$tumor_type<-x$`tissue:ch1`
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE41137")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE41137")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE41137"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\_.*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[14]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE41137",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE41137","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[15])
pData(x)<-pData(x)[,c(1,2,36:41)]
x$tissue<-x$`tissue:ch1`
x$tumor_type<-"clear cell renal cell carcinoma"
pdata<-pData(x)
sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE46699")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE46699")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE46699"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\_.*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[15]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE46699",".cel.gz",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE46699","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[16])
x<-x[[1]]
pData(x)<-pData(x)[,c(1,2,38:43)]
x$tissue<-ifelse(x$`tissue:ch1`=="ccRCC","tumor","normal")
x$tumor_type<-ifelse(x$`tissue:ch1`=="ccRCC","clear cell renal cell carcinoma","normal tissue")
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE47032")

gzfiles<-list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE47032",".CEL.gz",full.names=T)
sapply(gzfiles,gunzip)
cel<-list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE47032",".CEL",full.names=T)

raw<-read.celfiles(filenames = cel,
                   sampleNames = rownames(pdata))
dataset<-frma(raw,target="core")
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="huex10sttranscriptcluster.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[16]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE47032",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE47032","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[17])
pData(x)<-pData(x)[,c(1,2,33:34)]
x$tissue<-ifelse(x$`tissue:ch1`=="clear cell renal cell carcinoma","tumor","normal")
x$tumor_type<-x$`tissue:ch1`
pdata<-pData(x)
sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE53757")

cel<-list.celfiles("~/CDC3/data/kidney_datasets/ncbi_geo/GSE53757")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/kidney_datasets/ncbi_geo/GSE53757"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\_.*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[17]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE53757",".cel.gz",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE53757","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[18])
pData(x)<-pData(x)[,c(1,2,35:37)]
x<-x[,x$`diagnosis:ch1`!="--"]
x$tissue<-"tumor"
x$tumor_type<-x$`diagnosis:ch1`
pdata<-pData(x)
pdata$title<-as.vector(pdata$title)
pdata$title<-gsub("\\.","",gsub("\\/",".",gsub(" ","",pdata$title)))

getGEOSuppFiles("GSE65615",makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE65615")
gunzip("~/CDC3/data/kidney_datasets/ncbi_geo/GSE65615/GSE65615_RawData.txt.gz")
raw<-read.table("~/CDC3/data/kidney_datasets/ncbi_geo/GSE65615/GSE65615_RawData.txt",header=T,as.is=T,sep="\t",row.names=1,quote="",comment.char="")

exp.raw<-raw[,grep("Signal",colnames(raw))]
colnames(exp.raw)<-gsub("\\.","",gsub("^X","",gsub("\\.AVG.*","",gsub("\\..AVG.*","",colnames(exp.raw)))))
det.raw<-raw[,grep("Detection",colnames(raw))]
colnames(det.raw)<-colnames(exp.raw)
exp.raw<-exp.raw[,pdata$title]
det.raw<-det.raw[,pdata$title]
identical(colnames(exp.raw),as.vector(pdata$title))
colnames(exp.raw)<-rownames(pdata)
colnames(det.raw)<-rownames(pdata)

fdata<-data.frame(ILMN_ID=rownames(exp.raw),stringsAsFactors=F)
fdata$Symbol<-unlist(mget(fdata$ILMN_ID,illuminaHumanv4SYMBOL))
rownames(fdata)<-fdata$ILMN_ID

dataset<-ExpressionSet(assayData = as.matrix(exp.raw),
                       phenoData=new("AnnotatedDataFrame",pdata),
                       featureData=new("AnnotatedDataFrame",fdata))
detection(dataset)<-det.raw
dataset<-lumiN(lumiT(dataset,method="log2"),method="rsn")

dataset<-dataset[rowSums(detection(dataset)<0.01)>0,]
dataset<-dataset[!is.na(fData(dataset)$Symbol),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$Symbol,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[18]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE65615","RawData.txt",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE65615","RAW.tar",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE65615","series_matrix.RData",full.names=T))

#-------------------------------------------------------------------------------------------
load(files[19])
pData(x)<-pData(x)[,c(1,2,39:42)]
x$tissue<-"tumor"
x$tissue[x$`analysis group:ch1`=="Benign"]<-"benign"
x$tissue[x$`analysis group:ch1`=="normal"]<-"normal"
x$tumor_type<-x$`tissue:ch1`
x<-x[,!x$tumor_type%in%c("renal cyst","non-functional kidney")]
pdata<-pData(x)

sapply(x$geo_accession,getGEOSuppFiles,makeDirectory=F,baseDir="~/CDC3/data/kidney_datasets/ncbi_geo/GSE68417")

gzfiles<-list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE68417",".CEL.gz",full.names=T)
sapply(gzfiles,gunzip)
cel<-list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE68417",".CEL",full.names=T)

raw<-read.celfiles(filenames = cel,
                   sampleNames = rownames(pdata))
dataset<-frma(raw,target="core")
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hugene10sttranscriptcluster.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file=gsub(".RData","_collMaxMean.RData",files[19]))

file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE68417",".CEL",full.names=T))
file.remove(list.files("~/CDC3/data/kidney_datasets/ncbi_geo/GSE68417","series_matrix.RData",full.names=T))

#***************************
# TCGA-KIRC
#***************************

dir.create("~/CDC3/data/kidney_datasets/tcga_kirc")

q1<-GDCquery(project = "TCGA-KIRC",
         data.category = "Transcriptome Profiling",
         data.type = "Gene Expression Quantification",
         workflow.type = "HTSeq - FPKM")
GDCdownload(q1,directory="~/CDC3/data/kidney_datasets/tcga_kirc")

expdat <- GDCprepare(query = q1,
                     save = F,
                     summarizedExperiment=T,
                     directory=path.expand("~/CDC3/data/kidney_datasets/tcga_kirc"))

exp<-assays(expdat)$'HTSeq - FPKM'
pdata<-as.data.frame(colData(expdat))

exp.coll<-aggregate(exp,by=list(Symbol=rowData(expdat)$external_gene_name),sum)
fdata<-data.frame(Symbol=exp.coll$Symbol, stringsAsFactors = F)
rownames(fdata)<-fdata$Symbol

exp.coll<-as.matrix(exp.coll[,-1])

kirc<-ExpressionSet(assayData=as.matrix(exp.coll),
                    phenoData=new("AnnotatedDataFrame",pdata),
                    featureData=new("AnnotatedDataFrame",fdata))

save(kirc,file="~/CDC3/data/kidney_datasets/tcga_kirc/fpkm_tcga_kirc_collSum.RData")

#***************************
# TCGA-KIRP
#***************************

dir.create("~/CDC3/data/kidney_datasets/tcga_kirp")

q1<-GDCquery(project = "TCGA-KIRP",
             data.category = "Transcriptome Profiling",
             data.type = "Gene Expression Quantification",
             workflow.type = "HTSeq - FPKM")
GDCdownload(q1,directory="~/CDC3/data/kidney_datasets/tcga_kirp")

expdat <- GDCprepare(query = q1,
                     save = F,
                     summarizedExperiment=T,
                     directory=path.expand("~/CDC3/data/kidney_datasets/tcga_kirp"))

exp<-assays(expdat)$'HTSeq - FPKM'
pdata<-as.data.frame(colData(expdat))

exp.coll<-aggregate(exp,by=list(Symbol=rowData(expdat)$external_gene_name),sum)
fdata<-data.frame(Symbol=exp.coll$Symbol, stringsAsFactors = F)
rownames(fdata)<-fdata$Symbol

exp.coll<-as.matrix(exp.coll[,-1])

kirp<-ExpressionSet(assayData=as.matrix(exp.coll),
                    phenoData=new("AnnotatedDataFrame",pdata),
                    featureData=new("AnnotatedDataFrame",fdata))

save(kirp,file="~/CDC3/data/kidney_datasets/tcga_kirp/fpkm_tcga_kirp_collSum.RData")

#***************************
# TCGA-KICH
#***************************

dir.create("~/CDC3/data/kidney_datasets/tcga_kich")

q1<-GDCquery(project = "TCGA-KICH",
             data.category = "Transcriptome Profiling",
             data.type = "Gene Expression Quantification",
             workflow.type = "HTSeq - FPKM")
GDCdownload(q1,directory="~/CDC3/data/kidney_datasets/tcga_kich")

expdat <- GDCprepare(query = q1,
                     save = F,
                     summarizedExperiment=T,
                     directory=path.expand("~/CDC3/data/kidney_datasets/tcga_kich"))

exp<-assays(expdat)$'HTSeq - FPKM'
pdata<-as.data.frame(colData(expdat))

exp.coll<-aggregate(exp,by=list(Symbol=rowData(expdat)$external_gene_name),sum)
fdata<-data.frame(Symbol=exp.coll$Symbol, stringsAsFactors = F)
rownames(fdata)<-fdata$Symbol

exp.coll<-as.matrix(exp.coll[,-1])

kich<-ExpressionSet(assayData=as.matrix(exp.coll),
                    phenoData=new("AnnotatedDataFrame",pdata),
                    featureData=new("AnnotatedDataFrame",fdata))

save(kich,file="~/CDC3/data/kidney_datasets/tcga_kich/fpkm_tcga_kich_collSum.RData")

# download curated clinical data
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", destfile = "~/CDC3/data/1-s2.0-S0092867418302290-mmc1.xlsx")
clindat<-as.data.frame(read_xlsx("~/CDC3/data/1-s2.0-S0092867418302290-mmc1.xlsx", sheet=1))
clindat<-clindat[,-1]
write.table(clindat,file="~/CDC3/data/TCGA-CDR-SupplementalTableS1.txt",sep="\t",row.names=F,quote=F)

#***************************
# Cheval et al. Nephron dataset
#***************************

download.file("https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0046876.s003", destfile="~/CDC3/data/journal.pone.0046876.s003.XLS")
cheval<-as.data.frame(read_xls("~/CDC3/data/journal.pone.0046876.s003.XLS", skip=1))
cheval<-cheval[,-2]
write.table(cheval, file="~/CDC3/data/cheval_et_al_nephron_dataset.txt",sep="\t",row.names=F,quote=F)

#***************************
# GSE131685
#***************************

dir.create("~/CDC3/data/kidney_datasets/GSE131685")

getGEOSuppFiles("GSE131685", makeDirectory = F, baseDir = "~/CDC3/data/kidney_datasets/GSE131685")
untar(tarfile=path.expand("~/CDC3/data/kidney_datasets/GSE131685/GSE131685_RAW.tar"),exdir = path.expand("~/CDC3/data/kidney_datasets/GSE131685"))
file.remove("~/CDC3/data/kidney_datasets/GSE131685/GSE131685_RAW.tar")

dir.create("~/CDC3/data/kidney_datasets/GSE131685/GSM4145204")
dir.create("~/CDC3/data/kidney_datasets/GSE131685/GSM4145205")
dir.create("~/CDC3/data/kidney_datasets/GSE131685/GSM4145206")

file.move(list.files("~/CDC3/data/kidney_datasets/GSE131685","kidney1",full.names=T),"~/CDC3/data/kidney_datasets/GSE131685/GSM4145204")
file.rename(from=list.files("~/CDC3/data/kidney_datasets/GSE131685/GSM4145204",full.names=T),to=gsub("GSM4145204_kidney1_","",list.files("~/CDC3/data/kidney_datasets/GSE131685/GSM4145204",full.names=T)))
file.move(list.files("~/CDC3/data/kidney_datasets/GSE131685","kidney2",full.names=T),"~/CDC3/data/kidney_datasets/GSE131685/GSM4145205")
file.rename(from=list.files("~/CDC3/data/kidney_datasets/GSE131685/GSM4145205",full.names=T),to=gsub("GSM4145205_kidney2_","",list.files("~/CDC3/data/kidney_datasets/GSE131685/GSM4145205",full.names=T)))
file.move(list.files("~/CDC3/data/kidney_datasets/GSE131685","kidney3",full.names=T),"~/CDC3/data/kidney_datasets/GSE131685/GSM4145206")
file.rename(from=list.files("~/CDC3/data/kidney_datasets/GSE131685/GSM4145206",full.names=T),to=gsub("GSM4145206_kidney3_","",list.files("~/CDC3/data/kidney_datasets/GSE131685/GSM4145206",full.names=T)))

K1.data <- Read10X(data.dir = "~/CDC3/data/kidney_datasets/GSE131685/GSM4145204")
K1 <- CreateSeuratObject(counts = K1.data, project = "kidney1", min.cells = 8, min.features = 200)
K2.data <- Read10X(data.dir = "~/CDC3/data/kidney_datasets/GSE131685/GSM4145205")
K2 <- CreateSeuratObject(counts = K2.data, project = "kidney2", min.cells = 6, min.features = 200)
K3.data <- Read10X(data.dir = "~/CDC3/data/kidney_datasets/GSE131685/GSM4145206")
K3 <- CreateSeuratObject(counts = K3.data, project = "kidney3", min.cells = 10, min.features = 200)
kid <- merge(x = K1, y = list(K2, K3))

kid[["percent.mt"]] <- PercentageFeatureSet(kid, pattern = "^MT-")
kid <- subset(kid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
kid <- NormalizeData(kid, normalization.method = "LogNormalize", scale.factor = 10000)
kid <- FindVariableFeatures(kid, selection.method = "vst", nfeatures = 2000)
kid <- CellCycleScoring(kid, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
kid <- ScaleData(kid, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(kid))

kid <- RunPCA(kid, pc.genes = kid@var.genes, npcs = 20, verbose = FALSE)
kid <- kid %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
kid <- kid %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.25) %>% 
  identity()
new.cluster.ids <- c(1,2,2,3,4,5,6,7,8,9,10)
names(new.cluster.ids) <- levels(kid)
kid <- RenameIdents(kid, new.cluster.ids)
save(kid,file="~/CDC3/data/kidney_datasets/GSE131685/GSE131685_processed.RData")

#*************************
# RMC from Msaouel et al
#*************************

# Fastq files for RMC tumors were downloaded from the sequence read archive (SRA) using the SRA toolkit (BioProject: PRJNA605003)
# Fastq file were aligned to the human genome version hg19 using the following command line:
# STAR --runThreadN 20 /--genomeDir <path to genomeDir> --readFilesIn <path to fastq1> <path to fastq1> --sjdbGTFfile <path to gtf file> --quantMode GeneCounts
# Genome indices were generated as described in the STAR manual.
# The human genome fasta file and GTF file were downloaded from Gencode version 32 (files GRCh37.primary_assembly.genome.fa, gencode.v32lift37.annotation.gtf)
# Save STAR output in the following directory ~/CDC3/data/kidney_datasets/PRJNA605003

dir.create("~/CDC3/data/kidney_datasets/PRJNA605003")

files<-list.files("~/CDC3/data/kidney_datasets/PRJNA605003","out.tab",full.names=T)
                                                                          
if (length(files) > 0){                                                                          
  samplename<-gsub("ReadsPerGene.out.tab","",gsub(".*\\/","",files))
  pdata<-data.frame(Run=samplename,histology=c("Normal","Normal","Tumor","Tumor","Tumor","Normal","Tumor","Tumor","Tumor","Tumor","Normal","Tumor","Normal","Normal","Tumor","Tumor","Tumor"))

  exp<-read.table(files[1],sep="\t",as.is=T)
  exp<-exp[,c(1,4)]
  rownames(exp)<-exp$V1
  for(i in 2:length(files)){
    temp<-read.table(files[i],sep="\t",as.is=T,row.names=1)
    exp<-cbind(exp,temp$V4)
  }
  exp<-exp[grep("ENSG",rownames(exp)),-1]
  colnames(exp)<-samplename

  rownames(pdata)<-pdata$Run
  pdata<-pdata[colnames(exp),]
  identical(rownames(pdata),colnames(exp))

  exp<-data.frame(Ensembl_ID=gsub("\\..*","",rownames(exp)),exp, stringsAsFactors = F)

  annot<-read.table("~/CDC3/data/gene_annotation_gencode32_hg19.txt",header=T,sep="\t",as.is=T)
  annot<-annot[,c(1,3:4)]
  annot$Ensembl_ID=gsub("\\..*","",annot$gene_id)
  annot<-unique(annot[,-1])
  exp<-merge(annot,exp,by="Ensembl_ID",all.y=T)
  fdata<-aggregate(exp$gene_length,by=list(Symbol=exp$gene_name),max)
  colnames(fdata)[2]<-"gene_length"
  rownames(fdata)<-fdata$Symbol
  exp<-aggregate(exp[,-c(1:3)],by=list(Symbol=exp$gene_name),sum)
  rownames(exp)<-exp$Symbol
  exp<-as.matrix(exp[,-1])

  dataset<-ExpressionSet(assayData=exp,
                       featureData = new("AnnotatedDataFrame",fdata),
                       phenoData=new("AnnotatedDataFrame",pdata))
  save(dataset,file="~/CDC3/data/kidney_datasets/PRJNA605003/PRJNA605003_raw_counts_collSum.RData")
  } else {
    print("No files detected. Skipping dataset PRJNA605003.")
  }

#*************************
# GSE2109
#*************************

dir.create("~/CDC3/data/GSE2109")

x<-getGEO("GSE2109",destdir = "~/CDC3/data/GSE2109")
x<-x[[1]]
ren<-pData(x)[grep("renal pelvis|ureter",x$source_name_ch1,ignore.case=T),]
ren<-ren[!ren$source_name_ch1%in%c("Kidney and Ureter","Ureterovesicle junction"),]
pdata<-ren[,c(1,2,8)]
colnames(pdata)[3]<-"source"

sapply(rownames(ren),getGEOSuppFiles,makeDirectory = F, baseDir = "~/CDC3/data/GSE2109")
cel<-list.celfiles("~/CDC3/data/GSE2109")

raw<-ReadAffy(filenames = cel,
              compress = T,
              celfile.path = path.expand("~/CDC3/data/GSE2109"))
dataset<-frma(raw)
colnames(dataset)<-gsub("\\..*","",colnames(dataset))
dataset<-dataset[,rownames(pdata)]
identical(rownames(pdata),colnames(dataset))
pData(dataset)<-pdata

dataset<-annotateEset(dataset,columns = c("PROBEID", "ENTREZID", "SYMBOL"),x="hgu133plus2.db")
dataset<-dataset[!is.na(fData(dataset)$SYMBOL),]

collrow<-collapseRows(exprs(dataset),
                      rowID=rownames(dataset),
                      rowGroup=fData(dataset)$SYMBOL,
                      method="MaxMean")
dataset.coll<-dataset[collrow$selectedRow,]
save(dataset.coll,file="~/CDC3/data/GSE2109/GSE2109_UTUC_frma_collMaxMean.RData")

file.remove(list.files("~/CDC3/data/GSE2109","CEL.gz",full.names=T))
file.remove("~/CDC3/data/GSE2109/GSE2109_series_matrix.txt.gz")
file.remove("~/CDC3/data/GSE2109/GPL570.soft")

#*************************
# UTUC from cBioPortal
#*************************

dir.create("~/CDC3/data/utuc_cbm")

download.file("http://download.cbioportal.org/utuc_cornell_baylor_mdacc_2019.tar.gz",destfile = "~/CDC3/data/utuc_cbm/utuc_cornell_baylor_mdacc_2019.tar.gz")
gunzip(path.expand("~/CDC3/data/utuc_cbm/utuc_cornell_baylor_mdacc_2019.tar.gz"))
untar(tarfile=path.expand("~/CDC3/data/utuc_cbm/utuc_cornell_baylor_mdacc_2019.tar"), exdir=path.expand("~/CDC3/data/utuc_cbm"))

exp<-read.table("~/CDC3/data/utuc_cbm/data_RNA_Seq_expression_median.txt",header=T,sep="\t",as.is=T,row.names=1)
exp<-as.matrix(exp[,-1])
pdata<-read.table("~/CDC3/data/utuc_cbm/data_clinical_patient.txt",header=T,sep="\t",as.is=T)
rownames(pdata)<-gsub("-","\\.",pdata$PATIENT_ID)
pdata<-pdata[colnames(exp),]
identical(rownames(pdata),colnames(exp))
fdata<-data.frame(Symbol=rownames(exp),stringsAsFactors = F)
rownames(fdata)<-rownames(exp)
dataset<-ExpressionSet(assayData = exp,
                       phenoData=new("AnnotatedDataFrame",pdata),
                       featureData = new("AnnotatedDataFrame",fdata))
save(dataset,file="~/CDC3/data/utuc_cbm/data_RNA_Seq_mRNA_median.RData")

file.remove("~/CDC3/data/utuc_cbm/utuc_cornell_baylor_mdacc_2019.tar")
file.remove(list.files("~/CDC3/data/utuc_cbm","txt",full.names=T))
file.remove("~/CDC3/data/utuc_cbm/LICENSE")
unlink("~/CDC3/data/utuc_cbm/case_lists", recursive=T)

#*************************
# ICGC-RECA-EU
#*************************

# RNA-Seq and sample data for ICGC RECA-EU project were donloaded from the following link: https://dcc.icgc.org/releases/release_28/Projects/RECA-EU
# Put downloaded files in directory ~/CDC3/data/kidney_datasets/icgc_reca_eu

dir.create("~/CDC3/data/kidney_datasets/icgc_reca_eu")

download.file("https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/RECA-EU/exp_seq.RECA-EU.tsv.gz", destfile = "~/CDC3/data/kidney_datasets/icgc_reca_eu/exp_seq.RECA-EU.tsv.gz")
download.file("https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/RECA-EU/sample.RECA-EU.tsv.gz", destfile = "~/CDC3/data/kidney_datasets/icgc_reca_eu/sample.RECA-EU.tsv.gz")
download.file("https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/RECA-EU/donor.RECA-EU.tsv.gz", destfile = "~/CDC3/data/kidney_datasets/icgc_reca_eu/donor.RECA-EU.tsv.gz")
download.file("https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/RECA-EU/specimen.RECA-EU.tsv.gz", destfile = "~/CDC3/data/kidney_datasets/icgc_reca_eu/specimen.RECA-EU.tsv.gz")
files<-list.files("~/CDC3/data/kidney_datasets/icgc_reca_eu",full.names=T,"tsv.gz")
sapply(files,gunzip)

exp<-read.table("~/CDC3/data/kidney_datasets/icgc_reca_eu/exp_seq.RECA-EU.tsv", sep="\t", header=T,as.is=T,comment.char = "")
exp2<-dcast(exp, exp$gene_id ~ exp$submitted_sample_id, value.var="raw_read_count")
colnames(exp2)[1]<-"Ensembl_ID"

info_donor<- read.table("~/CDC3/data/kidney_datasets/icgc_reca_eu/donor.RECA-EU.tsv", sep="\t", header=T,as.is=T,comment.char = "")
sample<- read.table("~/CDC3/data/kidney_datasets/icgc_reca_eu/sample.RECA-EU.tsv", sep="\t", header=T,as.is=T,comment.char = "")
specimen<- read.table("~/CDC3/data/kidney_datasets/icgc_reca_eu/specimen.RECA-EU.tsv", sep="\t", header=T,as.is=T,comment.char = "")

annot<-read.table("~/CDC3/data/gene_annotation_gencode32_hg19.txt",header=T,sep="\t",as.is=T)
annot<-annot[,c(1,3:4)]
annot$Ensembl_ID<-gsub("\\..*","",annot$gene_id)
annot<-unique(annot[,-1])
exp3<-merge(annot,exp2,by="Ensembl_ID",all.y=T)
fdata<-exp3[,2:3]
fdata<-aggregate(fdata$gene_length,by=list(Symbol=fdata$gene_name),max)
colnames(fdata)[2]<-"gene_length"
rownames(fdata)<-fdata$Symbol
exp3<-aggregate(exp3[,-c(1:3)],by=list(Symbol=exp3$gene_name),sum)
rownames(exp3)<-exp3$Symbol
exp3<-as.matrix(exp3[,-1])
identical(rownames(fdata),rownames(exp3))

pdata<-specimen[match(colnames(exp3),specimen$submitted_specimen_id),c(4,1,5:7,21,22,25,29)]
rownames(pdata)<-pdata$submitted_specimen_id
identical(colnames(exp3),rownames(pdata))

dataset<-ExpressionSet(assayData=exp3,
                       phenoData=new("AnnotatedDataFrame",pdata),
                       featureData=new("AnnotatedDataFrame",fdata))

dataset$specimen_type<-ifelse(dataset$specimen_type=="Normal - tissue adjacent to primary", "Normal_Adjacent_Tissue","Primary_solid_tumor")

save(dataset,file="~/CDC3/data/kidney_datasets/icgc_reca_eu/icgc_reca_eu_raw_counts.RData")

file.remove("~/CDC3/data/kidney_datasets/icgc_reca_eu/exp_seq.RECA-EU.tsv")
file.remove("~/CDC3/data/kidney_datasets/icgc_reca_eu/sample.RECA-EU.tsv")
file.remove("~/CDC3/data/kidney_datasets/icgc_reca_eu/donor.RECA-EU.tsv")
file.remove("~/CDC3/data/kidney_datasets/icgc_reca_eu/specimen.RECA-EU.tsv")

#**********************************************************
# Assembly of a dataset with all kidney cancer histologies
#**********************************************************

# import TCGA datasets
load("~/CDC3/data/kidney_datasets/tcga_kirc/fpkm_tcga_kirc_collSum.RData")
pdata<-data.frame(sample=colnames(kirc),histology=ifelse(kirc$definition=="Solid Tissue Normal","Normal","ccRCC"),stringsAsFactors = F)
exprs(kirc)<-log2(exprs(kirc)+1)
load("~/CDC3/data/kidney_datasets/tcga_kirp/fpkm_tcga_kirp_collSum.RData")
temp<-data.frame(sample=colnames(kirp),histology=ifelse(kirp$definition=="Solid Tissue Normal","Normal","papillary"),stringsAsFactors = F)
exprs(kirp)<-log2(exprs(kirp)+1)
pdata<-rbind(pdata,temp)
load("~/CDC3/data/kidney_datasets/tcga_kich/fpkm_tcga_kich_collSum.RData")
temp<-data.frame(sample=colnames(kich),histology=ifelse(kich$definition=="Solid Tissue Normal","Normal","chromophobe"),stringsAsFactors = F)
exprs(kich)<-log2(exprs(kich)+1)
pdata<-rbind(pdata,temp)

# import RMC
load("~/CDC3/data/kidney_datasets/PRJNA605003/PRJNA605003_raw_counts_collSum.RData")
rmc<-dataset
temp<-data.frame(sample=colnames(rmc),histology=ifelse(rmc$histology=="Tumor","RMC","Normal"),stringsAsFactors = F)
pdata<-rbind(pdata,temp)
rmc<-rmc[!is.na(fData(rmc)$gene_length),]
colnames(fData(rmc))[2]<-"basepairs"
dds <- DESeqDataSetFromMatrix(countData = exprs(rmc),colData = pData(rmc),design=~1)
mcols(dds)<-DataFrame(mcols(dds),fData(rmc))
rmc<-fpkm(dds)
rmc<-log2(rmc+1)

# import microarray datasets
files<-list.files("~/CDC3/data/kidney_datasets/ncbi_geo","collMaxMean.RData",full.names=T,recursive=T)
arrays<-list()
genes.arrays<-list()
annot<-NULL
for(i in files){
  print(paste("loading dataset",gsub(".*\\/","",i)))
  load(i)
  annot<-c(annot,dataset.coll@annotation)
  if(max(exprs(dataset.coll))>50){
    exprs(dataset.coll)<-log2(exprs(dataset.coll))
  }
  dataset.coll<-dataset.coll[,grep("normal|clear-cell|conventional|clear cell|ccRCC|pap|chromo|oncocy|collecting",dataset.coll$tumor_type,ignore.case=T)]
  if(length(grep("sarcoma",dataset.coll$tumor_type))>0){
    dataset.coll<-dataset.coll[,-grep("sarcoma",dataset.coll$tumor_type,ignore.case=T)]
  }
  colnames(fData(dataset.coll))[grep("SYMBOL",colnames(fData(dataset.coll)),ignore.case = T)]<-"Symbol"
  arrays<-c(arrays,dataset.coll)
  genes.arrays<-c(genes.arrays,list(as.vector(fData(dataset.coll)$Symbol)))
}
cg.arrays<-Reduce(intersect,genes.arrays)

histology.array<-pData(arrays[[1]])[,colnames(pData(arrays[[1]]))%in%c("geo_accession","tumor_type")]
for(i in 2:length(arrays)){
  temp<-pData(arrays[[i]])[,colnames(pData(arrays[[i]]))%in%c("geo_accession","tumor_type")]
  histology.array<-rbind(histology.array,temp)
}
colnames(histology.array)<-c("sample","histology")
pdata<-rbind(pdata,histology.array)

# import UTUC dataset GSE2109
load("~/CDC3/data/GSE2109/GSE2109_UTUC_frma_collMaxMean.RData")
utuc<-dataset.coll
rownames(utuc)<-fData(utuc)$SYMBOL
histology.utuc<-data.frame(sample=colnames(utuc),histology="UTUC")
pdata<-rbind(pdata,histology.utuc)

# import UTUC dataset cbioportal
load("~/CDC3/data/utuc_cbm/data_RNA_Seq_mRNA_median.RData")
utuc2<-dataset
exprs(utuc2)<-log2(exprs(utuc2)+1)
histology.utuc2<-data.frame(sample=colnames(utuc2),histology="UTUC")
pdata<-rbind(pdata,histology.utuc2)
utuc2<-utuc2[rowSums(is.na(exprs(utuc2)))==0,]

# INT
load("~/CDC3/data/kidney_datasets/int/INT_processed.RData")
int<-dataset[fData(dataset)$SYMBOL!="NA",]
rownames(int)<-fData(int)$SYMBOL
temp<-data.frame(sample=colnames(int),histology=int$histology,stringsAsFactors = F)
pdata<-rbind(pdata,temp)

# GSE89122
load("~/CDC3/data/kidney_datasets/GSE89122/GSE89122_raw_counts.RData")
wang<-dataset[!is.na(fData(dataset)$gene_length),]
temp<-data.frame(sample=colnames(wang),histology=wang$histology,stringsAsFactors = F)
pdata<-rbind(pdata,temp)
colnames(fData(wang))[2]<-"basepairs"
dds <- DESeqDataSetFromMatrix(countData = exprs(wang),colData = pData(wang),design=~1)
mcols(dds)<-DataFrame(mcols(dds),fData(wang))
wang<-fpkm(dds)
wang<-log2(wang+1)

# WACH
load("~/CDC3/data/kidney_datasets/wach/wach_raw_counts.RData")
wach<-dataset
temp<-data.frame(sample=colnames(wach),histology=wach$histology,stringsAsFactors = F)
pdata<-rbind(pdata,temp)
rownames(wach)<-fData(wach)$Symbol
wach<-wach[!is.na(fData(wach)$gene_length),]
colnames(fData(wach))[2]<-"basepairs"
dds <- DESeqDataSetFromMatrix(countData = exprs(wach),colData = pData(wach),design=~1)
mcols(dds)<-DataFrame(mcols(dds),fData(wach))
wach<-fpkm(dds)
wach<-log2(wach+1)

# ICGC
load("~/CDC3/data/kidney_datasets/icgc_reca_eu/icgc_reca_eu_raw_counts.RData")
icgc<-dataset
temp<-data.frame(sample=colnames(icgc),histology=ifelse(icgc$specimen_type=="Primary_solid_tumor","ccRCC","Normal"),stringsAsFactors = F)
pdata<-rbind(pdata,temp)
rownames(icgc)<-fData(icgc)$Symbol
icgc<-icgc[!is.na(fData(icgc)$gene_length),]
colnames(fData(icgc))[2]<-"basepairs"
dds <- DESeqDataSetFromMatrix(countData = exprs(icgc),colData = pData(icgc),design=~1)
mcols(dds)<-DataFrame(mcols(dds),fData(icgc))
icgc<-fpkm(dds)
icgc<-log2(icgc+1)

all.platforms<-list(rownames(kirc),rownames(kirp),rownames(kich),cg.arrays,rownames(utuc),rownames(utuc2),rownames(int),rownames(wang),rownames(wach),rownames(icgc),rownames(rmc))
cg.all<-Reduce(intersect,all.platforms)

utuc<-utuc[cg.all,]
utuc2<-utuc2[cg.all,]
int<-int[cg.all,]
wang<-wang[cg.all,]
wach<-wach[cg.all,]
icgc<-icgc[cg.all,]
kirc<-kirc[cg.all,]
kich<-kich[cg.all,]
kirp<-kirp[cg.all,]
rmc<-rmc[cg.all,]
for(i in 1:length(arrays)){
  rownames(arrays[[i]])<-fData(arrays[[i]])$Symbol
  arrays[[i]]<-arrays[[i]][cg.all,]
}

datasets<-list(int,wang,wach,icgc,kirc,kirp,kich,rmc)
datasets<-c(datasets,arrays,utuc,utuc2)
names(datasets)<-c("INT","GSE89122","WACH","ICGC","TCGA-KIRC","TCGA-KIRP","TCGA-KICH","RMC-PRJNA605003",gsub("_.*","",gsub(".*\\/","",files)),"GSE2109-UTUC","UTUC-CORNELL-BAYLOR-MDACC")

# create expression set with all datasets

exp<-exprs(datasets[[1]])
batch=rep(names(datasets)[1],ncol(datasets[[1]]))
samplename<-colnames(datasets[[1]])
for(i in 2:length(datasets)){
  if(class(datasets[[i]])=="matrix"){
    exp<-cbind(exp,datasets[[i]])
  } else {
    exp<-cbind(exp,exprs(datasets[[i]]))
  }
  batch<-c(batch,rep(names(datasets)[i],ncol(datasets[[i]])))
  samplename<-c(samplename,colnames(datasets[[i]]))
}
pdata2<-data.frame(sample=samplename,dataset=batch,stringsAsFactors = F)
pdata<-merge(pdata,pdata2,by="sample")
pdata<-unique(pdata)
rownames(pdata)<-pdata$sample
pdata<-pdata[colnames(exp),]
unscaled<-ExpressionSet(assayData=exp,
                        phenoData=new("AnnotatedDataFrame",pdata))
table(unscaled$histology)
unscaled$histology[grep("normal|SANO",unscaled$histology,ignore.case=T)]<-"normal"
unscaled$histology[grep("pap",unscaled$histology,ignore.case=T)]<-"papillary"
unscaled$histology[grep("chromo",unscaled$histology,ignore.case=T)]<-"chromophobe"
unscaled$histology[grep("clear|CC|conventional",unscaled$histology,ignore.case=T)]<-"ccRCC"
unscaled$histology[grep("collecting",unscaled$histology,ignore.case=T)]<-"CDC"
unscaled$histology[grep("oncocy",unscaled$histology,ignore.case=T)]<-"oncocytoma"
save(unscaled,file="~/CDC3/data/kidney_datasets/all_kidney_histologies_expression_dataset.RData")

#**********************************************************
# Assembly of a dataset with CDC and normal samples
#**********************************************************

load("~/CDC3/data/kidney_datasets/int/INT_processed.RData")
int<-dataset
load("~/CDC3/data/kidney_datasets/ncbi_geo/GSE11151/GSE11151_series_matrix_collMaxMean.RData")
gse<-dataset.coll
load("~/CDC3/data/kidney_datasets/GSE89122/GSE89122_raw_counts.RData")
wang<-dataset
load("~/CDC3/data/kidney_datasets/wach/wach_raw_counts.RData")
wach<-dataset[rowSds(exprs(dataset))>0,]

cg<-intersect(intersect(fData(int)$SYMBOL,fData(gse)$SYMBOL),intersect(rownames(wang),rownames(wach)))

# INT
int<-int[,int$histology!="ccRCC"]
rownames(int)<-fData(int)$SYMBOL
int<-int[cg,]

# GSE11151
gse<-gse[,gse$tumor_type%in%c("collecting duct carcinoma","adult normal kidney")]
rownames(gse)<-fData(gse)$SYMBOL
gse<-gse[cg,]

# WANG
wang<-wang[cg,]
dge<-DGEList(counts=exprs(wang), genes=fData(wang),samples=pData(wang))
dge <- calcNormFactors(dge,method="TMM")
wang.cpm<-cpm(dge,normalized.lib.sizes = T,log = T,prior.count = 1)
exprs(wang)<-wang.cpm

# WACH
wach<-wach[cg,]
dge<-DGEList(counts=exprs(wach), genes=fData(wach),samples=pData(wach))
dge <- calcNormFactors(dge,method="TMM")
wach.cpm<-cpm(dge,normalized.lib.sizes = T,log = T,prior.count = 1)
exprs(wach)<-wach.cpm

exp.zscore<-cbind(t(scale(t(exprs(int)))),t(scale(t(exprs(gse)))),t(scale(t(exprs(wang)))),t(scale(t(exprs(wach)))))

#exp<-cbind(exprs(int),exprs(gse),exprs(wang),exprs(wach))
pdata<-data.frame(sampleID=colnames(exp.zscore),batch=c(rep("INT",ncol(int)),rep("GSE11151",ncol(gse)),rep("GSE89122",ncol(wang)),rep("WACH",ncol(wach))),histology=c(int$histology,gse$tumor_type,wang$histology,as.vector(wach$histology)),stringsAsFactors = F)
pdata$histology[grep("normal|sano",pdata$histology,ignore.case = T)]<-"Normal"
pdata$histology[pdata$histology!="Normal"]<-"CDC"
rownames(pdata)<-pdata$sampleID
fdata<-data.frame(Symbol=rownames(exp.zscore))
rownames(fdata)<-fdata$Symbol


dataset.zscore<-ExpressionSet(assayData=exp.zscore,
                              phenoData=new("AnnotatedDataFrame",pdata),
                              featureData = new("AnnotatedDataFrame",fdata))
dataset.zscore<-dataset.zscore[,dataset.zscore$histology=="CDC"]
exprs(dataset.zscore)<-normalizeBetweenArrays(exprs(dataset.zscore),method="quantile")
save(dataset.zscore,file="~/CDC3/data/kidney_datasets/metadataset_CDC_zscore_quantile.RData")

#**********************************************
# CCLE and CTRP data
#**********************************************

dir.create("~/CDC3/data/ccle_ctrp")

# Download RNA-Seq FPKM data (file CCLE_RNAseq_genes_rpkm_20180929.gct.gz) cell line annotation (file Cell_lines_annotations_20181226.txt)
# and gene annotations (file gencode.v19.genes.v7_model.patched_contigs.gtf.gz) 
# from https://portals.broadinstitute.org/ccle/data (login required).
# Save downloaded file in folder ~/CDC3/data/ccle_ctrp
                                                                          
isPresent<-sum(list.files("~/CDC3/data/ccle_ctrp") %in% c("CCLE_RNAseq_genes_rpkm_20180929.gct.gz", "Cell_lines_annotations_20181226.txt", "gencode.v19.genes.v7_model.patched_contigs.gtf.gz"))

if(isPresent == 3){
  gunzip(filename=path.expand("~/CDC3/data/ccle_ctrp/CCLE_RNAseq_genes_rpkm_20180929.gct.gz"),destname=path.expand("~/CDC3/data/ccle_ctrp/CCLE_RNAseq_genes_rpkm_20180929.gct"))
  gunzip(filename=path.expand("~/CDC3/data/ccle_ctrp/gencode.v19.genes.v7_model.patched_contigs.gtf.gz"),destname=path.expand("~/CDC3/data/ccle_ctrp/gencode.v19.genes.v7_model.patched_contigs.gtf"))

  exp<-read.gct("~/CDC3/data/ccle_ctrp/CCLE_RNAseq_genes_rpkm_20180929.gct")
  colnames(exp)<-gsub("^X","",colnames(exp))

  pdata<-read.table("~/CDC3/data/ccle_ctrp/Cell_lines_annotations_20181226.txt",header=T,sep="\t",as.is=T,quote="",comment.char="")
  rownames(pdata)<-pdata$CCLE_ID
  pdata<-pdata[colnames(exp),]
  identical(rownames(pdata),colnames(exp))

  fdata<-read.table("~/CDC3/data/ccle_ctrp/gencode.v19.genes.v7_model.patched_contigs.gtf",header=F,sep="\t",as.is=T)

  fdata$gene_id<-gsub(".* ","",sapply(strsplit(fdata$V9,";"),"[",1))
  fdata$gene_type<-gsub(".* ","",sapply(strsplit(fdata$V9,";"),"[",3))
  fdata$gene_status<-gsub(".* ","",sapply(strsplit(fdata$V9,";"),"[",4))
  fdata$gene_name<-gsub(".* ","",sapply(strsplit(fdata$V9,";"),"[",5))
  fdata2<-unique(fdata[,c(2,10:13)])
  rownames(fdata2)<-fdata2$gene_id
  common<-intersect(rownames(fdata2),rownames(exp))
  fdata2<-fdata2[common,]
  exp<-exp[common,]
  identical(rownames(fdata2),rownames(exp))

  exp.coll<-aggregate(exp,by=list(Symbol=fdata2$gene_name),sum)
  rownames(exp.coll)<-exp.coll$Symbol
  exp.coll<-as.matrix(exp.coll[,-1])
  colnames(exp.coll)<-gsub("^X","",colnames(exp.coll))
  identical(rownames(pdata),colnames(exp.coll))
  fdata<-data.frame(Symbol=rownames(exp.coll),stringsAsFactors = F)
  rownames(fdata)<-rownames(exp.coll)
  dataset<-ExpressionSet(assayData = exp.coll,
                       phenoData=new("AnnotatedDataFrame",pdata),
                       featureData = new("AnnotatedDataFrame",fdata))
  save(dataset,file="~/CDC3/data/ccle_ctrp/CCLE_RNAseq_genes_rpkm_collSum.RData")
  } else {
    print("Missing files for CCLE. Skipping it.")
  }

# Drug response data in the Cancer Therapeutic Response Portal version 2 (CTRP) dataset (Seashore-Ludlow et al., 2015)
# was downloaded from the supplementary information files of the corresponding main publication (PUBMED ID: 26482930)

download.file("https://cancerdiscovery.aacrjournals.org/highwire/filestream/36278/field_highwire_adjunct_files/1/145780_2_supp_3058746_nrhtdz.xlsx", destfile="~/CDC3/data/ccle_ctrp/145780_2_supp_3058746_nrhtdz.xlsx")

# drug annotation
df<-as.data.frame(read_xlsx("~/CDC3/data/ccle_ctrp/145780_2_supp_3058746_nrhtdz.xlsx",sheet=2))
write.table(df,file="~/CDC3/data/ccle_ctrp/ctrp_compounds.txt",sep="\t",row.names=F,quote=F)

# cell line annotation
df<-as.data.frame(read_xlsx("~/CDC3/data/ccle_ctrp/145780_2_supp_3058746_nrhtdz.xlsx",sheet=3))
write.table(df,file="~/CDC3/data/ccle_ctrp/ctrp_cell_lines.txt",sep="\t",row.names=F,quote=F)

# auc values
df<-as.data.frame(read_xlsx("~/CDC3/data/ccle_ctrp/145780_2_supp_3058746_nrhtdz.xlsx",sheet=4))
write.table(df,file="~/CDC3/data/ccle_ctrp/ctrp_drugs_auc.txt",sep="\t",row.names=F,quote=F)

file.remove(c("~/CDC3/data/ccle_ctrp/CCLE_RNAseq_genes_rpkm_20180929.gct",
              "~/CDC3/data/ccle_ctrp/Cell_lines_annotations_20181226.txt",
              "~/CDC3/data/ccle_ctrp/gencode.v19.genes.v7_model.patched_contigs.gtf",
              "~/CDC3/data/ccle_ctrp/145780_2_supp_3058746_nrhtdz.xlsx"))

#***************************************
# GDSC
#***************************************

dir.create("~/CDC3/data/gdsc")

# download cell line annotation, screened compound annotation and drugs auc values fro supplementary files of Iorio et al. 2016.

download.file("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx",destfile="~/CDC3/data/gdsc/TableS1E_cell_lines_details.xlsx")
download.file("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1F.xlsx",destfile="~/CDC3/data/gdsc/TableS1F_screened_compounds.xlsx")
sc<-as.data.frame(read_xlsx("~/CDC3/data/gdsc/TableS1F_screened_compounds.xlsx",skip=2))
colnames(sc)<-gsub(" ","_",colnames(sc))
colnames(sc)[8]<-"Targeted_pathway"
write.table(sc,file="~/CDC3/data/gdsc/TableS1F_screened_compounds.txt",sep="\t",row.names=F,quote=F)
download.file("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4B.xlsx",destfile="~/CDC3/data/gdsc/TableS4B_cell_lines_drug_auc.xlsx")
aucval<-as.data.frame(read_xlsx("~/CDC3/data/gdsc/TableS4B_cell_lines_drug_auc.xlsx",skip=4))
aucval<-aucval[-1,]
write.table(aucval,file="~/CDC3/data/gdsc/TableS4B_cell_lines_drug_auc.txt",sep="\t",row.names=F,quote=F)

# download gene expression data from ArrayExpress

rawset = ArrayExpress("E-MTAB-3610", path=path.expand("~/CDC3/data/gdsc"))

pdata<-pData(rawset)[,1:3]
pdata$Array.Data.File<-rownames(pdata)
colnames(pdata)[2:3]<-c("organism","cell_line")
info<-as.data.frame(read_xlsx("~/CDC3/data/gdsc/TableS1E_cell_lines_details.xlsx",skip=2))
info<-info[-1,-c(3:7)]
colnames(info)[1]<-"cell_line"
pdata<-merge(pdata,info,by="cell_line",all.x=T)
rownames(pdata)<-pdata$Array.Data.File

dataset<-rma(rawset)
pdata<-pdata[colnames(exprs(dataset)),]
identical(rownames(pdata),colnames(exprs(dataset)))

Entrez<-unlist(mget(rownames(exprs(dataset)),hgu219ENTREZID,ifnotfound=NA))
Symbol<-unlist(mget(rownames(exprs(dataset)),hgu219SYMBOL,ifnotfound=NA))

fdata<-data.frame(ProbeSet=rownames(exprs(dataset)),Entrez=Entrez,Symbol=Symbol,stringsAsFactors=F)
rownames(fdata)<-rownames(exprs(dataset))
fData(dataset)<-fdata
pData(dataset)<-pdata

# collapse
dataset<-dataset[!is.na(fData(dataset)$Symbol),]
x<-collapseRows(exprs(dataset),
                rowID=rownames(exprs(dataset)),
                rowGroup=fData(dataset)$Symbol,
                method="maxRowVariance"
)

dataset.coll<-dataset[x$selectedRow,]
save(dataset.coll,file="~/CDC3/data/gdsc/E_MTAB_3610_collMaxVar.RData")
