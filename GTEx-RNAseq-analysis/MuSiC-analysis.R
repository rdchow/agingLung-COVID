#!/usr/bin/Rscript
library(MuSiC)
library(Biobase)
library(xbioc)
library(Seurat)
data<-readRDS("hlca-merged.2.rds") # Combined Seurat RDS from the Human Lung Cell Atlas
counts<-data[['RNA']]@counts
counts<-as.matrix(counts)
metadata<-data@meta.data
phenoData<-new("AnnotatedDataFrame",data=metadata)
sc.eset<-ExpressionSet(assayData=counts,phenoData=phenoData)

library(data.table)
bulk<-data.frame(fread("GTEX-counts.Lung.geneSum.txt",sep="\t"),row.names=1)
bulkm<-as.matrix(bulk)
bulk.eset<-ExpressionSet(assayData=bulkm)


results<-music_prop(bulk.eset=bulk.eset,sc.eset=sc.eset,clusters="free_annotation",samples="orig.ident",verbose=TRUE)
saveRDS(results,"GTEx-music-analysis.rds")
