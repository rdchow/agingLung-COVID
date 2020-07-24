#!/usr/bin/Rscript
library(Biobase)
library(xbioc)
library(SingleR)
library(scater)
library(Seurat)

data<-readRDS("hlca-merged.2.rds")
counts<-data[['RNA']]@counts
counts<-as.matrix(counts)
metadata<-data@meta.data

pseudoMat<-sumCountsAcrossCells(counts,ids=metadata$free_annotation)
head(pseudoMat)

hpca.se = BlueprintEncodeData()
hpca.se


pred.sc = SingleR(test=pseudoMat,ref=hpca.se,labels=hpca.se$label.fine)
write.table(pred.sc,"hlca-pseudobulk-SingleR.out.txt",sep="\t",row.names=TRUE)
