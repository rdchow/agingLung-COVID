#!/usr/bin/Rscript
library(Seurat)
data<-readRDS("GSE145926-merged.all.rds")
meta<-data[[]]

patient <- rownames(meta)
patient2 <- vapply(strsplit(patient,"_"), `[`, 1, FUN.VALUE=character(1))


data2<-AddMetaData(data,metadata=patient2,col.name="patient")

library(scater)
data.sce <- as.SingleCellExperiment(data2)
pseudoMat<-sumCountsAcrossCells(data.sce,ids=data.sce$patient)

write.table(as.matrix(pseudoMat),"GSE145926-pseudoBulk.counts.txt",sep="\t",row.names=TRUE)