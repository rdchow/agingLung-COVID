#!/usr/bin/Rscript
library(data.table)

data = data.frame(fread("GTEX-counts.filterTissue.geneSum.txt"))
info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE)
infof<-info[info[,1] %in% colnames(data),]

all(infof[,1] == colnames(data)[2:ncol(data)])

myTissue<-"Lung"
infof2<-infof[infof$simpleTissue == myTissue,]
dataf<-cbind(data[,1],data[,colnames(data) %in% infof2[,1]])
colnames(dataf)[1]<-"Name"

all(infof2[,1] == colnames(dataf)[2:ncol(dataf)])
table(cbind(infof2[5:6]))


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=dataf, colData=infof2, design=~Age, tidy = TRUE)
dds <- DESeq(dds,test="LRT",reduced= ~1)

res<-results(dds,tidy=TRUE)
res<-res[order(res$padj),]
head(res)
summary(res)
write.table(res,"GTEx-lung-DEseq2-results.LRT.age.txt",sep="\t",row.names=TRUE)
