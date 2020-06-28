#!/usr/bin/Rscript
library(data.table)

#####
data = data.frame(fread("GSE145926-pseudoBulk.counts.txt"))
info<-read.table("GSE145926-sample-info.txt",sep="\t",header=TRUE)

infof2<-info
dataf<-data
all(infof2[,1] == colnames(dataf)[2:ncol(dataf)])

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=dataf, colData=infof2, design=~Type, tidy = TRUE)
dds <- DESeq(dds)

resultsNames(dds)

resSev<-results(dds,tidy=TRUE,name="Type_severe_vs_ctrl")
resSev<-resSev[order(resSev$padj),]
head(resSev)
summary(resSev)
write.table(resSev,"GSE145926-severe_COVID-vs-Control.de.txt",sep="\t",row.names=TRUE)
