#!/usr/bin/Rscript
library(data.table)

#####
data = data.frame(fread("GSE147507_RawReadCounts_Human.tsv"))
info<-read.table("GSE147507-sample-info.txt",sep="\t",header=TRUE)
infof<-info[info[,1] %in% colnames(data),]

all(infof[,1] == colnames(data)[2:ncol(data)])

infof2<-infof[infof$SARS2involved == "yes",]
dataf<-cbind(data[,1],data[,colnames(data) %in% infof2[,1]])
colnames(dataf)[1]<-"Name"

all(infof2[,1] == colnames(dataf)[2:ncol(dataf)])

groups<-unique(sort(infof2$ComparisonGroup))
library(DESeq2)

#loop through each comparison group, compare the effect of virus treatment
for (c in groups){
    myinfo = infof2[infof2$ComparisonGroup == c,]
    mydata = dataf[,c(TRUE,infof2$ComparisonGroup == c)]
    dds <- DESeqDataSetFromMatrix(countData=mydata, colData=myinfo, design=~Treatment, tidy = TRUE)
    dds <- DESeq(dds)
    res<-results(dds,tidy=TRUE)
    res<-res[order(res$padj),]
    head(res)
    summary(res)
    write.table(res,paste(c,"-SARS2-vs-Ctrl.DE.wald.txt",sep=""),sep="\t",row.names=TRUE)
}