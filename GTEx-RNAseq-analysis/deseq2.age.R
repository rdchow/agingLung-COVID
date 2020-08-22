#!/usr/bin/Rscript
library(stringr)
library(data.table)

#####

#####
data = data.frame(fread("GTEX-counts.filterTissue.geneSum.txt"))
info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1) ## note that the detailed clinical annotations are controlled access, so this code will not work with the provided sample table. An anonymized/randomized example of the full sample table is provided: "sample-info.filterTissue.anonymized.full.txt".
dataf<-cbind(data[,1], data[,colnames(data) %in% rownames(info)])

info2<-info[match(colnames(dataf)[2:ncol(dataf)],rownames(info)),]
all(rownames(info2) == colnames(dataf)[2:ncol(dataf)])

colnames(dataf)[1]<-"Name"
table(cbind(info2[4:5]))

infoR = info2[!is.na(info2$Smoking),]
infoR = infoR[!is.na(infoR$Hardy),] 

infoR$AgeBin<-str_replace(infoR$AgeBin,"-","_")
table(cbind(infoR[4:5]))

datafR<-cbind(dataf[,1],dataf[,colnames(dataf) %in% rownames(infoR)])
all(rownames(infoR) == colnames(datafR)[2:ncol(datafR)])

library(DESeq2)

infoR$Hardy = as.factor(infoR$Hardy)
table(infoR$Hardy)
str(infoR)

dds <- DESeqDataSetFromMatrix(countData=datafR, colData=infoR, design=~Sex+Smoking+Hardy+AgeBin, tidy = TRUE) # note that these other variables are controlled acess, so not present in the provided sample info file
dds = DESeq(dds,test="LRT",reduced= ~Sex+Smoking+Hardy)
saveRDS(dds,"DESeq2.LRT.agebin.rds")
resultsNames(dds)


