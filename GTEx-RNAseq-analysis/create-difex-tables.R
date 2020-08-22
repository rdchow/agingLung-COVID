#!/usr/bin/Rscript
library(DESeq2)
library(limma)
library(stringr)
library(data.table)

dds = readRDS("DESeq2.LRT.agebin.rds")
resultsNames(dds)
res = results(dds,contrast=c("AgeBin","70_79","60_69"),independentFiltering=TRUE,alpha=0.05,pAdjustMethod="BH") # change the ages to get all the other pairwise comparisons as desired
res<-res[order(res$padj),]
head(res)
summary(res)
write.table(res,"GTEx-lung-DEseq2-results.LRT.timeBin.7-6.SexSmokingIschemia.txt",sep="\t",row.names=TRUE) # change the file name for other pairwise comparisons


# VST and "batch correction" for downstream visualization of adjusted expression values
vsd = varianceStabilizingTransformation(dds,blind=FALSE)
plotPCA(vsd,"Sex")
plotPCA(vsd,"Smoking")
plotPCA(vsd,"Hardy")
plotPCA(vsd,"AgeBin")

data = data.frame(fread("GTEX-counts.filterTissue.geneSum.txt")) # raw count ile
info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1) ## note that the detailed clinical annotations are controlled access, so this code will not work with the provided sample table. An anonymized/randomized example of the full sample table is provided: "sample-info.filterTissue.anonymized.full.txt".
dataf<-cbind(data[,1], data[,colnames(data) %in% rownames(info)])

info2<-info[match(colnames(dataf)[2:ncol(dataf)],rownames(info)),]
all(rownames(info2) == colnames(dataf)[2:ncol(dataf)])

colnames(dataf)[1]<-"Name"
table(cbind(info2[4:5]))

infoR<-info2[!is.na(info2$Smoking),]
infoR = infoR[!is.na(infoR$Hardy),]

infoR$AgeBin<-str_replace(infoR$AgeBin,"-","_")
table(cbind(infoR[4:5]))

datafR<-cbind(dataf[,1],dataf[,colnames(dataf) %in% rownames(infoR)])
all(rownames(infoR) == colnames(datafR)[2:ncol(datafR)])
colnames(vsd)
all(rownames(infoR) == colnames(vsd))

infoR$Hardy = as.factor(infoR$Hardy)


assay(vsd) = limma::removeBatchEffect(assay(vsd),covariates=model.matrix(~Sex+Smoking+Hardy,data=infoR),design = model.matrix(~AgeBin,data=infoR))
write.table(format(as.data.frame(assay(vsd)),digits=3),"GTEX-vst.adjusted.txt",sep="\t",row.names=TRUE) # vst transformed and adjusted expression matrix

