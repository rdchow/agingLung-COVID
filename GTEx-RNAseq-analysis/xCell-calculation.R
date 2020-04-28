#!/usr/bin/Rscript
########### 
# Run xCell
library(xCell)
library(data.table)
exprMatrix = data.frame(fread("GTEX-TPM.filterTissue.geneSum.txt"),row.names=1)
results<-xCellAnalysis(exprMatrix)
results1<-cbind(rownames(results),results)
colnames(results1)<-c("cellType",colnames(results))
write.table(results1,"xCell-results.filterTissue.txt",sep="\t",row.names=FALSE)

###############
# Statistical analysis of the xCell table
#runs Kruskal-Wallis test over each row of the xCell results table
info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
data<-read.table("xCell-results.filterTissue.txt",sep="\t",header=TRUE,row.names=1)
infof<-info[rownames(info)%in%colnames(data),]
all(rownames(infof)==colnames(data))

myTissue<-"Lung"
infof2<-infof[infof$simpleTissue == myTissue,]
dataf<-data[,colnames(data)%in%rownames(infof2)]
all(rownames(infof2)==colnames(dataf))
table(cbind(infof2[4:5]))

mydf<-matrix(nrow=nrow(dataf),ncol=2)
for (i in 1:nrow(dataf)){
    c<-rownames(dataf)[i]
    myCell<-c
    mydata<-as.numeric(dataf[rownames(dataf)==myCell,])
    res<-kruskal.test(mydata~ Age ,data=infof2[4:5])
    mydf[i,1]<-c
    mydf[i,2]<-res$p.value
}
colnames(mydf)<-c("Cell","KW_pval")
write.table(mydf,"lung-kruskal-results.txt",sep="\t",row.names=FALSE)