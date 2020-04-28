#!/usr/bin/Rscript

#sums up rows with the same gene annotation
install.packages("data.table")
library(data.table)
data<-data.frame(fread("GTEX-TPM.filterTissue.txt",sep="\t"))
sumdata<-aggregate(data[-1],by=list(data$Name),FUN=sum,na.rm=TRUE)
write.table(sumdata,"GTEX-TPM.filterTissue.geneSum.txt",sep="\t",row.names=FALSE)


data<-data.frame(fread("GTEX-counts.filterTissue.txt",sep="\t"))
sumdata<-aggregate(data[-1],by=list(data$Name),FUN=sum,na.rm=TRUE)
write.table(sumdata,"GTEX-counts.filterTissue.geneSum.txt",sep="\t",row.names=FALSE)