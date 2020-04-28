#!/usr/bin/Rscript
library(Seurat)
data<-readRDS("lung_ts.rds")
metadata<-data[[]]

## SARS-CoV-2 PPI genes, with age-association
wantgenes<-read.table("PPI-annotated-intersection-genes.txt",sep="\t",header=TRUE)

## SARS-CoV siRNA screen genes, with age-association 
#wantgenes<-read.table("SARS-siRNA-screen-annotated-intersection-genes.txt",sep="\t",header=TRUE)

wantgenes<-wantgenes[complete.cases(wantgenes),]
wantgenes1<-wantgenes[wantgenes$Age_cluster == 1,]
wantgenes2<-wantgenes[wantgenes$Age_cluster == 2,]
coolgenes<-c(as.vector(wantgenes1$gene),as.vector(wantgenes2$gene))

nyd<-FetchData(object=data,vars=coolgenes) 

library(NMF)
all(rownames(metadata)==rownames(nyd))

###### % based version
#calculate % of cells expressing each gene

bindata<-nyd
bindata[bindata>0]<-1
tbindata<-t(bindata)
celltypes<-unique(metadata$Celltypes)

mydf<-matrix(nrow=nrow(tbindata),ncol=length(celltypes))
for (i in 1:nrow(tbindata)){
    mydata<-tbindata[i,]
    mygene<-rownames(tbindata)[i]
    stats<-t(table(mydata,metadata$Celltypes))

    pctg<-stats[,2]/(stats[,1]+stats[,2])*100
    mydf[i,]<-pctg
}
rownames(mydf)<-rownames(tbindata)
colnames(mydf)<-rownames(stats)

mydf1<-mydf[rownames(mydf) %in% wantgenes1$gene,]
mydf2<-mydf[rownames(mydf) %in% wantgenes2$gene,]

dim(mydf1)
dim(mydf2)

pdf("PPI-ageUp-intersect-heatmap.pdf",height=8,width=7,onefile=FALSE)
aheatmap(t(mydf1),color="RdPu:100")
dev.off()

pdf("PPI-ageDown-intersect-heatmap.pdf",height=8,width=7,onefile=FALSE)
aheatmap(t(mydf2),color="RdPu:100")
dev.off()