#!/usr/bin/Rscript
library(Seurat)

data<-readRDS("lung_ts.rds")
metadata<-data[[]]

#read in the table of age-associated genes, from earlier
agegenes<-read.table("degreport-p0.05-cov-genes-table.upC.downC.txt",sep="\t",header=TRUE)
nyd<-FetchData(object=data,vars=agegenes[,1])

#236 genes not matched, 1049 genes remaining
all(rownames(metadata)==rownames(nyd))

library(NMF)

###### %-based heatmaps and tables
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

### Write output tables of the % expressing
# percentages
dfout<-mydf
dfout<-cbind(rownames(mydf),mydf)
colnames(dfout)[1]<-"gene"

upData<-dfout[rownames(dfout) %in% agegenes[agegenes$cluster == 1,1],]
downData<-dfout[rownames(dfout) %in% agegenes[agegenes$cluster == 2,1],]
write.table(upData,"pctg-expressing.ageUp.avgCelltype.txt",sep="\t",row.names=FALSE)
write.table(downData,"pctg-expressing.ageDown.avgCelltype.txt",sep="\t",row.names=FALSE)

#scaled percentages 
mydfs<-t(scale(t(mydf)))
dfouts<-mydfs
dfouts<-cbind(rownames(mydfs),mydfs)
colnames(dfouts)[1]<-"gene"

upDatas<-dfouts[rownames(dfouts) %in% agegenes[agegenes$cluster == 1,1],]
downDatas<-dfouts[rownames(dfouts) %in% agegenes[agegenes$cluster == 2,1],]
write.table(upDatas,"pctg-expressing.Scaled.ageUp.avgCelltype.txt",sep="\t",row.names=FALSE)
write.table(downDatas,"pctg-expressing.Scaled.ageDown.avgCelltype.txt",sep="\t",row.names=FALSE)



### scaled % heatmaps
#split into upreg and downreg gene clusters

mydfsUp<-mydfs[rownames(mydfs) %in% agegenes[agegenes$cluster == 1,1],] #460 genes
mydfsDown<-mydfs[rownames(mydfs) %in% agegenes[agegenes$cluster == 2,1],] #589 genes

pdf("age-genes-pctgExpressing-zScore-heatmap.prescale.Up.pdf",height=7,width=10,onefile = FALSE)
aheatmap(mydfsUp,color="-RdBu:100",distfun="manhattan",hclustfun="ward")
dev.off()

pdf("age-genes-pctgExpressing-zScore-heatmap.prescale.Down.pdf",height=7,width=10,onefile = FALSE)
aheatmap(mydfsDown,color="-RdBu:100",distfun="manhattan",hclustfun="ward")
dev.off()

### no scaling
#mydfUp<-mydf[rownames(mydf) %in% agegenes[agegenes$cluster == 1,1],]
#mydfDown<-mydf[rownames(mydf) %in% agegenes[agegenes$cluster == 2,1],]

#hmnU<-aheatmap(t(mydfUp),color="-RdBu:100",distfun="manhattan",hclustfun="ward")
#hmnD<-aheatmap(t(mydfDown),color="-RdBu:100",distfun="manhattan",hclustfun="ward")
