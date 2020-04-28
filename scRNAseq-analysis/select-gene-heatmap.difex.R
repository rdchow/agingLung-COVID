#!/usr/bin/Rscript
library(Seurat)
setwd("C:/CRISPR/COVID19/human-lung-scRNA/tenover/")

data<-readRDS("../lung_ts.rds")
metadata<-data[[]]
#write.table(metadata,"seurat-metadata.txt",sep="\t",row.names=TRUE)

#nyd<-FetchData(object = data, vars = c("ACE2","TMPRSS2","CTSL"))
#all(rownames(metadata)==rownames(nyd))
#nydf<-cbind(nyd,metadata)
#write.table(nydf,"entry-factor-expr.annot.txt",sep="\t",row.names=TRUE)

#pdf("Entry-factors-violin.0pt.pdf",height=5,width=20,useDingbats=FALSE)
#VlnPlot(data,features=c("ACE2","TMPRSS2","CTSL"),pt.size=0)
#dev.off()

wantgenes<-read.table("Difex-annotated-intersection-genes.txt",sep="\t",header=TRUE)
wantgenes<-wantgenes[complete.cases(wantgenes),]
ga<-wantgenes[wantgenes$category == "a",]
gb<-wantgenes[wantgenes$category == "b",]
gc<-wantgenes[wantgenes$category == "c",]
gd<-wantgenes[wantgenes$category == "d",]

coolgenes<-c(as.vector(ga$Gene),as.vector(gb$Gene),as.vector(gc$Gene),as.vector(gd$Gene))




nyd<-FetchData(object=data,vars=wantgenes[,1]) 

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

scalemydf<-scale(t(mydf))

mydf1<-mydf[rownames(mydf) %in% ga$Gene,]
mydf2<-mydf[rownames(mydf) %in% gb$Gene,]
mydf3<-mydf[rownames(mydf) %in% gc$Gene,]
mydf4<-mydf[rownames(mydf) %in% gd$Gene,]

dim(mydf1)
dim(mydf2)
dim(mydf3)
dim(mydf4)

library(viridis)
mycols = magma(100)
d1<-aheatmap(t(mydf1),color=mycols,distfun="manhattan",hclustfun="ward")
d2<-aheatmap(t(mydf2),color=mycols,distfun="manhattan",hclustfun="ward")
d3<-aheatmap(t(mydf3),color=mycols,distfun="manhattan",hclustfun="ward")
d4<-aheatmap(t(mydf4),color=mycols,distfun="manhattan",hclustfun="ward")

d1order<-rownames(mydf1[d1$colInd,])
d2order<-rownames(mydf2[d2$colInd,])
d3order<-rownames(mydf3[d3$colInd,])
d4order<-rownames(mydf4[d4$colInd,])

geneorder<-c(d1order,d2order,d3order,d4order)
orderInd<-match(geneorder,rownames(mydf))

#all together
library(ggsci)
Var1=pal_lancet()(4)
Var1=c("#00AFBB","#161A92","#FC4E07","#782F0D")
#Var1=c( "#DF8F44FF", "#374E55FF",  "#B24745FF", "#00A1D5FF")
names(Var1) = c("a", "b","c","d")
ann_colors = list(Var1 = Var1)
annotation = data.frame(Var1 = wantgenes$category)
pdf("difex-age-intersect-heatmap.pdf",height=5,width=12,onefile=FALSE)
aheatmap(t(mydf),color=mycols,Colv=orderInd,annCol=annotation,annColors=ann_colors)
dev.off()
