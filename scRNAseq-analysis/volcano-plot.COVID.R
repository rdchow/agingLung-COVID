#!/usr/bin/Rscript

data<-read.table("GSE145926-severe_COVID-vs-Control.de.Aging.txt",sep="\t",header=TRUE,row.names=1)

#1025 genes in both aging list and DEseq results, before filtering for pvalue
dim(data[data$AgingCluster != 0,])

data<-data[complete.cases(data$padj), ]
#data[data$padj==0,]$padj<-10e-300


data$nlp<-log10(data$padj)*-1
data$AgingCluster<-as.factor(data$AgingCluster)
downdata<-data[data$padj < 0.05 & data$AgingCluster == 2,]
updata<-data[data$padj < 0.05 & data$AgingCluster == 1,]
normdata<-data[data$padj >= 0.05 | data$AgingCluster == 0,]
data$genes<-rownames(data)
dim(downdata)
dim(updata)


library(ggpubr)
library(ggrepel)



colScale <- scale_colour_manual(values=c("#00AFBB","#FC4E07","gray77"))
pdf("GSE145926-volcano-severeSARS2-vs-ctrl.nolab.pdf",height=5.5,width=6,useDingbats=FALSE)
g<-ggplot(data,aes(data$log2FoldChange,data$nlp))+geom_point(pch=19,alpha=0.7,cex=1,data=normdata,aes(x=normdata$log2FoldChange,y=normdata$nlp,color="gray77"))+geom_point(pch=19,alpha=0.85,cex=1,data=downdata,aes(x=downdata$log2FoldChange,y=downdata$nlp,color="#FC4E0"))+geom_point(pch=19,alpha=0.85,cex=1,data=updata,aes(x=updata$log2FoldChange,y=updata$nlp,color="#00AFBB"))+colScale+theme_classic() +theme_pubr()+ theme(legend.position="none")
g
dev.off()


# stats calculation
overlapdata<-data[data$padj < 0.05 & data$AgingCluster != 0,]
nrow(overlapdata)

sub1<-overlapdata[overlapdata$log2FoldChange >0 & overlapdata$AgingCluster == 1,]
nrow(sub1)

sub2<-overlapdata[overlapdata$log2FoldChange >0 & overlapdata$AgingCluster == 2,]
nrow(sub2)

sub3<-overlapdata[overlapdata$log2FoldChange <0 & overlapdata$AgingCluster == 1,]
nrow(sub3)

sub4<-overlapdata[overlapdata$log2FoldChange <0 & overlapdata$AgingCluster == 2,]
nrow(sub4)

sigdata<-data[data$padj < 0.05,]
dim(sigdata[sigdata$log2FoldChange >0,])
dim(sigdata[sigdata$log2FoldChange <0,])

write.table(overlapdata,"severe-overlap-genes.txt",sep="\t",row.names=TRUE)
write.table(sigdata,"severe-significant-genes.txt",sep="\t",row.names=TRUE)
