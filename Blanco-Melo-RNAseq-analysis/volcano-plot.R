#!/usr/bin/Rscript
library(ggpubr)

data<-read.table("c-SARS2-vs-Ctrl.DE.wald.Aging.txt",sep=""),sep="\t",header=TRUE,row.names=1)

data<-data[complete.cases(data$padj), ]
data[data$padj==0,]$padj<-10e-300
data$nlp<-log10(data$padj)*-1
data$AgingCluster<-as.factor(data$AgingCluster)
downdata<-data[data$padj < 0.05 & data$AgingCluster == 2,]
updata<-data[data$padj < 0.05 & data$AgingCluster == 1,]
normdata<-data[data$padj >= 0.05 | data$AgingCluster == 0,]
data$genes<-rownames(data)
dim(downdata)
dim(updata)

colScale <- scale_colour_manual(values=c("#00AFBB","#FC4E07","gray77"))
pdf("c-volcano-SARS2-vs-Mock.nolab.pdf",sep=""),height=5.5,width=6,useDingbats=FALSE)
g<-ggplot(data,aes(data$log2FoldChange,data$nlp))+geom_point(pch=19,alpha=0.7,cex=1,data=normdata,aes(x=normdata$log2FoldChange,y=normdata$nlp,color="gray77"))+geom_point(pch=19,alpha=0.85,cex=1,data=downdata,aes(x=downdata$log2FoldChange,y=downdata$nlp,color="#FC4E0"))+geom_point(pch=19,alpha=0.85,cex=1,data=updata,aes(x=updata$log2FoldChange,y=updata$nlp,color="#00AFBB"))+colScale+theme_classic() +theme_pubr()+theme(legend.position="none")
g
dev.off()

# overlap counts for Venn/statistics
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