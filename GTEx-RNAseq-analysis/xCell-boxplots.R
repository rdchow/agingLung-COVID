#!/usr/bin/Rscript
setwd("C:/CRISPR/COVID19")
#runs two-way anova over each row of the xCell results table
info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
data<-read.table("xCell-results.filterTissue.txt",sep="\t",header=TRUE,row.names=1)
infof<-info[rownames(info)%in%colnames(data),]
all(rownames(infof)==colnames(data))

myTissue<-"Lung"
infof2<-infof[infof$simpleTissue == myTissue,]
dataf<-data[,colnames(data)%in%rownames(infof2)]
all(rownames(infof2)==colnames(dataf))
table(cbind(infof2[4:5]))

###specific boxplots
library(ggpubr)
library(tidyr)
library(RColorBrewer)
library(ggsci)

celltypes<-c("Fibroblasts","Epithelial cells","B-cells","CD4+ memory T-cells","CD4+ naive T-cells","CD4+ T-cells","CD4+ Tcm","CD4+ Tem","CD8+ naive T-cells","CD8+ T-cells", "CD8+ Tcm","CD8+ Tem","aDC","cDC","DC","iDC","pDC","Class-switched memory B-cells","CLP","CMP","ly Endothelial cells","Endothelial cells","Macrophages M1","Macrophages M2","Macrophages","Memory B-cells","naive B-cells","Neurons","NK cells","NKT","Plasma cells","Tgd cells","Th1 cells","Th2 cells","Tregs")

for (mycell in celltypes){
    myCell<-mycell
    mydata<-as.numeric(dataf[rownames(dataf)==myCell,])
    mydf<-cbind(infof2[4:5],mydata)
    mydf$combined<-unite(mydf[1:2],col="combined",Age,Gender,sep="_",remove=TRUE)[,1]
    mydf$combined<-as.factor(mydf$combined)
    mydf$logdata<-log2(10*mydf$mydata+1)
    mydf$Age<-as.factor(mydf$Age)
    mycolors<-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#de956e", "#B15928")

    #combined
    p<-ggplot(mydf,aes( x = mydf$Age, y = mydf$mydata))+geom_boxplot(outlier.shape=NA,aes(color=mydf$Age))+scale_fill_brewer("Set1")+geom_jitter(aes(color=mydf$Age,shape=mydf$Gender),width=0.25,height=0)+ylab(paste(myCell," abundance",sep=""))+scale_color_nejm()+theme_pubr()

    pdf(paste(myCell,"-lung-boxplots.pdf",sep=""),height=5,width=8,useDingbats=FALSE)
    print(ggpar(p,legend="none")+stat_compare_means())
    dev.off()
}