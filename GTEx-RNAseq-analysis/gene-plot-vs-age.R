#!/usr/bin/Rscript

library(data.table)
info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
data<-data.frame(fread("GTEX-TPM.filterTissue.geneSum.txt"),row.names=1)
infof<-info[rownames(info)%in%colnames(data),]
all(rownames(infof)==colnames(data))

myTissue<-"Lung"
infof2<-infof[infof$simpleTissue == myTissue,]
dataf<-data[,colnames(data)%in%rownames(infof2)]
all(rownames(infof2)==colnames(dataf))
table(cbind(infof2[4:5]))

##Generate specific boxplots
library(ggpubr)
library(tidyr)
library(RColorBrewer)
library(ggsci)

g<-list()

wantgenes<-c("ACE2","TMPRSS2","CTSL")
for (i in 1:length(wantgenes)){
    myGene<-wantgenes[i]
    mydata<-as.numeric(dataf[rownames(dataf)==myGene,])
    mydf<-cbind(infof2[4:5],mydata)
    mydf$combined<-unite(mydf[1:2],col="combined",Age,Gender,sep="_",remove=TRUE)[,1]
    mydf$combined<-as.factor(mydf$combined)
    mydf$logdata<-log2(mydf$mydata+1)
    mydf$Age<-as.factor(mydf$Age)
    mycolors<-c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#de956e", "#B15928")

    #Boxplots, split by age (gender with different point shapes)
    p<-ggplot(mydf,aes_string( x = mydf$Age, y = mydf$logdata))+geom_boxplot(outlier.shape=NA,aes(color=mydf$Age),size=1.25)+scale_fill_brewer("Set1")+geom_jitter(aes(color=mydf$Age,shape=mydf$Gender),width=0.25,height=0,size=2,alpha=0.5)+ylab(paste(myGene," expression (log2 TPM)",sep=""))+scale_color_nejm()+theme_pubr()

    #Statistical analysis
    my_comparisons <- list( c("20-29", "30-39"),c("20-29", "40-49"),c("20-29", "50-59"), c("20-29", "60-69"), c("20-29", "70-79") )


    #version with pairwise comparison and Kruskal-Wallis test for all ages
    #g[[i]]<-ggpar(p,legend="none")+stat_compare_means(comparisons = my_comparisons,method="wilcox.test")+stat_compare_means()

    #version without pairwise comparisons
    g[[i]]<-ggpar(p,legend="none")+stat_compare_means()
}

library(gridExtra)

pdf(paste("SARS2-hostFactor-boxplots.pdf",sep=""),height=5,width=20,useDingbats=FALSE
gt<-arrangeGrob(g[[1]],g[[2]],g[[3]],ncol=3)
as_ggplot(gt)
dev.off()