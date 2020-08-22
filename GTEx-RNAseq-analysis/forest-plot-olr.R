#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
data=read.table("CIBERSORTx-sBatch-OLR-results.nonNA.txt",sep="\t",header=TRUE,row.names=1)  # OLR output table
data=data[order(data$pval),]
data$cell = factor(rownames(data),levels=rev(rownames(data)))
data$nlp = log(data$pval)*-1
updata = data[data$OLR_estimate > 0 & data$pval < 0.05,]
downdata = data[data$OLR_estimate < 0 & data$pval < 0.05,]
otherdata = data[data$OLR_estimate ==0 | data$pval > 0.05,]
data$color = rep("black",nrow(data))

data[data$cell %in% updata$cell,"color"] = "red3"
data[data$cell %in% downdata$cell,"color"] = "dodgerblue3"

pdf("OLR-CIBERSORTx-sBatch.forest.pdf",height=10,width=10,useDingbats=FALSE)
ggplot(data,aes(x=OLR_estimate,y=cell))+geom_vline(xintercept=0,linetype=2,color="gray44")+geom_errorbar(aes(xmin=lowerCI,xmax=upperCI),width=0.2)+geom_point(aes(color=color,size=nlp))+theme_bw()+scale_color_manual(values=c("black","dodgerblue3","red3")) + scale_size_continuous(name = "-log pvalue",breaks = c(2,4,6,8,12),labels=c("2","4","6","8","12"),limits=c(0,10),range=c(1,8))
dev.off()
