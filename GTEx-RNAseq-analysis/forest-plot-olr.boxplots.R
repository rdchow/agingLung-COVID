#!/usr/bin/Rscript
library(ggplot2)
library(ggpubr)
data=read.table("CIBERSORTx-OLR-age-results.COVID.final.nonNA.txt",sep="\t",header=TRUE,row.names=1)
#data=read.table("music-krasnow-OLR-results.nonNA.txt",sep="\t",header=TRUE,row.names=1)
data=data[order(data$pval),]
data$cell = factor(rownames(data),levels=rev(rownames(data)))
data$nlp = log(data$pval)*-1
updata = data[data$OLR_estimate > 0 & data$pval < 0.05,]
downdata = data[data$OLR_estimate < 0 & data$pval < 0.05,]
otherdata = data[data$OLR_estimate ==0 | data$pval > 0.05,]
data$color = rep("black",nrow(data))

data[data$cell %in% updata$cell,"color"] = "#00AFBB"
data[data$cell %in% downdata$cell,"color"] = "#FC4E07"

pdf("OLR-CIBERSORTx-lung.age.v2.pdf",height=10,width=10,useDingbats=FALSE)
ggplot(data,aes(x=OLR_estimate,y=cell))+geom_vline(xintercept=0,linetype=2,color="gray44")+geom_errorbar(aes(xmin=lowerCI,xmax=upperCI),width=0.2)+geom_point(aes(color=color,size=nlp))+theme_bw()+scale_color_manual(values=c("#00AFBB","#FC4E07","black")) + scale_size_continuous(name = "-log pvalue",range=c(1,5))
dev.off()

####################################
# box plot of estimated abundances, no filtering
library(robustbase)
results=read.table("CIBERSORTx_Adjusted.txt",sep="\t",header=TRUE,row.names=1)
df=results[,1:(length(results)-3)]
df = df[,order(colMedians(as.matrix(df)),colMeans(df),decreasing=TRUE)]

library(dplyr)
library(reshape2)
meltdata = melt(df)
meltdata = meltdata %>% group_by(variable) %>% mutate(outlier = value > median(value) + IQR(value)*1.5 | value < median(value) - IQR(value)*1.5) 
#mycols = data$color
#mycols = gsub("black","gray44",mycols)
bp = ggplot(meltdata,aes(x=variable,y=value))+geom_jitter(data=meltdata[meltdata$outlier==TRUE,],aes(color=variable),size=0.5,width=0.2,pch=19)+geom_boxplot(aes(fill=variable),color="black",outlier.shape=NA,size=0.25)+theme_bw() + theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+scale_fill_manual(values=mycols)+scale_color_manual(values=mycols)

# bar plot of % samples > 0 proportion
prop0 = as.data.frame((colSums(df != 0)/nrow(df))*100)
prop0$cell = rownames(prop0)
prop0$cell = factor(prop0$cell,levels=colnames(df))
prop0$data = prop0[,1]
pbp = ggplot(prop0,aes(x=cell,y=data))+geom_bar(stat="identity",aes(fill=cell))+theme_bw() + theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

library(gridExtra)
pdf("CIBERSORTx-lung-estimated-proportions.boxplot.pdf",height=9,width=12,useDingbats=FALSE)
grid.arrange(bp,pbp,nrow=2,heights=c(1.5,1))
dev.off()