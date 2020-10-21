#!/usr/bin/Rscript
library(MASS)
library(ordinal)

results<-read.table("CIBERSORTx_Adjusted.txt",sep="\t",header=TRUE,row.names=1)
df<-results[,1:(length(results)-3)]

#filter to cell types that are > 0 abundane in >50% samples
filtered = colnames(df[,(colSums(df != 0)/nrow(df)) < 0.5])
df = df[,(colSums(df != 0)/nrow(df)) >= 0.5]


info<-read.table("full-sample-info.2.txt",sep="\t",header=TRUE,row.names=1) # not provided, controlled access data
data<-t(df)
infof<-info[rownames(info)%in%colnames(data),]
infof = infof[!is.na(infof$Hardy) & !is.na(infof$Smoking) & !is.na(infof$Diabetes) & !is.na(infof$HTN),]
infof$Hardy = as.factor(infof$Hardy)
infof$Sex = as.factor(infof$Sex)
infof$Smoking = as.factor(infof$Smoking)
infof$Diabetes = as.factor(infof$Diabetes)
infof$HTN = as.factor(infof$HTN)

dataf = data[,colnames(data)%in% rownames(infof)]
infof = infof[match(rownames(infof),colnames(dataf)),]
data = dataf[rowSums(dataf) > 0,]
all(rownames(infof)==colnames(data))

data = as.data.frame(t(data))
data$ageBin = infof$AgeBin

#adjust for covariates
mydf = as.data.frame(matrix(nrow=nrow(data),ncol=ncol(data)-1))
library(plyr)
for (i in 1:(ncol(data)-1)){
    c<-colnames(data)[i]
    myCell<-c
    mydata<-as.matrix(data[,colnames(data)==myCell])
    cres<-glm(mydata[,1] ~ (Sex+Smoking+Hardy),data=infof)
    mydf[,i] = as.vector(cres$residuals)
    
}
colnames(mydf)<-c(colnames(data)[-ncol(data)])
mydf$ageBin = data$ageBin

#get cell type order from OLR results
results = read.table("CIBERSORTx-OLR-age-results.COVID.final.nonNA.txt",sep="\t",header=TRUE,row.names=1)
celltypes = rownames(results)
plot_list = list()

#unadjusted boxplots
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr)
for (i in 1:length(celltypes)){
    mycell = celltypes[i]
    mydata = mydf[,c(mycell,"ageBin")]
    g = ggplot(data=mydata,aes_string(x="ageBin",y=mycell))+geom_boxplot(aes(fill=ageBin))+scale_fill_nejm()+theme_pubr()+theme(legend.position="none")+ggtitle(label=mycell)+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),plot.title = element_text(hjust = 0.5))
    plot_list[[i]]= g
}

pdf("adjusted-CIBERSORTx-boxplots-agebin.pdf",height=20,width=20,useDingbats=FALSE)
do.call(grid.arrange,plot_list)
dev.off()
