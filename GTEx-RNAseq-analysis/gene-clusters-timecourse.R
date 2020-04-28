#!/usr/bin/Rscript
#!/usr/bin/Rscript
setwd("C:/CRISPR/COVID19")
library(data.table)
library(ggpubr)
library(tidyr)
library(RColorBrewer)
library(ggsci)
#runs two-way anova over each row of the xCell results table
info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
data<-data.frame(fread("GTEX-TPM.filterTissue.geneSum.txt"),row.names=1)
infof<-info[rownames(info)%in%colnames(data),]
all(rownames(infof)==colnames(data))

results<-read.table("GTEx-lung-DEseq2-results.LRT.age.txt",sep="\t",header=TRUE)
results<-results[complete.cases(results),]
filtresults<-results[results$padj<0.0001,]
wantgenes<-as.vector(filtresults$row)

myTissue<-"Lung"
infof2<-infof[infof$simpleTissue == myTissue,]
dataf<-data[,colnames(data)%in%rownames(infof2)]
dataf2<-dataf[rownames(dataf)%in%wantgenes,]
dataf3<-log2(dataf2+1)
all(rownames(infof2)==colnames(dataf3))
table(cbind(infof2[4:5]))

library(DEGreport)
clusters <- degPatterns(dataf2, metadata = infof2[5], time = "Age")
filtclusters<-clusters$normalized
filtclusters<-filtclusters[filtclusters$cluster == "1" | filtclusters$cluster == "2",]

pdf("degreport-p0.0001-genes-color-boxplot.lines.pdf",height=3,width=8,useDingbats=FALSE)
degPlotCluster(filtclusters,time="Age",color="Age")+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()
dev.off()
myTable<-clusters$df
myTable<-myTable[myTable$cluster == "1" | myTable$cluster == "2",]
write.table(myTable,"degreport-p0.0001-genes-table.c1-c2-only.txt",sep="\t",row.names=FALSE)