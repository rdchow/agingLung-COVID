#!/usr/bin/Rscript
#!/usr/bin/Rscript
library(data.table)
library(ggpubr)
library(tidyr)
library(RColorBrewer)
library(ggsci)

info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
data<-data.frame(fread("GTEX-vst.adjusted.txt"),row.names=1) # this is the VST and limma adjusted expression matrix
infof<-info[rownames(info)%in%colnames(data),]
all(rownames(infof)==colnames(data))

results<-read.table("GTEx-lung-DEseq2-results.LRT.timeBin.combined.SexSmokingIschemia.all.txt",sep="\t",header=TRUE) #the DESeq2 differential expression table. Note that the p-values are the same for each of the pairwise contrasts, so this does not need to be the combined DE table
results<-results[complete.cases(results),]
filtresults<-results[results$padj<0.05,]
wantgenes<-as.vector(filtresults$row)
dataf = data[rownames(data) %in% wantgenes,]

library(DEGreport)
clusters <- degPatterns(dataf, metadata = infof, time = "Age")
filtclusters<-clusters$normalized
upc = c("1","3","6")
downc = c("12","15")
filtclusters2<-filtclusters[filtclusters$cluster %in% upc | filtclusters$cluster %in% downc,]
filtclusters2[filtclusters2$cluster %in% upc,"cluster"] = 1
filtclusters2[filtclusters2$cluster %in% downc,"cluster"] = 2


pdf("degreport-p0.05-cov-genes-color-boxplot.lines.pdf",height=3,width=8,useDingbats=FALSE)
degPlotCluster(filtclusters2,time="Age",color="Age")+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()
dev.off()

mydf = clusters$df
write.table(mydf,"degreport-p0.05-cov-genes-table.all.txt",sep="\t",row.names=FALSE)

mydf3 = mydf[mydf$cluster %in% upc | mydf$cluster %in% downc,]
mydf3[mydf3$cluster %in% upc,2] = "1"
mydf3[mydf3$cluster %in% downc,2] = "2"

write.table(mydf3,"degreport-p0.05-cov-genes-table.upC.downC.txt",sep="\t",row.names=FALSE)