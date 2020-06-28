#!/usr/bin/Rscript
results<-readRDS("GTEx-music-analysis.rds")

df<-results$Est.prop.weighted

#boxplot(df)
write.table(df,"gtex-music-table.txt",sep="\t",row.names=TRUE)


info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
data<-t(df)
infof<-info[rownames(info)%in%colnames(data),]
all(rownames(infof)==colnames(data))


table(cbind(infof[4:5]))

mydf<-matrix(nrow=nrow(data),ncol=2)
for (i in 1:nrow(data)){
    c<-rownames(data)[i]
    myCell<-c
    mydata<-as.numeric(data[rownames(data)==myCell,])
    res<-kruskal.test(mydata~ Age ,data=infof[4:5])
    mydf[i,1]<-c
    mydf[i,2]<-res$p.value
}
colnames(mydf)<-c("Cell","KW_pval")
write.table(mydf,"music-kruskal-results.txt",sep="\t",row.names=FALSE)