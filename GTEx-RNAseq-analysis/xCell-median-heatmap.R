#!/usr/bin/Rscript

info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
data<-read.table("xCell-results.filterTissue.txt",sep="\t",header=TRUE,row.names=1)
infof<-info[rownames(info)%in%colnames(data),]
all(rownames(infof)==colnames(data))

myTissue<-"Lung"
infof2<-infof[infof$simpleTissue == myTissue,]
dataf<-data[,colnames(data)%in%rownames(infof2)]
all(rownames(infof2)==colnames(dataf))
table(cbind(infof2[4:5]))

scaledata<-t(scale(t(dataf)))

aggdata<-t(aggregate(t(scaledata),by=list(Age=infof2[,5]),FUN=median))
write.table(aggdata,"xcell-analysis_median-zscore-by-age.txt",sep="\t",row.names=TRUE)

library(superheat)
finaldata<-read.table("xcell-analysis_median-zscore-by-age.txt",sep="\t",header=TRUE,row.names=1)

celltypes<-c("Fibroblasts","Epithelial cells","Neutrophils","Eosinophils","B-cells","CD4+ memory T-cells","CD4+ naive T-cells","CD4+ T-cells","CD4+ Tcm","CD4+ Tem","CD8+ naive T-cells","CD8+ T-cells", "CD8+ Tcm","CD8+ Tem","aDC","cDC","DC","iDC","pDC","CLP","CMP","Endothelial cells","Macrophages M1","Macrophages M2","Macrophages","Memory B-cells","naive B-cells","NK cells","NKT","Plasma cells","Tgd cells","Th1 cells","Th2 cells","Tregs","Pericytes")
finaldataf<-finaldata[rownames(finaldata)%in%celltypes,]

cellorder<-read.table("xcell-order.txt",sep="\t",header=FALSE)
finaldatafr<-finaldataf[rev(match(cellorder[,1],rownames(finaldataf))),]


png("lung-xcell-heatmap.png",height=4000,width=2400,res=450)
superheat(finaldatafr,scale = FALSE, extreme.values.na = FALSE,left.label.size=1,left.label.text.size=3,bottom.label.text.size=4,grid.vline = FALSE,grid.hline = FALSE, left.label.col = "white", bottom.label.col = "white", left.label="variable")
dev.off()


