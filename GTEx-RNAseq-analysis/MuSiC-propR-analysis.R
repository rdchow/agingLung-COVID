#!/usr/bin/Rscript
library(propr)
library(NMF)
library(ggplot2)
library(ggpubr)

results<-readRDS("GTEX-music-analysis.rds")
df<-results$Est.prop.weighted
df = df[,colSums(df) != 0] #remove cell types with proportion = 0 across all samples

info<-read.table("sample-info.filterTissue.rn.txt",sep="\t",header=TRUE,row.names=1)
infof<-info[rownames(info)%in%rownames(df),]
all(rownames(infof)==rownames(data))

pd.d = propd(df,infof$Age,alpha=0.5,weighted=TRUE,p=100)
pd.d <- updateCutoffs(pd.d, cutoff = seq(0.9, 1, .01))
pd.d <- updateF(pd.d, moderated = TRUE,ivar = "clr")

tab <- getResults(pd.d)
write.table(tab, "propR-analysis-table.txt",sep="\t",row.names=FALSE)

mat = getMatrix(pd.d)
matSig = mat
matSig[matSig > 0.95] = NA # filter out pairs with theta > 0.95
matSig[matSig == 0] = NA # filter out homotypic pairs since these are not informative
ord = read.table("propr-row-order.txt",sep="\t",header=FALSE)
ord = as.vector(ord[ord[,1] %in% colnames(matSig),])
matSigR = matSig[match(ord,rownames(matSig)),]
matSigR = matSigR[,match(ord,colnames(matSigR))]

aheatmap(matSigR,color="-YlOrRd:100",Colv=NA,Rowv=NA,border="gray77",width=3,height=3)


