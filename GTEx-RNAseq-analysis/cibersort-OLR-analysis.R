#!/usr/bin/Rscript
library(ordinal)

results<-read.table("CIBERSORTx_sBatch.txt",sep="\t",header=TRUE,row.names=1) # CIBERSORTx output matrix, with S-mode batch correction
df<-results[,1:57] # keep the relevant columns
boxplot(df)

info<-read.table("sample-info.filterTissue.anonymized.full.txt",sep="\t",header=TRUE,row.names=1) # need the detailed clinical annotation file; the provided file is anonymized.
data<-t(df)
infof<-info[rownames(info)%in%colnames(data),]
infof = infof[!is.na(infof$Hardy) & !is.na(infof$Smoking),]
infof$Hardy = as.factor(infof$Hardy)

dataf = data[,colnames(data)%in% rownames(infof)]
infof = infof[match(rownames(infof),colnames(dataf)),]
data = dataf[rowSums(dataf) > 0,] # remove cell types that were "0" across all samples

all(rownames(infof)==colnames(data))


table(cbind(infof[4:5]))

mydf<-matrix(nrow=nrow(data),ncol=5)
for (i in 1:nrow(data)){
    c<-rownames(data)[i]
    myCell<-c
    mydata<-as.numeric(data[rownames(data)==myCell,])
    mydataO = ordered(mydata,levels=c(sort(unique(mydata))))
    
    cres<-clm(mydataO ~ (Age+Sex+Smoking+Hardy),data=infof,link="logit",Hess=TRUE)
    csum = summary(cres)$coef
    mydf[i,1]<-c
    mydf[i,2:5]<-csum["Age",]
}
colnames(mydf)<-c("Cell","OLR_estimate","StdError","z","pval")
write.table(mydf,"CIBERSORTx-sBatch-OLR-results.txt",sep="\t",row.names=FALSE)

mydf2 = as.data.frame(mydf[complete.cases(mydf),])
mydf2$OLR_estimate = as.numeric(as.character(mydf2$OLR_estimate))
mydf2$StdError = as.numeric(as.character(mydf2$StdError))

mydf2$lowerCI = mydf2$OLR_estimate - 1.96*mydf2$StdError
mydf2$upperCI = mydf2$OLR_estimate + 1.96*mydf2$StdError
write.table(mydf2,"CIBERSORTx-sBatch-OLR-results.nonNA.txt",sep="\t",row.names=FALSE)