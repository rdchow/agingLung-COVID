#!/usr/bin/Rscript
library(MASS)
library(ordinal)
results<-read.table("CIBERSORTx_Adjusted.txt",sep="\t",header=TRUE,row.names=1)
df<-results[,1:(length(results)-3)]

#filter to cell types that are > 0 abundance in >50% samples
filtered = colnames(df[,(colSums(df != 0)/nrow(df)) < 0.5])
df = df[,(colSums(df != 0)/nrow(df)) >= 0.5]

boxplot(df)


info<-read.table("full-sample-info.2.txt",sep="\t",header=TRUE,row.names=1) # this is controlled access information, not included here
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


table(cbind(infof[5:6]))

mydf<-matrix(nrow=nrow(data),ncol=5)
for (i in 1:nrow(data)){
    c<-rownames(data)[i]
    myCell<-c
    mydata<-as.numeric(data[rownames(data)==myCell,])
    mydataO = ordered(mydata,levels=c(sort(unique(mydata))))
  
    cres<-clm(mydataO ~ (Age+Sex+Smoking+Hardy),data=infof,link="logit",Hess=TRUE)
    csum = summary(cres)$coef
    mydf[i,1]<-c
    #mydf[i,2:5]<-csum["Diabetesyes",]
    mydf[i,2:5]<-csum["Age",]
}
colnames(mydf)<-c("Cell","OLR_estimate","StdError","z","pval")
#write.table(mydf,"CIBERSORTx-OLR-age-results.final.txt",sep="\t",row.names=FALSE)

mydf2 = as.data.frame(mydf[complete.cases(mydf),])
mydf2$OLR_estimate = as.numeric(as.character(mydf2$OLR_estimate))
mydf2$StdError = as.numeric(as.character(mydf2$StdError))

mydf2$lowerCI = mydf2$OLR_estimate - 1.96*mydf2$StdError
mydf2$upperCI = mydf2$OLR_estimate + 1.96*mydf2$StdError
mydf2$adjp = p.adjust(mydf2$pval,method="BH")
write.table(mydf2[order(mydf2$pval),],"CIBERSORTx-OLR-age-results.COVID.final.nonNA.txt",sep="\t",row.names=FALSE)
