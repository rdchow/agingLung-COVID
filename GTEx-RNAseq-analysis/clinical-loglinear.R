#!/usr/bin/Rscript
library(dplyr)

data = read.table("anonymized-clinical-annot.txt",sep="\t",header=TRUE,row.names=1,na.strings=c("","NA"))

#take the categorical variables
dataf = data[complete.cases(data),]


dataf2 = as.data.frame(dataf %>% group_by(AgeBin,Sex,Obesity,HTN,T1D,T2D,Vent,Smoking) %>% tally())

#log-linear model, no interaction terms
model0 = glm(n ~ (AgeBin + Sex + Obesity + HTN + T1D + T2D + Vent + Smoking), data = dataf2, family="poisson")
summary(model0)
pchisq(deviance(model0), df = df.residual(model0), lower.tail = F)

##log-linear model, with interaction terms involving age
model2 = glm(n ~ (AgeBin*Sex + AgeBin*Obesity + AgeBin*HTN + AgeBin*T1D + AgeBin*T2D + AgeBin*Vent + AgeBin*Smoking), data = dataf2, family="poisson")
summary(model2)
pchisq(deviance(model2), df = df.residual(model2), lower.tail = F)
