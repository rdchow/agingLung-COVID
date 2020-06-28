#!/usr/bin/Rscript

library(Seurat)
c141<-Read10X_h5("./h5-files/GSM4339769_C141_filtered_feature_bc_matrix.h5")
c142<-Read10X_h5("./h5-files/GSM4339770_C142_filtered_feature_bc_matrix.h5")
c143<-Read10X_h5("./h5-files/GSM4339771_C143_filtered_feature_bc_matrix.h5")
c144<-Read10X_h5("./h5-files/GSM4339772_C144_filtered_feature_bc_matrix.h5")
c145<-Read10X_h5("./h5-files/GSM4339773_C145_filtered_feature_bc_matrix.h5")
c146<-Read10X_h5("./h5-files/GSM4339774_C146_filtered_feature_bc_matrix.h5")

c51<-Read10X_h5("./h5-files/GSM4475048_C51_filtered_feature_bc_matrix.h5")
c52<-Read10X_h5("./h5-files/GSM4475049_C52_filtered_feature_bc_matrix.h5")
c100<-Read10X_h5("./h5-files/GSM4475050_C100_filtered_feature_bc_matrix.h5")
c148<-Read10X_h5("./h5-files/GSM4475051_C148_filtered_feature_bc_matrix.h5")
c149<-Read10X_h5("./h5-files/GSM4475052_C149_filtered_feature_bc_matrix.h5")
c152<-Read10X_h5("./h5-files/GSM4475053_C152_filtered_feature_bc_matrix.h5")

c51s<-CreateSeuratObject(c51)
#saveRDS(c51s,"c51.seurat.rds")
c52s<-CreateSeuratObject(c52)
c100s<-CreateSeuratObject(c100)
c148s<-CreateSeuratObject(c148)
c149s<-CreateSeuratObject(c149)
c152s<-CreateSeuratObject(c152)
c141s<-CreateSeuratObject(c141)
c142s<-CreateSeuratObject(c142)
c143s<-CreateSeuratObject(c143)
c144s<-CreateSeuratObject(c144)
c145s<-CreateSeuratObject(c145)
c146s<-CreateSeuratObject(c146)

test<-merge(x=c141s,y=list(c142s,c143s,c144s,c145s,c146s,c148s,c149s,c152s,c51s,c52s,c100s),add.cell.ids=list("C141","C142","C143","C144","C145","C146","C148","C149","C152","C51","C52","C100"))
test[["percent.mt"]] <- PercentageFeatureSet(test,pattern="^MT-")

test<-subset(test,subset=nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 10 & nCount_RNA >= 1000)
saveRDS(test,file="GSE145926-merged.all.rds")

