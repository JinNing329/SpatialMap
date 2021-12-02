
rm(list=ls())
library(Seurat)
library(xlsx)

####################### STEP1 :load ##################
counts <- Read10X("MERFISH/GSE113576/")
counts <- as.matrix(counts)
scMeta <- read.xlsx("aau5324_Moffitt_Table-S1.xlsx",sheetIndex = 1)
colnames(scMeta)=c("Cell_name", "Sex", "Repulicate_number", "Cell_class", "Non_neuron_class", "Neuron_class")
scMeta=scMeta[-1,]
## Remove cell with class "Ambiguous" or "Unstable"
ind1=which(scMeta$Cell_class=="Ambiguous"|scMeta$Cell_class=="Unstable")
scMeta=scMeta[-ind1,]
counts=counts[,-ind1]
## Specify the class
scMeta$Class=as.character(scMeta$Non_neuron_class)
ind2=which(is.na(scMeta$Class))
scMeta$Class[ind2]=as.character(scMeta$Neuron_class[ind2])
rownames(scMeta)=scMeta$Cell_name
## filtering gene
counts=counts[which(rowSums(counts)>100),]

####################### STEP2 :CreateSeuratObject ##################
seurat <- CreateSeuratObject(counts = counts, project = 'Moffit', meta.data=scMeta)

####################### STEP3 :Pipeline ##################
seurat <- NormalizeData(object = seurat)
seurat <- FindVariableFeatures(object = seurat, nfeatures = 2000)
seurat <- ScaleData(object = seurat)
seurat <- RunPCA(object = seurat, npcs = 50, verbose = FALSE)
seurat <- RunUMAP(object = seurat, dims = 1:50, nneighbors = 5)
saveRDS(seurat, file="Seurat/MERFISH/Moffit.rds")

