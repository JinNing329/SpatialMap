
rm(list=ls())
library(Seurat)
######################### STEP 1: load ######################
osmFISH <- read.csv("osmFISH_data.csv")
osmFISHMeta<- read.csv("osmFISH_meta.csv")
rownames(osmFISH)=osmFISH[,1]
osmFISH=osmFISH[,-1]
osmFISH=t(osmFISH)
index1<-which(rowSums(osmFISH)<200)
osmFISH=osmFISH[-index1,]
osmFISHMeta=osmFISHMeta[-index1,]
rownames(osmFISH)=paste0("sites", 1:nrow(osmFISH))
rownames(osmFISHMeta)=paste0("sites", 1:nrow(osmFISHMeta))
seurat <- CreateSeuratObject(counts=t(osmFISH), project = "osmFISH",assay = "FISH",meta.data = osmFISHMeta)

## original clustering visualization
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat,selection.method = "vst")
plot1 <- VariableFeaturePlot(seurat)
LabelPoints(plot = plot1,points = rownames(seurat),repel = TRUE)
seurat <- ScaleData(seurat)
seurat <-RunPCA(seurat)
seurat <- RunUMAP(seurat, dims=1:10)
seurat <- FindNeighbors(seurat, dims=1:10)
seurat <- FindClusters(seurat, resolution=0.7)

p1 <- DimPlot(seurat, reduction = "umap", pt.size = 2)
p2 <-DimPlot(seurat, reduction = "umap", group.by='Region', pt.size = 2)
pdf("./Visualize/Original.pdf", width = 14, height = 7)
p1+p2
dev.off()

## Advanced clustering visualization
para=c(paste0("sig",1:5),paste0("car",1:9))[paraindex]
prediction=read.csv(paste0("./Prediction/",para,"_ReferGene.csv"))
rownames(prediction)=prediction$X
prediction=prediction[,-1]

newcounts=rbind(prediction, t(osmFISH))
seuratNew <- CreateSeuratObject(counts=newcounts, project = "osmFISH",assay = "FISH",meta.data = osmFISHMeta)
seuratNew <- NormalizeData(seuratNew)
seuratNew <- FindVariableFeatures(seuratNew,selection.method = "vst")
plot3 <- VariableFeaturePlot(seuratNew)
LabelPoints(plot = plot3,points = rownames(seuratNew),repel = TRUE)
seuratNew <- ScaleData(seuratNew)
seuratNew <-RunPCA(seuratNew)
seuratNew <- RunUMAP(seuratNew, dims=1:5)
seuratNew <- FindNeighbors(seuratNew, dims=1:5)
seuratNew <- FindClusters(seuratNew, resolution=0.2)

p1 <- DimPlot(seuratNew, reduction = "umap", pt.size = 2)
p2 <-DimPlot(seuratNew, reduction = "umap", group.by='Region', pt.size = 2)
pdf(paste0("./Visualize/",para,".pdf"), width = 14, height = 7)
p1+p2
dev.off()

