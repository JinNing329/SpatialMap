
rm(list = ls())
library(data.table)
library(Seurat)
library(methods)
library(data.table)
library(ggplot2)
library(Matrix)
####################### STEP 1: load #######################
STARmap <-read.csv("cell_barcode_count.csv", header = F)
gene.names <- read.csv("genes.csv",header = F)
coordinates <- read.table("centroids.tsv", sep = '\t')
colnames(STARmap) = gene.names$V1
STARmap=as.matrix(STARmap)

####################### STEP 2: SeuratObject #######################
STARmap=t(STARmap)
colnames(STARmap)=paste0("site",1:ncol(STARmap))
rownames(coordinates)=paste0("site",1:ncol(STARmap))
colnames(coordinates)=c("x","y")
starmap <- CreateSeuratObject(counts = STARmap, project = "STARmap", assay = "RNA", 
                              meta.data = coordinates)
####################### STEP 3: pipeline #######################
total.counts = colSums(x = as.matrix(starmap@assays$RNA@counts))
starmap <- NormalizeData(starmap, scale.factor = median(x = total.counts))
starmap <- FindVariableFeatures(object = starmap, nfeatures = 2000)
starmap <- ScaleData(starmap)
starmap <- RunPCA(object = starmap, npcs = 50, verbose = FALSE)

####################### STEP 4: integration #######################
allen <- readRDS("AllenVISp.rds")
#Project on allen labels
i2 <- FindTransferAnchors(reference = allen,query = starmap, features = rownames(starmap),
                          reduction = 'cca',reference.assay = 'RNA',query.assay = 'RNA')

refdata <- GetAssayData(object = allen,assay = 'RNA',slot = 'data')

imputation <- TransferData(anchorset = i2,refdata = refdata,
                           weight.reduction = 'cca')
starmap[['ss2']] <- imputation

####################### STEP 5:Imputation #######################
genes.leaveout =read.csv("Random100.csv")
genes.leaveout=genes.leaveout$x
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(starmap@assays$RNA)[2])
rownames(Imp_genes) <- genes.leaveout
anchor_time <- vector(mode= "numeric")
Transfer_time <- vector(mode= "numeric")

run_imputation <- function(ref.obj, query.obj, feature.remove) {
  message(paste0('removing ', feature.remove))
  features <- setdiff(rownames(query.obj), feature.remove)
  DefaultAssay(ref.obj) <- 'RNA'
  DefaultAssay(query.obj) <- 'RNA'
  
  start_time <- Sys.time()
  anchors <- FindTransferAnchors(
    reference = ref.obj,
    query = query.obj,
    features = features,
    dims = 1:30,
    reduction = 'cca'
  )
  end_time <- Sys.time()
  anchor_time <<- c(anchor_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  
  refdata <- GetAssayData(
    object = ref.obj,
    assay = 'RNA',
    slot = 'data'
  )
  
  start_time <- Sys.time()
  imputation <- TransferData(
    anchorset = anchors,
    refdata = refdata,
    weight.reduction = 'cca'
  )
  query.obj[['seq']] <- imputation
  end_time <- Sys.time()
  Transfer_time <<- c(Transfer_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  return(query.obj)
}

for(i in 1:length(genes.leaveout)) {
  cat(paste(i,genes.leaveout[i]),'\n')
  imputed.ss2 <- run_imputation(ref.obj = allen, query.obj = starmap, feature.remove = genes.leaveout[[i]])
  starmap[['ss2']] <- imputed.ss2[, colnames(starmap)][['seq']]
  Imp_genes[genes.leaveout[[i]],] = as.vector(starmap@assays$ss2[genes.leaveout[i],])
}
write.csv(Imp_genes,file = 'STARmap/Prediction/Random100.csv')
