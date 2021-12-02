rm(list=ls())
library(Seurat)
library(methods)
library(data.table)
library(ggplot2)
####################### STEP 1: load #######################
SeqFISH <- fread("cortex_svz_counts.csv")
SeqFISHMeta <- fread("cortex_svz_cellcentroids.csv")
SeqFISH=as.data.frame(SeqFISH)
SeqFISHMeta=as.data.frame(SeqFISHMeta)

####################### STEP 2: filtering #######################
SeqFISH=t(SeqFISH)
colnames(SeqFISH)=paste0("site",1:ncol(SeqFISH))
rownames(SeqFISHMeta)=paste0("site",1:ncol(SeqFISH))
seqFISH <-CreateSeuratObject(counts=SeqFISH, project='seqFISH', assay='RNA', meta.data=SeqFISHMeta)
rm(SeqFISH)
####################### STEP 3: pipeline #######################
total.counts = colSums(x = as.matrix(seqFISH@assays$RNA@counts))
seqFISH <- NormalizeData(seqFISH, scale.factor = median(x = total.counts))
seqFISH <- FindVariableFeatures(object = seqFISH, nfeatures = 2000)
seqFISH <- ScaleData(seqFISH)
seqFISH <- RunPCA(object = seqFISH, npcs = 50, verbose = FALSE)
####################### STEP 4: integration #######################
allen <- readRDS("AllenVISp.rds")
#Project on allen labels
i2 <- FindTransferAnchors(reference = allen,query = seqFISH, 
                          features = rownames(seqFISH), reduction = 'cca',
                          reference.assay = 'RNA',query.assay = 'RNA')

refdata <- GetAssayData(object = allen,assay = 'RNA',slot = 'data')

imputation <- TransferData(anchorset = i2,refdata = refdata,weight.reduction = 'cca')
seqFISH[['ss2']] <- imputation

####################### STEP 5:Imputation #######################
genes.leaveout =read.csv("Random100.csv")$x
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(seqFISH@assays$RNA)[2])
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
  imputed.ss2 <- run_imputation(ref.obj = allen, query.obj = seqFISH, feature.remove = genes.leaveout[[i]])
  seqFISH[['ss2']] <- imputed.ss2[, colnames(seqFISH)][['seq']]
  Imp_genes[genes.leaveout[[i]],] = as.vector(seqFISH@assays$ss2[genes.leaveout[i],])
  write.csv(Imp_genes,file = 'SeqFISH/Prediction/Random100.csv')
}
