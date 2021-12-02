
rm(list = ls())
library(Seurat)
library(methods)
library(data.table)
library(ggplot2)
library(Matrix)

####################### STEP 1: load #######################
library(data.table)
MERFISH = fread("Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv")
MERFISH = as.data.frame(MERFISH)
MERFISH = MERFISH[MERFISH['Bregma']==0.11&MERFISH['Animal_ID']==18&MERFISH['Cell_class']!='Ambiguous',]
MERFISHMeta = MERFISH[,c(1:9)]
MERFISH = MERFISH[,c(10:170)]
drops = c('Blank_1','Blank_2','Blank_3','Blank_4','Blank_5','Fos')
MERFISH = MERFISH[ , !(colnames(MERFISH) %in% drops)]
MERFISH = ceiling(MERFISH)

####################### STEP 2: CreateSeuratObject #######################
MERFISH=t(MERFISH)
MERFISH_seurat <- CreateSeuratObject(counts = MERFISH, project = 'MERFISH', assay = 'RNA', meta.data = MERFISHMeta)

####################### STEP 3: pipeline  #######################
total.counts = colSums(x = as.matrix(MERFISH_seurat@assays$RNA@counts))
MERFISH_seurat <- NormalizeData(MERFISH_seurat, scale.factor = median(x = total.counts))
MERFISH_seurat <- ScaleData(MERFISH_seurat)

####################### STEP 4: integration #######################
Moffit <- readRDS("./Moffit.rds")
i2 <- FindTransferAnchors(reference = Moffit,query = MERFISH_seurat, features = rownames(MERFISH_seurat),
                          reduction = 'cca',reference.assay = 'RNA',query.assay = 'RNA')

refdata <- GetAssayData(object = Moffit,assay = 'RNA',slot = 'data')

imputation <- TransferData(anchorset = i2,refdata = refdata,
                           weight.reduction = 'cca')
MERFISH_seurat[['ss2']] <- imputation

####################### STEP 5:Imputation #######################
genes.leaveout =read.csv("ReferGene.csv")$x
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(MERFISH_seurat@assays$RNA)[2])
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
  imputed.ss2 <- run_imputation(ref.obj = Moffit, query.obj = MERFISH_seurat, feature.remove = genes.leaveout[[i]])
  MERFISH_seurat[['ss2']] <- imputed.ss2[, colnames(MERFISH_seurat)][['seq']]
  Imp_genes[genes.leaveout[[i]],] = as.vector(MERFISH_seurat@assays$ss2[genes.leaveout[i],])
}
write.csv(Imp_genes,file = 'Seurat/MERFISH/Prediction/ReferGene.csv')


