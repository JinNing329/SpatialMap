rm(list=ls())
library(Seurat)

################### STEP 1: load ###################
counts <- read.table(file = "mouse_VISp_2018-06-14_exon-matrix.csv",
                    row.names = 1,sep = ',', stringsAsFactors = FALSE, header = TRUE)
counts <- as.matrix(x = counts)
genes <- read.table(file = "mouse_VISp_2018-06-14_genes-rows.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
rownames(x = counts) <- make.unique(names = genes$gene_symbol)
meta.data <- read.csv(file = "mouse_VISp_2018-06-14_samples-columns.csv",
                      row.names = 1, stringsAsFactors = FALSE)

################### STEP 2: Create SeuratObject ###############
seurat <- CreateSeuratObject(counts = counts, project = 'VISp', meta.data = meta.data, min.cells = 10)
removeCells = colnames(seurat)[seurat$class %in% c('Low Quality', 'No Class')]
seurat = subset(seurat, cells = setdiff(colnames(seurat),removeCells))

################### STEP 3: Pipeline ###############
seurat <- NormalizeData(object = seurat)
seurat <- FindVariableFeatures(object = seurat, nfeatures = 2000)
seurat <- ScaleData(object = seurat)
seurat <- RunPCA(object = seurat, npcs = 50, verbose = FALSE)
seurat <- RunUMAP(object = seurat, dims = 1:50, nneighbors = 5)
seurat <- RunCCA(seurat)
saveRDS(seurat, file ="AllenVISp.rds")