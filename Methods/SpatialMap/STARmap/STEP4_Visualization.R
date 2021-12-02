
rm(list = ls())
source("funcs.R")
para=c(paste0("sig",1:5),paste0("car",1:9))[index1]

########################### STEP 1: load information ###########################
ST=readRDS("ST.rds")
coordinates =readRDS("coordinates.rds")
# SpatialMap
Prediction=read.csv(paste0(para, "_ReferGene.csv"))
rownames(Prediction)=Prediction$X
Prediction=Prediction[,-1]
# Seurat
Seurat <- read.csv(paste0("Seurat/",tech, "/Prediction/ReferGene.csv"))
rownames(Seurat)=Seurat$X
Seurat=Seurat[,-1]
# SpaGE
SpaGE <- read.csv(paste0("SpaGE/",tech, "/objects/ReferGene.csv"))
rownames(SpaGE)=SpaGE$X
SpaGE=SpaGE[,-1]
########################### STEP 2: Plot ###########################
xy=paste(coordinates[,2], coordinates[,1], sep='x')
for(Gene in rownames(Prediction)){
  cat(paste0("-------- Visualizing: ",Gene,"--------"),'\n')
  truth=ST[,Gene]
  pred=c(t(Prediction[Gene,]))
  spage=as.numeric(SpaGE[,Gene])
  seurat=as.numeric(Seurat[Gene,])
  dat1<-data.frame(xy=xy, V1=log2((truth+1)*1e4))
  dat2<-data.frame(xy=xy, V1=log2((pred+1)*1e4))
  dat3<-data.frame(xy=xy, V1=log2((spage+1)*1e4))
  dat4<-data.frame(xy=xy, V1=log2((seurat+1)*1e4))
  datalist=list(dat1=dat1,dat2=dat2,dat3=dat3, dat4=dat4)
  tmplist <- datalist
  for (i in 1:length(datalist)) {
    tmplist[[i]][, 2] <- relative_func(datalist[[i]][, 2])
  }
  df <- setNames(cbind.data.frame(tmplist[[1]][, 1], 
                                  do.call(cbind, sapply(tmplist, "[", 2))), 
                 c("xy", tech, "SpatialMap","SpaGE","Seurat"))
  pp <- lapply(1:length(tmplist), function(x) {
    pattern_plot2(df, x, xy = F, main = T, titlesize = 1.5, pointsize =2)
  })
  grid.arrange(grobs = pp, ncol = 4)
}





























