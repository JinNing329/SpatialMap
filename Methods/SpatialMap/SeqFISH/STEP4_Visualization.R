
rm(list = ls())
source("funcs.R")
para=c(paste0("sig",1:5),paste0("car",1:9))[index1]

############################# Random100 ###################

########################### STEP 1: load information ###########################
ST=readRDS("ST.rds")
coordinates =readRDS("coordinates.rds")
# SpatialMap
Prediction=read.csv(paste0("./Prediction/", para, "_SPARKSelected.csv"))
rownames(Prediction)=Prediction$X
Prediction=Prediction[,-1]
# Seurat
Seurat <- read.csv(paste0("./Seurat/",tech, "/Prediction/SPARKSelected.csv"))
rownames(Seurat)=Seurat$X
Seurat=Seurat[,-1]
# SpaGE
SpaGE <- read.csv(paste0("./SpaGE/",tech, "/objects/SPARKSelected.csv"))
rownames(SpaGE)=SpaGE$X
SpaGE=SpaGE[,-1]
########################### STEP 2: Plot ###########################
xy=paste(coordinates[,1], coordinates[,2], sep='x')
plotNames=c("Add2","Agap2","Cdk16","Cplx1","Dlgap3","Eef1a2",
            "Fscn1","Kcnab2","Mag","Mapk8ip1","Mapk8ip2","Napa",
            "R3hdm4","Shank1","Trim9")

for(Gene in plotNames){
  cat(paste0("-------- Visualizing: ",Gene,"--------"),'\n')
  truth=ST[,Gene]
  pred=c(t(Prediction[Gene,]))
  cat(paste0("-------- Correlation: ",cor(pred,truth),"--------"),'\n')
  seurat=as.numeric(Seurat[Gene,])
  spage=as.numeric(SpaGE[,Gene])
  dat1<-data.frame(xy=xy, V1=log(((truth)+1)*1e2))
  dat2<-data.frame(xy=xy, V1=log(((pred)+1)*1e2))
  dat3<-data.frame(xy=xy, V1=log(((spage)+1)*1e2))
  dat4<-data.frame(xy=xy, V1=log(((seurat)+1)*1e2))
  datalist=list(dat1=dat1,dat2=dat2,dat3=dat3, dat4=dat4)
  tmplist <- datalist
  for (i in 1:length(datalist)) {
    tmplist[[i]][, 2] <- relative_func(datalist[[i]][, 2])
  }
  df <- setNames(cbind.data.frame(tmplist[[1]][, 1], 
                                  do.call(cbind, sapply(tmplist, "[", 2))), 
                 c("xy", tech, "SpatialMap","SpaGE","Seurat"))
  pp <- lapply(1:length(tmplist), function(x) {
    pattern_plot2(df, x, xy = F, main = T, titlesize = 1.5, pointsize = 3)
  })
  grid.arrange(grobs = pp, ncol = 4)
}





























