
rm(list = ls())
source("funcs.R")
para=c(paste0("sig",1:5),paste0("car",1:9))[index1]

########################### STEP 1: load information ###########################
coordinates =readRDS("coordinates.rds")
# SpatialMap
Prediction=read.csv(paste0(para, "_ReferGene.csv"))
rownames(Prediction)=Prediction$X
Prediction=Prediction[,-1]
########################### STEP 2: Plot ###########################
xy=paste(coordinates[,1], coordinates[,2], sep='x')
for(Gene in rownames(Prediction)){
  cat(paste0("-------- Visualizing: ",Gene,"--------"),'\n')
  pred=c(t(Prediction[Gene,]))
  dat2<-data.frame(xy=xy, V1=log2((pred+1)*1e4))
  datalist=list(dat2=dat2)
  tmplist <- datalist
  for (i in 1:length(datalist)) {
    tmplist[[i]][, 2] <- relative_func(datalist[[i]][, 2])
  }
  df <- setNames(cbind.data.frame(tmplist[[1]][, 1], 
                                  do.call(cbind, sapply(tmplist, "[", 2))), 
                 c("xy", "SpatialMap"))
  pp <- lapply(1:length(tmplist), function(x) {
    pattern_plot2(df, x, xy = F, main = T, titlesize = 1.5, pointsize =3)
  })
  grid.arrange(grobs = pp, ncol = 1)
}

