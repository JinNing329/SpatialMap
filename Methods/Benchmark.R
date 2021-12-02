

######################### Random 100 ############################
rm(list = ls())
library(MASS)
para=c(paste0("sig",1:5),paste0("car",1:9))

######### STEP 1: load ###########
ST=readRDS("./objects/ST.rds")
SpatialMap=list()
for (i in para) {
  cat(paste0("----------", i, "---------"),'\n')
  rfile=paste0('./Prediction/',i,"_Random100.csv")
  pred=read.csv(rfile)
  rownames(pred)=pred$X
  pred=pred[,-1]
  SpatialMap[[i]]=pred
  rm(pred)
}
Seurat <- read.csv("Prediction/Random100.csv")
rownames(Seurat)=Seurat$X
Seurat=Seurat[,-1]
# SpaGE
SpaGE <- read.csv("objects/Random100.csv")
SpaGE=SpaGE[,-1]

########## STEP 2: correlation ########
GeneSet=rownames(SpatialMap[[1]])
Truth=ST[,GeneSet]
SpatialMapCor=lapply(SpatialMap, function(x){diag(cor(t(x),Truth))})
SpaGECor=diag(cor(SpaGE, Truth))
SeuratCor=diag(cor(t(Seurat), Truth))

######### STEP 3: BoxPlot ###########
require(ggplot2)
library(reshape2)
df <- do.call(rbind.data.frame, SpatialMapCor)
colnames(df)=GeneSet
df=rbind(df, SeuratCor,SpaGECor)
para=c(paste0("Sigma ",1:5),paste0("CAR ",1:9))
df$Methods=factor(c(para,"Seurat","SpaGE"))
df <- melt(df, id.vars="Methods", measure.vars =colnames(df)[1:100])

## SpatialMap: All kernel #####
df0=df[c(grep("Sigma", df$Methods),grep("CAR", df$Methods)),]
p0 <- ggplot(df0, aes(x=Methods, y=value, fill=Methods))+geom_violin(alpha=0.2, aes(linetype=NA))+
  geom_jitter(shape=21, aes(fill=Methods), position = position_jitter(width = 0.2))+
  xlab("Gaussian kernels with different parameter")+ylab("Pearson Correlation")+
  theme_bw()+theme(legend.position = "none")
pdf("KernelRandom100.pdf", width = 7, height = 3)
p0
dev.off()



## SpatialMap: GaussianKernel #####
df1=df[grep("Sigma", df$Methods),]
p1<- ggplot(df1, aes(x=Methods, y=value, fill=Methods))+geom_violin(alpha=0.2, aes(linetype=NA))+
  geom_jitter(shape=21, aes(fill=Methods), position = position_jitter(width = 0.2))+
  xlab("Gaussian kernels with different parameter")+ylab("Pearson Correlation")+
  theme_bw()+theme(legend.position = "none")
pdf("GaussianRandom100.pdf", width = 7, height = 3)
p1
dev.off()



## SpatialMap: CAR Kernel #####
df2=df[grep("CAR", df$Methods),]
p2 <- ggplot(df2, aes(x=Methods, y=value, fill=Methods))+geom_violin(alpha=0.2, aes(linetype=NA))+
  geom_jitter(shape=21, aes(fill=Methods), position = position_jitter(width = 0.2))+
  xlab("CAR model with different parameter")+ylab("Pearson Correlation")+
  theme_bw()+theme(legend.position = "none")
pdf("CARRandom100.pdf", width = 7, height = 3)
p2
dev.off()

## Benchmark #####
index=tapply(c("Sigma 3","SpaGE","Seurat"), c("Sigma 3","SpaGE","Seurat"),function(x){grep(x, df$Methods)})
df3=df[unlist(index),]
df3$Methods=gsub(pattern = "Sigma 3", replacement = "SpatialMap", x=as.character(df3$Methods))
df3$Methods=factor(df3$Methods, levels=c("SpatialMap", "SpaGE","Seurat"), ordered = T)
p3 <- ggplot(df3, aes(x=Methods, y=value, fill=Methods))+geom_violin(alpha=0.2, aes(linetype=NA))+
  geom_jitter(shape=21, aes(fill=Methods), position = position_jitter(width = 0.2))+
  xlab("Methods")+ylab("Pearson Correlation")+
  theme_bw()+theme(legend.position = "none")
pdf("BenchmarkRandom100.pdf", width = 7, height = 3)
p3
dev.off()









############################## ReferGene ##########################
rm(list = ls())
library(MASS)
para=c(paste0("sig",1:5), paste0("car",1:2))

######### STEP 1: load ###########
ST=readRDS("./objects/ST.rds")
SpatialMap=list()
for (i in para) {
  cat(paste0("----------", i, "---------"),'\n')
  rfile=paste0('./Prediction/',i,"_ReferGene.csv")
  pred=read.csv(rfile)
  rownames(pred)=pred$X
  pred=pred[,-1]
  SpatialMap[[i]]=pred
  rm(pred)
}
Seurat <- read.csv("./Prediction/ReferGene.csv")
rownames(Seurat)=Seurat$X
Seurat=Seurat[,-1]
# SpaGE
SpaGE <- read.csv("./objects/ReferGene.csv")
SpaGE=SpaGE[,-1]

########## STEP 2: correlation ########
GeneSet=rownames(SpatialMap[[1]])
Truth=ST[,GeneSet]
SpatialMapCor=lapply(SpatialMap, function(x){diag(cor(t(x),Truth))})
SpaGECor=diag(cor(SpaGE, Truth))
SeuratCor=diag(cor(t(Seurat), Truth))

######### STEP 3: VlnPlot ###########
library(ggplot2)
library(reshape2)
df <- do.call(rbind.data.frame, SpatialMapCor)
colnames(df)=GeneSet
df=rbind(df, SeuratCor,SpaGECor)
para=c(paste0("Sigma ",1:5), paste0("CAR ", 1:2))
df$Methods=factor(c(para,"Seurat","SpaGE"))
df <- melt(df, id.vars="Methods", measure.vars =colnames(df)[-ncol(df)])

## SpatialMap: GaussianKernel #####
df1=df[grep("Sigma", df$Methods),]
pdf("./Visualize/GaussianSPARK.pdf", height = 3, width = 7)
ggplot(df1, aes(x=Methods, y=value, fill=Methods))+geom_violin(alpha=0.2, aes(linetype=NA))+
  geom_jitter(shape=21, aes(fill=Methods), position = position_jitter(width = 0.2))+
  xlab("Gaussian kernels with different parameter")+ylab("Pearson Correlation")+
  theme_bw()+theme(legend.position = "none")
dev.off()
## SpatialMap: CAR Kernel #####
df2=df[grep("CAR", df$Methods),]
pdf("./Visualize/CARSPARK.pdf", height = 3, width = 7)
ggplot(df2, aes(x=Methods, y=value, fill=Methods))+geom_violin(alpha=0.2, aes(linetype=NA))+
  geom_jitter(shape=21, aes(fill=Methods), position = position_jitter(width = 0.2))+
  xlab("CAR model with different parameter")+ylab("Pearson Correlation")+
  theme_bw()+theme(legend.position = "none")
dev.off()
## Benchmark #####
index=tapply(c("Sigma 2","SpaGE","Seurat"), c("Sigma 2","SpaGE","Seurat"),function(x){grep(x, df$Methods)})
df3=df[unlist(index),]
df3$Methods=gsub(pattern = "Sigma 2", replacement = "SpatialMap", x=as.character(df3$Methods))
df3$Methods=factor(df3$Methods, levels=c("SpatialMap", "SpaGE","Seurat"), ordered = T)
pdf("./Visualize/BenchmarkSPARK.pdf", height = 3, width = 7)
ggplot(df3, aes(x=Methods, y=value, fill=Methods))+geom_violin(alpha=0.2, aes(linetype=NA))+
  geom_jitter(shape=21, aes(fill=Methods), position = position_jitter(width = 0.2))+
  xlab("Methods")+ylab("Pearson Correlation")+
  theme_bw()+theme(legend.position = "none")
dev.off()

