
rm(list=ls())
library(Seurat)
library(xlsx)

####################### STEP1 :preprocessing ##################

### scRNAseq: Moffit
counts <- Read10X("GSE113576/")
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
## filtering gene
counts=counts[which(rowSums(counts)>100),]
#saveRDS(counts, "scFull.rds")

### ST: MERFISH
library(data.table)
MERFISH = fread("Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv")
MERFISH = as.data.frame(MERFISH)
MERFISH = MERFISH[MERFISH['Bregma']==0.11&MERFISH['Animal_ID']==18&MERFISH['Cell_class']!='Ambiguous',]

MERFISHMeta = MERFISH[,c(1:9)]
MERFISH = MERFISH[,c(10:170)]
drops = c('Blank_1','Blank_2','Blank_3','Blank_4','Blank_5','Fos')
MERFISH = MERFISH[ , !(colnames(MERFISH) %in% drops)]
MERFISH = ceiling(MERFISH)


############### match spatial with sc
scGeneName = rownames(counts)
STGeneName = colnames(MERFISH)
ShareGene =intersect(scGeneName, STGeneName)
MERFISH=MERFISH[,match(ShareGene,STGeneName)]
counts=counts[match(ShareGene, scGeneName),]
saveRDS(counts, file="sc.rds")
saveRDS(MERFISH, file = "ST.rds")
saveRDS(MERFISHMeta, file="Meta.rds")
####################### STEP2 :Allocating cell types to spatial sites ##################
# scRNA-seq Mean by celltype annotations
# Normalization
data = apply(counts, 1, function(x){log((x+1)*10e4)})
scale.data = apply(data, 1, scale)
# mean
matrix1 =as.data.frame(t(scale.data))
matrix1$subclass=scMeta$Class
countsMean=aggregate(matrix1[,1:nrow(counts)], list(matrix1$subclass), mean)
rownames(countsMean)=countsMean$Group.1
countsMean=countsMean[,-1]

#####################################################
# allocating cell types to spatial sites
CorrelationMatrix = cor(x=t(countsMean),y=t(MERFISH))
Match <- apply(CorrelationMatrix, 2, function(x){which.max(x)})
rm(countsMean)
#sc <- readRDS("scFull.rds")
#data = apply(sc, 1, function(x){log((x+1)*10e4)})
#scale.data =apply(data, 1, scale)
#matrix1 =as.data.frame(t(scale.data))
#matrix1$subclass=scMeta$Class
#countsMean=aggregate(matrix1[,1:nrow(scale.data)], list(matrix1$subclass), mean)
#rownames(countsMean)=countsMean$Group.1
#countsMean=countsMean[,-1]
#colnames(countsMean)=rownames(sc)
#saveRDS(countsMean, file="scFull_mean.rds")
countsMean <-readRDS("scFull_mean.rds")
X=countsMean[Match,]
saveRDS(X, file="X.rds")

####################### STEP3 :Compute covariance matrix ##################
# kernel function
library("KRLS")
coordinates <- matrix(c(MERFISHMeta$Centroid_X, MERFISHMeta$Centroid_Y), ncol=2, byrow = F)
saveRDS(coordinates, file = "coordinates.rds")
coordinates <- apply(coordinates, 2, scale)
if(TRUE){
  # Gaussian Kernel
  D=as.matrix(dist(coordinates, p=2))
  para=exp(seq(from=log(min(D[D>0])/2), to=log(max(D)*2), length.out = 10))[3:7]
  for (f in 1:length(para)) {
    para1=para[f]
    K=gausskernel(coordinates, sigma=para1)
    # write in format
    STR=character()
    n=nrow(K)
    for(i in 1:n)
    {
      str=character()
      for(j in 1:n)
      {
        str=paste(str,as.character(K[i,j]),sep='\t')
      }
      STR=paste(STR,str,sep='\n')
    }
    STR=sub('^..','',STR)
    Kfile=paste0("./objects/Ksig",f,".txt")
    writeLines(STR,con=Kfile,sep="\n")
    cat(f,'\n')
  }
}

if(TRUE){
  D=as.matrix(dist(coordinates,p=2))
  para1=exp(seq(from=log(min(D[D>0])/2), to=log(max(D)*2), length.out = 10))[5]
  G1=gausskernel(coordinates, sigma=para1)
  W=G1
  W[G1>=0.9]=1
  W[G1<0.9]=0
  D=as.matrix(diag(rowSums(W)))
  diag(W)=0
  CARK <- list()
  for (f in 1:5) {
    alpha=f*0.1
    K=solve(D-alpha*W)
    K=K/sum(diag(K))*nrow(K)
    CARK[[f]]=K
    
    # write in format
    STR=character()
    n=nrow(K)
    for(i in 1:n)
    {
      str=character()
      for(j in 1:n)
      {
        str=paste(str,as.character(K[i,j]),sep='\t')
      }
      STR=paste(STR,str,sep='\n')
    }
    STR=sub('^..','',STR)
    Kfile=paste0("./objects/Kcar",f,".txt")
    writeLines(STR,con=Kfile,sep="\n")
    cat(f,'\n')
  }
}


######################## STEP 4: Write file in required format
total <- apply(MERFISH, 1, sum)
n<-length(total) ## n is the number of spots
str='site'
for(i in 1:n)
{
  str1=paste('idv',as.character(i),sep='')
  str=paste(str,str1,sep='\t')
}
STR=str
str='site1'
for(i in 1:n)
{
  str=paste(str,as.character(total[i]),sep='\t')
}
STR=paste(STR,str,sep='\n')
tfile=paste0('tfile.txt')
writeLines(STR,con=tfile,sep="\n")





