
rm(list=ls())
library(Seurat)
library(methods)
library(data.table)
####################### STEP1 :load ##################
### scRNAseq:AllenVISP
colMeta<-as.data.frame(fread("mouse_VISp_2018-06-14_samples-columns.csv"))
rowMeta<- as.data.frame(fread("mouse_VISp_2018-06-14_genes-rows.csv"))
counts<- as.matrix(fread("mouse_VISp_2018-06-14_exon-matrix.csv"))
rownames(counts)=as.character(counts[,1])
counts=counts[,-1]
### remove cells with class "No Class" or "Low Quality"
index1 <- which(colMeta$class=="No Class"|colMeta$class=="Low Quality")
colMeta=colMeta[-index1,]
counts=counts[,-index1]
readDepth=rowSums(counts)
counts=counts[readDepth>10,]
rowMeta=rowMeta[readDepth>10,]
rownames(counts)=rowMeta$gene_symbol[match(rownames(counts), rowMeta$gene_entrez_id)]

### ST: osmFISH
osmFISH <- read.csv("osmFISH_data.csv")
osmFISHMeta<- read.csv("osmFISH_meta.csv")
rownames(osmFISH)=osmFISH[,1]
osmFISH=osmFISH[,-1]
osmFISH=t(osmFISH)
index1<-which(rowSums(osmFISH)<200)
osmFISH=osmFISH[-index1,]
osmFISHMeta=osmFISHMeta[-index1,]

############### match spatial with sc
scGeneName = rowMeta$gene_symbol
STGeneName = colnames(osmFISH)
ShareGene =intersect(scGeneName, STGeneName)
osmFISH=osmFISH[,match(ShareGene,STGeneName)]
counts=counts[match(ShareGene, scGeneName),]
rowMeta=rowMeta[match(ShareGene, scGeneName),]
saveRDS(osmFISH, file = "./objects/ST.rds")
####################### STEP2 :Allocating cell types to spatial sites ##################
# scRNA-seq Mean by celltype annotations
# Normalization
data = apply(counts, 1, function(x){log((x+1)*10e4)})
scale.data =apply(data, 1, scale)
# mean
matrix1 =as.data.frame(t(scale.data))
matrix1$subclass=colMeta$subclass
countsMean=aggregate(matrix1[,1:nrow(counts)], list(matrix1$subclass), mean)
rownames(countsMean)=countsMean$Group.1
countsMean=countsMean[,-1]

#####################################################
# allocating cell types to spatial sites
CorrelationMatrix = cor(x=t(countsMean),y=t(osmFISH))
Match <- apply(CorrelationMatrix, 2, function(x){which.max(x)})
rm(countsMean)
#sc <- readRDS("scFull.rds")
#data = apply(sc, 1, function(x){log((x+1)*10e4)})
#scale.data =apply(data, 1, scale)
#matrix1 =as.data.frame(t(scale.data))
#matrix1$subclass=colMeta$subclass
#countsMean=aggregate(matrix1[,1:nrow(scale.data)], list(matrix1$subclass), mean)
#rownames(countsMean)=countsMean$Group.1
#countsMean=countsMean[,-1]
#colnames(countsMean)=rownames(sc)
#saveRDS(countsMean, file="scFull_mean.rds")
countsMean=readRDS("scFull_mean.rds")
X=countsMean[Match,]
saveRDS(X, file="X.rds")
####################### STEP3 :Compute covariance matrix ##################
# kernel function
library("KRLS")
coordinates <- matrix(c(osmFISHMeta$X, osmFISHMeta$Y), ncol=2, byrow = F)
saveRDS(coordinates, file = "coordinates.rds")
coordinates <- apply(coordinates, 2, scale)

if(FALSE){
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

if(FALSE){
  D=as.matrix(dist(coordinates,p=2))
  para1=exp(seq(from=log(min(D[D>0])/2), to=log(max(D)*2), length.out = 10))[5]
  G1=gausskernel(coordinates, sigma=para1)
  W=G1
  W[G1>0.8]=1
  W[G1<=0.8]=0
  D=as.matrix(diag(rowSums(W)))
  diag(W)=0
  CARK <- list()
  for (f in 5:9) {
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
total <- apply(osmFISH, 1, sum)
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
tfile=paste0('./objects/','tfile.txt')
writeLines(STR,con=tfile,sep="\n")





