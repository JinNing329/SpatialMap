
rm(list=ls())
library(doParallel)
library(Matrix)
library(Rcpp)
source("PQLseq.R")
Rcpp::sourceCpp("PQLseq.cpp")

sig=c(paste0("sig",1:5),paste0("car",1:9))[index1]

cat(paste0("--------------",sig,"-------------",'\n'))
  ######################### read data #################
  y <- readRDS("./ST.rds")[,((index2-1)*11+1):(index2*11)]
  y <-t(y)
  genotype <- readRDS("./X.rds")[,((index2-1)*11+1):(index2*11)]
  genotype <-t(genotype)
  genename=rownames(y)
  m=nrow(y)
  n=ncol(y)
  #########  set format
  str="idv1"
  for (i in 2:n) {
    str1=paste("idv",as.character(i), sep="")
    str=c(str, str1)
  }
  col_name=str
  row_name=rownames(y)
  ###### write total file
  tfile="tfile.txt"
  total=read.table(tfile, header = T)
  total=total[,-1]
  total=as.matrix(total)
  ###### write K file 
  Kfile=paste0("K", sig,".txt")
  K <- read.table(Kfile)
  
  result=data.frame()
  for (i in 1:m) {
    ###### write count file
    count=t(y[i,])
    colnames(count)=col_name
    rownames(count)=row_name[i]
    count=as.data.frame(count)
    ###### write predictor file
    predictor=genotype[i,]
    
    names(predictor)=NULL
    res <- pqlseq(RawCountDataSet=count, Phenotypes=predictor,
                      RelatednessMatrix=K, LibSize=total, filtering=F,
                      fit.tol=1e-5,fit.model="PMM",numCore=2)
    result=rbind(result, res)
  }
  path <- paste0(sig,"_",index2,".rds")
  saveRDS(result,path) 

