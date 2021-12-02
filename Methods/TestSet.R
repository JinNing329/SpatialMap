# this file is used to set the test sets

rm(list=ls())
sc=readRDS("./scFull.rds")
######################## SeqFISH ######################## 
## 1. Random set 100 gene
ST=readRDS("./SeqFISH/objects/ST.rds")
Random100=sample(intersect(colnames(ST),rownames(sc)),100)
write.csv(Random100, file="./SeqFISH/objects/Random100.csv")
rm(ST)
## 2. SPARK Selected
library(SPARK)
ST=t(readRDS("./SeqFISH/objects/ST.rds"))
colnames(ST)=paste0("site",1:ncol(ST))
coordinates =readRDS("./SeqFISH/objects/coordinates.rds")
info <- cbind.data.frame(x=as.numeric(coordinates[,1]),
                         y=as.numeric(coordinates[,2]),
                         total_counts=apply(ST,2,sum))
rownames(info)=paste0("site",1:ncol(ST))
spark <- CreateSPARKObject(counts=ST,location=info[,1:2],
                           percentage = 0.1, min_total_counts = 10)
spark@lib_size <- apply(spark@counts, 2, sum)
spark <- spark.vc(spark,covariates = NULL, lib_size = spark@lib_size, 
                  num_core = 5,verbose = F)
spark <- spark.test(spark,check_positive = T,verbose = F)
SPARKSelected=rownames(spark@res_mtest)[which(spark@res_mtest$adjusted_pvalue<0.05)]
write.csv(SPARKSelected, file="./SeqFISH/objects/SPARKSelected.csv")
rm(ST)
######################## STARmap ########################
## 1. Random set 100 gene
ST=readRDS("./STARmap/objects/ST.rds")
Random100=sample(colnames(ST),100)
write.csv(Random100, file="./STARmap/objects/Random100.csv")
## 2. Gene from 
ReferGene <- c("Slc17a7", "Nov","Cux2","Rorb","Sulf2","Pcp4","Ctgf","Gad1","Pvalb","Sst","Npy","Vip")
write.csv(ReferGene, file="./STARmap/objects/ReferGene.csv")
rm(ST)

######################## OsmFISH ########################
## 1. Random set 10 gene
ST=readRDS("./osmFISH/objects/ST.rds")
Random10=sample(setdiff(rownames(sc), colnames(ST)),10)
ReferGene=c("Cux2","Tmem215","Pvrl3","Wfs1","Adam33","Rspo1","Tesc","Tox","Foxp2","Tle4")
write.csv(Random10, file="./osmFISH/objects/Random10.csv")
write.csv(ReferGene, file="./osmFISH/objects/ReferGene.csv")
rm(ST)

######################## MERFISH ########################
## Gene from 
ST=readRDS("./MERFISH/objects/ST.rds")
ReferGene=c("Gad1","Mbp","Cd24a","Myh11")
write.csv(ReferGene, file="./MERFISH/objects/ReferGene.csv")
rm(ST)
