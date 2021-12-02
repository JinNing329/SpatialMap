
rm(list = ls())
library(MASS)
library(Matrix)

para=c(paste0("sig",1:5),paste0("car",1:9))[index1]

######################## STEP 1: load information ######################## 
ST=readRDS("ST.rds")
X=readRDS("X.rds")
total=rowSums(ST)
n=nrow(ST)
sc=readRDS("scFull.rds")
K=read.table(paste0("K",para,".txt"))

filelist=dir("./Estimation")
file=filelist[grep(para, filelist)]
cat(file)
Estimation = data.frame()
for (i in 1:length(file)) {
  est=readRDS(paste0("./Estimation/",file[i]))
  Estimation=rbind(Estimation, est)
  cat(paste0("----",i,"----"),'\n')
}
## Remove NaN
index2=which(is.nan(Estimation$h2))
if(length(index2)==0){cat("No need to remove NaN")}else{Estimation=Estimation[-index2,]}
## Select by p-value
Estimation=Estimation[Estimation$pvalue<0.05,]

########################  STEP 2: Prepare reference ######################## 
TestNameList=read.csv("Random100.csv")
TestNameList=TestNameList$x
TestSC=sc[TestNameList,]
ReferNameList=setdiff(intersect(rownames(Estimation),rownames(sc)),TestNameList)
ReferSC=sc[ReferNameList,]
rm(sc)

########################  STEP 3: Prediction ######################## 
Prediction=data.frame()
for(l in 1:length(TestNameList)){
  Gene=TestNameList[l]
  cat(paste0("----- Predicting: ", Gene," -------"),'\n')
  Correlation=cor(TestSC[Gene,],t(ReferSC))
  ReferName=colnames(Correlation)[which.max(Correlation)]
  
  estimates=Estimation[ReferName,]
  alpha=estimates$alpha
  beta=estimates$beta
  h2=estimates$h2
  sigma2=estimates$sigma2
  SigmaG=h2*sigma2*K
  SigmaE=diag(rep((1-h2)*sigma2,n))
  
  g=mvrnorm(n=100, mu=rep(0,n),SigmaG)
  g=apply(g, 2, mean)
  e=mvrnorm(n=100, mu=rep(0,n),SigmaE)
  e=apply(e, 2, mean)
  
  lambda=exp(alpha+beta*X[,Gene]+g+e)
  lambda=lambda*10
  p=apply(matrix(lambda*total),1, function(x){round(mean(rpois(100,x)))})
  Prediction=rbind(Prediction,p)
}

rownames(Prediction)=TestNameList
colnames(Prediction)=paste0("sites",1:n)
####################### STEP 4: Record ######################
ofile=paste0("Prediction/",para,"_Random100.csv")
write.csv(Prediction, file=ofile)

