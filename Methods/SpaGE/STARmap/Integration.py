import os
os.chdir('/home/ningjin/project/SpaGE/STARmap/')
import pickle
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import time as tm

######################## STEP 1: load ########################
with open ('./objects/STARmap.pkl', 'rb') as f:
	datadict = pickle.load(f)

Starmap_data = datadict['Starmap_data']
Starmap_data_scaled = datadict['Starmap_data_scaled']
labels = datadict['labels']
qhulls = datadict['qhulls']
coords = datadict['coords']
del datadict

with open ('/home/ningjin/project/SpaGE/AllenVISp.pkl', 'rb') as f:
	datadict = pickle.load(f)

RNA_data = datadict['RNA']
RNA_data_scaled = datadict['RNAScaled']
del datadict

######################## STEP 2: Leave One Out Validation ########################
Common_data = RNA_data_scaled[np.intersect1d(Starmap_data_scaled.columns,RNA_data_scaled.columns)]
df=pd.read_csv('/home/ningjin/project/SpatialMap/STARmap/objects/Random100.csv') 	
Remove_Genes = df['x']
Remove_Genes=list(Remove_Genes)


Imp_Genes = pd.DataFrame(columns=Remove_Genes)
precise_time = []
knn_time = []

for i in Remove_Genes:
	print(i)
	start = tm.time()
	from principal_vectors import PVComputation
	n_factors = 50
	n_pv = 50
	dim_reduction = 'pca'
	dim_reduction_target = 'pca'
	pv_FISH_RNA = PVComputation(n_factors = n_factors,n_pv = n_pv,dim_reduction = dim_reduction,dim_reduction_target = dim_reduction_target)
	pv_FISH_RNA.fit(Common_data.drop(i,axis=1),Starmap_data_scaled[Common_data.columns].drop(i,axis=1))
	S = pv_FISH_RNA.source_components_.T
	Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > 0.3)
	S = S[:,0:Effective_n_pv]
	Common_data_t = Common_data.drop(i,axis=1).dot(S)
	FISH_exp_t = Starmap_data_scaled[Common_data.columns].drop(i,axis=1).dot(S)
	precise_time.append(tm.time()-start)
	start = tm.time()
	nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',metric = 'cosine').fit(Common_data_t)
	distances, indices = nbrs.kneighbors(FISH_exp_t)
	Imp_Gene = np.zeros(Starmap_data.shape[0])
	for j in range(0,Starmap_data.shape[0]):
		weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
		weights = weights/(len(weights)-1)
		Imp_Gene[j] = np.sum(np.multiply(RNA_data[i][indices[j,:][distances[j,:] < 1]],weights))
	Imp_Gene[np.isnan(Imp_Gene)] = 0
	Imp_Genes[i] = Imp_Gene
	knn_time.append(tm.time()-start)


Imp_Genes.to_csv('./objects/Random100.csv')
#precise_time = pd.DataFrame(precise_time)
#knn_time = pd.DataFrame(knn_time)
#precise_time.to_csv('Results/SpaGE_PreciseTime_SVG.csv', index = False)
#knn_time.to_csv('Results/SpaGE_knnTime_SVG.csv', index = False)

