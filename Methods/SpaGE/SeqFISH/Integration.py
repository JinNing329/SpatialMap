import os
os.chdir('./SpaGE/SeqFISH/')
import pickle
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import time as tm

################### STEP 1: load information ###############
with open ('SeqFISH.pkl', 'rb') as f:
	datadict = pickle.load(f)

seqFISH_data = datadict['seqFISH']
seqFISH_data_scaled = datadict['seqFISHScaled']
seqFISH_meta= datadict['seqFISHMeta']
del datadict

with open ('AllenVISp.pkl', 'rb') as f:
	datadict = pickle.load(f)

RNA_data = datadict['RNA']
RNA_data_scaled = datadict['RNAScaled']
del datadict

################### STEP 2: Leave One Out Validation ###############
Common_data = RNA_data_scaled[np.intersect1d(seqFISH_data_scaled.columns,RNA_data_scaled.columns)]
df=pd.read_csv('SPARKSelected.csv') 	
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
	
	pv_FISH_RNA = PVComputation(n_factors = n_factors,n_pv = n_pv,
		dim_reduction = dim_reduction,
		dim_reduction_target = dim_reduction_target)
    
	pv_FISH_RNA.fit(Common_data.drop(i,axis=1),seqFISH_data_scaled[Common_data.columns].drop(i,axis=1))
	S = pv_FISH_RNA.source_components_.T
	Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > 0.3)
	S = S[:,0:Effective_n_pv]
	Common_data_t = Common_data.drop(i,axis=1).dot(S)
	FISH_exp_t = seqFISH_data_scaled[Common_data.columns].drop(i,axis=1).dot(S)
	precise_time.append(tm.time()-start)
	start = tm.time()    
    
	nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',metric = 'cosine').fit(Common_data_t)
	distances, indices = nbrs.kneighbors(FISH_exp_t)
	Imp_Gene = np.zeros(seqFISH_data.shape[0])
    
	for j in range(0,seqFISH_data.shape[0]):
		weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
		weights = weights/(len(weights)-1)
		Imp_Gene[j] = np.sum(np.multiply(RNA_data[i][indices[j,:][distances[j,:] < 1]],weights))
    
	Imp_Gene[np.isnan(Imp_Gene)] = 0
	Imp_Genes[i] = Imp_Gene
	knn_time.append(tm.time()-start)

################### STEP 3: Record ###############
Imp_Genes.to_csv('SPARKSelected.csv')
