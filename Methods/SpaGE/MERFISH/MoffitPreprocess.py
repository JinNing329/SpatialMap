import os
os.chdir('Moffit_RNA/')
import numpy as np
import pandas as pd
import scipy.stats as st
import pickle
import scipy.io as io

############################## STEP 1: load ############################
genes = pd.read_csv('GSE113576/genes.tsv',sep='\t',header=None)
barcodes = pd.read_csv('GSE113576/barcodes.tsv',sep='\t',header=None)

genes = np.array(genes.loc[:,1])
barcodes = np.array(barcodes.loc[:,0])
RNA_data = io.mmread('GSE113576/matrix.mtx')
RNA_data = RNA_data.todense()
RNA_data = pd.DataFrame(RNA_data,index=genes,columns=barcodes)

############################## STEP 2: filtering ############################
Genes_count = np.sum(RNA_data > 0, axis=1)
RNA_data = RNA_data.loc[Genes_count >=100,:]
del Genes_count
############################## STEP 3: normalization and scaling ############################
def Log_Norm(x):
    return np.log(((x/np.sum(x))*10000) + 1)

RNA_data = RNA_data.apply(Log_Norm,axis=0)
RNA_data_scaled = pd.DataFrame(data=st.zscore(RNA_data.T),index = RNA_data.columns,columns=RNA_data.index)
############################## STEP 4:record ############################
datadict = dict()
datadict['RNA_data'] = RNA_data.T
datadict['RNA_data_scaled'] = RNA_data_scaled

with open('Moffit.pkl','wb') as f:
    pickle.dump(datadict, f, protocol=4)
