import os
os.chdir('/home/ningjin/project/SpaGE/')
import numpy as np
import pandas as pd
import scipy.stats as st
import pickle

######################## STEP 1: load information ###############
RNA_data = pd.read_csv('/home/sqsunsph/data/SpatialData/SpaGE/scRNAseq/Allen_VISp/mouse_VISp_2018-06-14_exon-matrix.csv',
                       header=0,index_col=0,sep=',')
meta_data = pd.read_csv('/home/sqsunsph/data/SpatialData/SpaGE/scRNAseq/Allen_VISp/mouse_VISp_2018-06-14_samples-columns.csv',
                        header=0,sep=',')
Genes = pd.read_csv('/home/sqsunsph/data/SpatialData/SpaGE/scRNAseq/Allen_VISp/mouse_VISp_2018-06-14_genes-rows.csv',
                        header=0,sep=',')
RNA_data.index = Genes.gene_symbol
del Genes

####################### STEP 2: filtering ###############
HighQualityCells = (meta_data['class'] != 'No Class') & (meta_data['class'] != 'Low Quality')
RNA_data = RNA_data.iloc[:,np.where(HighQualityCells)[0]]
del meta_data, HighQualityCells

Genes_count = np.sum(RNA_data, axis=1)
RNA_data = RNA_data.loc[Genes_count>10,:]
del Genes_count

####################### STEP 3: Normalization and scaling ###############
def Log_Norm(x):
    return np.log(((x/np.sum(x))*10000) + 1)

RNA_data = RNA_data.apply(Log_Norm,axis=0)
RNA_data_scaled = pd.DataFrame(data=st.zscore(RNA_data.T),index = RNA_data.columns,columns=RNA_data.index)

datadict = dict()
datadict['RNA'] = RNA_data.T
datadict['RNAScaled'] = RNA_data_scaled

with open('./AllenVISp.pkl','wb') as f:
    pickle.dump(datadict, f, protocol=4)


