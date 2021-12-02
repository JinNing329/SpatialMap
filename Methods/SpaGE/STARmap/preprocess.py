import os
os.chdir('./SpaGE/STARmap')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle
from viz import GetQHulls

###################### STEP 1: load ####################
counts = np.load('cell_barcode_count.npy')
Genes = pd.read_csv('genes.csv',header=None)
Genes = (Genes.iloc[:,0])
counts = pd.DataFrame(data=counts,columns=Genes)
Starmap_data = counts.T
del Genes, counts

###################### STEP 2: normalization and scaling ####################
cell_count = np.sum(Starmap_data,axis=0)
def Log_Norm(x):
    return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

Starmap_data = Starmap_data.apply(Log_Norm,axis=0)
Starmap_data_scaled = pd.DataFrame(data=st.zscore(Starmap_data.T),index = Starmap_data.columns,columns=Starmap_data.index)

###################### STEP 3: record ####################
labels = np.load('labels.npz')["labels"]
qhulls,coords = GetQHulls(labels)

datadict = dict()
datadict['Starmap_data'] = Starmap_data.T
datadict['Starmap_data_scaled'] = Starmap_data_scaled
datadict['labels'] = labels
datadict['qhulls'] = qhulls
datadict['coords'] = coords

with open('STARmap.pkl','wb') as f:
    pickle.dump(datadict, f)


