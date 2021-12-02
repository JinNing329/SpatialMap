import os
import numpy as np
import pandas as pd
import scipy.stats as st
import pickle

############# STEP 1: load information
seqFISH = pd.read_csv('cortex_svz_counts.csv',header=0)
seqFISHMeta = pd.read_csv('cortex_svz_cellcentroids.csv',header=0)

############# STEP 2: Normalization and scaling
seqFISH = seqFISH.T
cell_count = np.sum(seqFISH,axis=0)
def Log_Norm(x):
    return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

seqFISH = seqFISH.apply(Log_Norm,axis=0)
seqFISHScaled = pd.DataFrame(data=st.zscore(seqFISH.T),index = seqFISH.columns,columns=seqFISH.index)

datadict = dict()
datadict['seqFISH'] = seqFISH.T
datadict['seqFISHScaled'] = seqFISHScaled
datadict['seqFISHMeta'] = seqFISHMeta

with open('SeqFISH.pkl','wb') as f:
    pickle.dump(datadict, f)
