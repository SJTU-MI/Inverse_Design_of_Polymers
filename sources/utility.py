#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.random import default_rng
import pandas as pd
import torch
import math
from Cal_TC import chromosome_ids_to_smiles


df=pd.read_csv(r'smi.csv')
chromosome=df["chromosome"]

def Initial_X (var_dim,n_Random,seed):
    np.random.seed(seed)
    seed_arr=np.random.randint(0,10000,n_Random)
    np.random.seed(seed_arr[0])
    X_list=np.random.randint(0,2,var_dim)
    for seed_ in seed_arr[1:]:
        np.random.seed(seed_)
        X_list=np.vstack([X_list, np.random.randint(0,2,var_dim)])
    return X_list


# In[6]:


def cover_torch (data_array):
    if type(data_array)==torch.Tensor:
        temp_array=data_array
    else:
        temp_array=torch.tensor(data_array)        
    return temp_array
def cover_numpy (data_array):
    if type(data_array)==torch.Tensor:
        temp_array=data_array.numpy()
    else:
        temp_array=data_array
    return temp_array


# In[8]:


def Binary_num_arr(arr_):
    a0=arr_.reshape(-1,5)
    ser=[]
    for i in range(a0.shape[0]):
        num=0
        for j in range(a0.shape[1]):
            num=num+math.pow(2,j)*a0[i][::-1][j]
        ser.append (int(num))
    return ser
def X_info(chrom_X):
        chorm=np.array(Binary_num_arr(chrom_X))
        rng=default_rng(seed=12357)
        smi=chromosome_ids_to_smiles(chorm,chromosome,rng)
        return list(chorm),smi





