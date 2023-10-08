#!/usr/bin/env python
# coding: utf-8

# In[47]:


from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
import math
from numpy.random import default_rng
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})


# In[31]:


import tensorflow as tf
from tensorflow.python.ops import math_ops
from tensorflow.python.framework import ops
tf.keras.backend.set_floatx('float64')
from tensorflow.keras.models import Sequential, save_model, load_model
from tensorflow.keras.layers import Dense, Activation, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow import keras


# In[32]:


def longest_smiles(smiles):
    """Returns longest chain of polymer with more than two stars.
        
    Done so only a linear chain is returned.
    """
    smiles = smiles.replace('*', 'Bi')
    mol_temp = Chem.MolFromSmiles(smiles)
    # Find endatom index
    at_idx = []  # Indices of Bi atoms in this fragment
    try:
        num_atom_m = mol_temp.GetNumAtoms()  # Number of atoms in this fragment
    except:
        return('') # This smiles is not wrong in syntax
    atom_i = 0
    for atom in mol_temp.GetAtoms():
        atom_i = atom_i + 1
        if atom.GetSymbol() == 'Bi':
            at_idx.append(atom.GetIdx())
        if atom_i == num_atom_m:
            break

    # Make all combination of 2
    star_pairs = list(itertools.combinations(at_idx, 2))
    # Calculate length
    idx1 = 0
    idx2 = 0
    len_path_max = -999
    for star_pair in star_pairs:
        path_temp = Chem.GetShortestPath(mol_temp, star_pair[0], star_pair[1])
        path_temp_ = list(path_temp)

        len_path = len(path_temp_)
        if len_path > len_path_max:
            len_path_max = len_path
            idx1 = star_pair[0]
            idx2 = star_pair[1]

    for idx in at_idx:
        if idx == idx1 or idx == idx2:
            mol_temp.GetAtomWithIdx(idx).SetAtomicNum(82)

    smiles = Chem.MolToSmiles(mol_temp)

    smiles = smiles.replace('([Bi])', '')
    smiles = smiles.replace('[Bi]', '')
    smiles = smiles.replace('[Pb]', '[*]')

    return smiles

def chromosome_ids_to_smiles(chromosome_ids: list, chromosomes: dict,
        rng: default_rng, **kwargs) -> str:
    """Combined chromosome ids to create new polymer smiles
       
    Combines chromosomes sequentially according to index, but combines
    joints randomly.

    Args:
        chromosome_ids (list): 
            list of chromosome ids
        chromosomes (dict):
            dict with key of chromosome id and value of the chromosome
        rng (default_rng):
            random number generator
            

    Returns (str):
        smiles string of combined chromosomes
    """
    chromosomes = [chromosomes[chromosome_id] for 
                   chromosome_id in chromosome_ids]
    mols_of_chromosomes = list()
    for chromosome in chromosomes:
        m = Chem.MolFromSmiles(chromosome)
        if m is not None:
            mols_of_chromosomes.append(m)
        # Error occurs
        else:
            return None

    # Get Ids of 'Bi' atoms in each fragment
    # edtg stands for symbols used in ladder polymers
    mols_edtg = list()
    mols_at_idx = list()

    # Iterate for each fragment (m) in the list ms
    for m, i in zip(mols_of_chromosomes, range(len(mols_of_chromosomes))):
        at_idx = []  # Indices of Bi atoms in this fragment
        num_atom_m = m.GetNumAtoms()  # Number of atoms in this fragment
        atom_i = 0
        for atom in m.GetAtoms():
            atom_i = atom_i + 1
            if atom.GetSymbol() == 'Bi':
                at_idx.append(atom.GetIdx())
            if atom_i == num_atom_m:
                break
        edtg = [''] * len(at_idx)
        temp_list_idx = list(range(len(at_idx)))

        temp_rnd_e = rng.choice(temp_list_idx)
        edtg[temp_rnd_e] = 'e'
        temp_list_idx.remove(temp_rnd_e)

        # Select random index and set the element of the index as 't'
        temp_rnd_t = rng.choice(temp_list_idx)
        edtg[temp_rnd_t] = 't'
        temp_list_idx.remove(temp_rnd_t)

        mols_at_idx.append(at_idx)
        mols_edtg.append(edtg)

    # Convert edtg to endatom symbols
    endatom = {
               'e': ['Sb', 51], 't': ['Po', 84], 
              }

    for m, i in zip(mols_of_chromosomes, range(len(mols_of_chromosomes))):
        for key in endatom.keys():
            try:
                m.GetAtomWithIdx(
                    mols_at_idx[i][
                        mols_edtg[i].index(key)
                                ]
                                ).SetAtomicNum(endatom[key][1])
            except ValueError as e:
                pass

    # Time to connect fragments.

    L_mol = mols_of_chromosomes[0]
    L_n_atoms = mols_of_chromosomes[0].GetNumAtoms()
    L_ms_edtg = mols_edtg[0]
    L_ms_at_idx = mols_at_idx[0]

    for i in range(1, len(mols_of_chromosomes)):
        R_mol = mols_of_chromosomes[i]
        R_ms_edtg = mols_edtg[i]
        R_ms_at_idx = mols_at_idx[i]

        # if 't' in L_ms_edtg and 'e' in R_ms_edtg:
        L_index_of_t = L_ms_at_idx[L_ms_edtg.index('t')]
        R_index_of_e = R_ms_at_idx[R_ms_edtg.index('e')]

        combo = Chem.CombineMols(L_mol, R_mol)
        edcombo = Chem.EditableMol(combo)

        edcombo.AddBond(L_index_of_t, R_index_of_e + L_n_atoms, order=Chem.rdchem.BondType.SINGLE)

        L_mol = edcombo.GetMol()
        L_ms_edtg = [''] * len(L_ms_edtg)
        L_ms_edtg.extend(mols_edtg[i])
        L_ms_at_idx.extend([x + L_n_atoms for x in mols_at_idx[i]])
        L_n_atoms = L_mol.GetNumAtoms()

    SMILES_connected = Chem.MolToSmiles(L_mol)
    SMILES_connected = SMILES_connected.replace("[Po][Sb]", "")
    SMILES_connected = SMILES_connected.replace("[Sb][Po]", "")
    SMILES_connected = SMILES_connected.replace("[Po]", "[*]")
    SMILES_connected = SMILES_connected.replace("[Sb]", "[*]")
    SMILES_connected = SMILES_connected.replace("[Bi]", "[*]")

    if SMILES_connected.count('*') > 2:
        SMILES_connected = longest_smiles(SMILES_connected)

    return SMILES_connected


# In[33]:


def Binary_num_arr(arr_):
    a0=arr_.reshape(-1,5)
    ser=[]
    for i in range(a0.shape[0]):
        num=0
        for j in range(a0.shape[1]):
            num=num+math.pow(2,j)*a0[i][::-1][j]
        ser.append (int(num))
    return ser


# In[34]:


def pre_TC(chrom_X,model,selected_keys,chromosome,Smi=True):
    chorm=np.array(Binary_num_arr(chrom_X))
    rng=default_rng(seed=12357)
    smi=chromosome_ids_to_smiles(chorm,chromosome,rng)
    new_length=len(selected_keys)
    mol=Chem.MolFromSmiles(smi)
    fp=AllChem.GetMorganFingerprint(mol, radius=3)
    fp_n=fp.GetNonzeroElements()
    my_finger = [0] * new_length
    for key in fp_n.keys():
        if key in selected_keys:
            index = np.where(selected_keys == key)[0][0]
            my_finger[index] = fp_n[key]
    if Smi:
        return np.round(model.predict(pd.DataFrame(my_finger).T, verbose=0)[0][0],3),chorm,smi
    else:
        return np.round(model.predict(pd.DataFrame(my_finger).T, verbose=0)[0][0],3)


# In[35]:


x_arr=np.random.randint(0,2,25)


# In[42]:


TC_model=load_model(r"TC_DNN.h5",compile=False)


# In[45]:


df=pd.read_csv(r'smi.csv')
chromosome=df["chromosome"]
selected_cor=pd.read_csv(r'select_keys.csv',index_col=0)
selected_keys=np.array(selected_cor["key"])
selected_Corr_df=np.array(selected_cor["Index"])


# In[48]:


pre_TC(x_arr,TC_model,selected_keys,chromosome,Smi=True)

