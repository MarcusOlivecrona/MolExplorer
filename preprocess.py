#!/usr/bin/env python
from rdkit import Chem
from rdkit import rdBase
from rdkit import DataStructs
from rdkit.Chem import AllChem
from sklearn.decomposition import PCA
import sys
import numpy as np
import random

rdBase.DisableLog('rdApp.error')

"""Adds first three dimensions of a PCA based on ECFP6 fingerprints for the structures
   The output file will retain previous data columns"""

def fingerprints_from_mols(mols):
        fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 3) for mol in mols] 
        np_fps = []
        for fp in fps:
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, arr)
            np_fps.append(arr)
        return np_fps

def fit_pca(mols):
    fps = fingerprints_from_mols(mols)
    pca = PCA(n_components=3)
    pca.fit(fps)
    pca = pca.transform(fps)
    return pca

def read_mols(fname):
    headers = []
    mols = []
    smiles = []
    properties = []
    with open(fname, 'r') as f:
        headers = f.readline().split()
        for line in f:
            line = line.split()
            mol = Chem.MolFromSmiles(line[0])
            if mol is not None:
                mols.append(mol)
                smiles.append(line[0])
                properties.append([prop for prop in line[1:]])
    return mols, smiles, headers, properties

def write_mols_and_data(smiles, headers, properties, pca, fname):
    with open(fname, 'w') as f:
        f.write("\t".join(headers) + "\tPCA_1\tPCA_2\tPCA_3\n")
        for smiles, properties, pca in zip(smiles, properties, pca):
            f.write(smiles + "\t" + "\t".join(properties) + "\t" + str(pca[0]) + "\t" + str(pca[1]) + "\t" + str(pca[2]) + "\n")  

def main():
    in_f = sys.argv[1]
    out_f = "pca.smi"
    mols, smiles, headers, properties = read_mols(in_f)
    pca = fit_pca(mols)
    write_mols_and_data(smiles, headers, properties, pca, out_f)
        
if __name__=="__main__":
    main()

