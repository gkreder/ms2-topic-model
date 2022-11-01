########################################################################
# Required packages
########################################################################
import argparse
import sys
import os
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
import tomotopy as tp
from pyteomics import mgf, auxiliary
from rdkit import Chem
import scipy.sparse
########################################################################
# Parse arguments
########################################################################
parser = argparse.ArgumentParser() 
def file_choices(choices,fname, argname):
    ext = os.path.splitext(fname)[1]
    if ext not in choices:
       parser.error(f"{argname} must be one of the following filetypes ({choices})")
    return fname
parser.add_argument('--train_mgf', required = True)
parser.add_argument('--test_mgf', required = True)
parser.add_argument('--validate_mgf', required = True)
parser.add_argument('--df_substructs', required = True)
parser.add_argument('--embedded_spectra_out', type=lambda s:file_choices((".npz"),s, '--embedded_spectra_out'), required = True)
parser.add_argument('--df_labels_out', type=lambda s:file_choices((".tsv"),s, '--df_labels_out'), required = True)
parser.add_argument('--df_ids_out', type=lambda s:file_choices((".tsv"),s, '--df_ids_out'), required = True)
args = parser.parse_args()
########################################################################
# Functions for checking ground truth labels
########################################################################
df_substructs = pd.read_csv(args.df_substructs, sep = '\t')

def get_substruct_vec(s):
	mol = Chem.MolFromSmiles(s)
	vec = np.array([mol.HasSubstructMatch(Chem.MolFromSmarts(patt)) for patt in df_substructs['smarts']]).astype(int)
	return(vec)
def get_substruct_labels(s):
	v = get_substruct_vec(s)
	return(df_substructs.iloc[np.where(v > 0)]['index'].values)    
########################################################################
# Consolidate spectra and IDs
########################################################################
train_spectra = [x for x in mgf.MGF(args.train_mgf)]
test_spectra = [x for x in mgf.MGF(args.test_mgf)]
validate_spectra = [x for x in mgf.MGF(args.validate_mgf)]
train_ids = [s['params']['id'] + '_train' for s in train_spectra]
test_ids = [s['params']['id'] + '_test' for s in test_spectra]
validate_ids = [s['params']['id'] + '_validate' for s in validate_spectra]

ids = train_ids + test_ids + validate_ids
spectra = train_spectra + test_spectra + validate_spectra
all_smiles = [s['params']['smiles'] for s in spectra]
########################################################################
# Make ground truth label matrix
########################################################################
lab_mat = np.zeros((len(spectra), len(df_substructs)))
print('Labeling spectra using RDKit...')
for i_spec, spec in enumerate(tqdm(spectra)):
    smiles = spec['params']['smiles']
    try:
        labs = get_substruct_labels(smiles)
    except:
        labs = np.array([])
    lab_mat[i_spec, :] = df_substructs['index'].isin(labs).astype(int).values
df_labs = pd.DataFrame(lab_mat, index = ids, columns = df_substructs['index'].values)
df_labs.to_csv(args.df_labels_out, sep = '\t')
########################################################################
# Make ids - smile reference file
########################################################################
df_ids = pd.DataFrame({'id' : ids, 'smiles' : all_smiles})
df_ids.to_csv(args.df_ids_out, sep = '\t', index = False)
########################################################################
# Make embedded binarized spectra for use in KNN model
########################################################################
bins = np.arange(0.1, 1000.1, step = 0.1)
bin_names = []
for i, cbin in enumerate(bins):
    cbin = cbin.round(2)
    if i == 0:
        bn = "[0, 0.1)"
    else:
        bn = f"[{bins[i - 1].round(2)}, {cbin})"
    bin_names.append(bn)
mat = np.zeros((len(spectra), len(bins)))
print('Binarizing and embedding spectra...')
for i_spec, spec in enumerate(tqdm(spectra)):
    mza = spec['m/z array'][spec['m/z array'] < 1000]
    mzs_binned = np.digitize(mza, bins = bins)
    mat[i_spec, :][mzs_binned] = 1.0
smat = scipy.sparse.csc_matrix(mat)
scipy.sparse.save_npz(args.embedded_spectra_out, smat)


