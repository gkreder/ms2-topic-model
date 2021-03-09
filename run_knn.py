########################################################################
# Required packages
########################################################################
import argparse
import sys
import os
import scipy.sparse as sparse
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from pyteomics import mgf, auxiliary
from rdkit import Chem
import sklearn.metrics
########################################################################
# Parse arguments
########################################################################
parser = argparse.ArgumentParser() 
def file_choices(choices,fname, argname):
    ext = os.path.splitext(fname)[1]
    if ext not in choices:
       parser.error(f"{argname} must be one of the following filetypes ({choices})")
    return fname
parser.add_argument('--df_ids', type=lambda s:file_choices((".tsv"),s, '--df_ids'), required = True)
parser.add_argument('--df_labels', required = True)
parser.add_argument('--k', required = True, type = int)
parser.add_argument('--embedded_spectra', required = True)
parser.add_argument('--out_dir', required = True)
args = parser.parse_args()
os.system(f"mkdir -p {args.out_dir}")

with open(os.path.join(args.out_dir, 'log.txt'), 'w') as f:
    print(f"{' '.join(sys.argv)}\n", file = f)
    print(args, file = f)
########################################################################
# Load input data
########################################################################
df_ids = pd.read_csv(args.df_ids, sep = '\t')
df_labels = pd.read_csv(args.df_labels, sep = '\t', index_col = 0)
mat = sparse.load_npz(args.embedded_spectra)
test_ids = np.where(df_ids['id'].str.contains('_test'))[0]
train_ids = np.where(df_ids['id'].str.contains('_train'))[0]
########################################################################
# Calculate distances and t3 metrics
########################################################################
dists = sklearn.metrics.pairwise.cosine_similarity(mat[test_ids], mat[train_ids], dense_output=True)
df_sims = pd.DataFrame(dists, index = df_ids.iloc[test_ids]['id'], columns = df_ids.iloc[train_ids]['id'])
pred_rows = []
for i, fname in enumerate(tqdm(df_sims.index)):
    pred_rows.append(df_labels.loc[df_sims.loc[fname].sort_values(ascending = False).iloc[0 : args.k].index].mean().values)
df_preds = pd.DataFrame(pred_rows, index = df_sims.index, columns = df_labels.columns)
df_preds.to_csv(os.path.join(args.out_dir, 'df_preds.tsv'), sep = '\t')
t3s = []
for fname in tqdm(df_preds.index):
    t = df_labels.loc[fname].loc[df_preds.loc[fname].sort_values(ascending = False).iloc[0 : 3].index].values.sum()
    t3s.append(t)
t3s = np.array(t3s)
df_t3s = pd.DataFrame({'id' : df_preds.index, 't3' : t3s})
# df_t3s.to_csv(os.path.join(args.out_dir, 'df_t3s.tsv'), sep = '\t', index = False)


with open(os.path.join(args.out_dir, 'log.txt'), 'a') as f:
	print(f"\nt3>=1 - {np.sum(t3s >= 1)} / {len(test_ids)}", file = f)
	print(f"\nt3>=2 - {np.sum(t3s >= 2)} / {len(test_ids)}", file = f)






