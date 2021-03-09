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
import sklearn.metrics
import matplotlib.pyplot as plt
import seaborn as sns
########################################################################
# Parse arguments
########################################################################
parser = argparse.ArgumentParser() 
def file_choices(choices,fname, argname):
    ext = os.path.splitext(fname)[1]
    if ext not in choices:
       parser.error(f"{argname} must be one of the following filetypes ({choices})")
    return fname
parser.add_argument('--df_labels', type=lambda s:file_choices((".tsv"),s, '--df_labels'), required = True)
parser.add_argument('--df_preds_knn', type=lambda s:file_choices((".tsv"),s, '--df_preds_knn'), required = True)
parser.add_argument('--df_preds_llda', type=lambda s:file_choices((".tsv"),s, '--df_preds_llda'), required = True)
parser.add_argument('--df_substructs', type=lambda s:file_choices((".tsv"),s, '--df_substructs'), required = True)
parser.add_argument('--out_dir', required = True)
args = parser.parse_args()

os.system(f"mkdir -p {args.out_dir}")
########################################################################
# Run calculations
########################################################################
df_labels = pd.read_csv(args.df_labels, sep =  '\t', index_col = 0)
df_substructs = pd.read_csv(args.df_substructs, sep =  '\t')
df_labels_test = df_labels[df_labels.index.str.contains('_test')]
df_labels_train = df_labels[df_labels.index.str.contains('_train')]

df_preds_llda = pd.read_csv(args.df_preds_llda, sep = '\t', index_col = 0)
df_preds_llda.columns = [f"{x}_test" for x in df_preds_llda.columns]
df_preds_llda = df_preds_llda.T
df_preds_knn = pd.read_csv(args.df_preds_knn, sep = '\t', index_col = 0)

def get_rscores(df):
    rscores = []
    pscores = []
    for c in df_labels_test.columns:
        true_c = df_labels_test[c]
        if np.sum(true_c) == 0.0:
            rscore = np.nan
            pscore = np.nan
        elif c not in df.columns:
            rscore = np.nan
            pscore = np.nan
        else:
            preds_c = df[c]
            rscore = sklearn.metrics.roc_auc_score(true_c, preds_c)
            pscore = sklearn.metrics.average_precision_score(true_c, preds_c)
        rscores.append(rscore)
        pscores.append(pscore)
    sout = pd.DataFrame({'auc_roc' : rscores, 'auc_prec' : pscores}, index = df_labels_test.columns)
    return(sout)

rscores_knn = get_rscores(df_preds_knn)
rscores_knn.columns = rscores_knn.columns + '_knn'
rscores_llda = get_rscores(df_preds_llda)
rscores_llda.columns = rscores_llda.columns + '_llda'
df_rscores = pd.concat([rscores_knn, rscores_llda], axis = 1)

num_atoms = [Chem.rdchem.Mol.GetNumAtoms(Chem.MolFromSmarts(x)) for x in df_substructs['smarts']]
df_substructs['num_atoms'] = num_atoms
df_rscores['num_atoms'] = df_substructs.set_index('index').loc[df_rscores.index]['num_atoms']
df_rscores['train_appearance'] = df_labels_train[df_rscores.index].sum()
df_rscores['test_appearance'] = df_labels_test[df_rscores.index].sum()

df_rscores.to_csv(os.path.join(args.out_dir, 'df_rscores.tsv'), sep = '\t')

fig,ax = plt.subplots()
plt.scatter(df_rscores['auc_roc_knn'], df_rscores['auc_roc_llda'], c = df_labels_train.sum().loc[df_rscores.index], alpha = 0.5)
plt.plot(np.linspace(0, 1), np.linspace(0, 1), 'k--')
plt.xlabel('ROC knn')
plt.ylabel('ROC LLDA')
plt.colorbar()
plt.title(f'Mean LLDA {round(df_rscores["auc_roc_llda"].mean(), 3)}\nMean knn {round(df_rscores["auc_roc_knn"].mean(), 3)}')
plt.savefig(os.path.join(args.out_dir, 'roc_aucs.pdf'), bbox_inches = 'tight')

fig,ax = plt.subplots()
plt.scatter(df_rscores['auc_prec_knn'], df_rscores['auc_prec_llda'], c = df_labels_train.sum().loc[df_rscores.index], alpha = 0.5)
plt.plot(np.linspace(0, 1), np.linspace(0, 1), 'k--')
plt.xlabel('ROC knn')
plt.ylabel('ROC LLDA')
plt.colorbar()
plt.title(f'Mean LLDA {round(df_rscores["auc_prec_llda"].mean(), 3)}\nMean knn {round(df_rscores["auc_prec_knn"].mean(), 3)}')
plt.savefig(os.path.join(args.out_dir, 'precs.pdf'), bbox_inches = 'tight')


