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
import scipy.spatial
########################################################################
# Parse arguments
########################################################################
parser = argparse.ArgumentParser() 
def file_choices(choices,fname, argname):
	ext = os.path.splitext(fname)[1]
	if ext not in choices:
	   parser.error(f"{argname} must be one of the following filetypes ({choices})")
	return fname
parser.add_argument('--Q', required = True, type = float)
parser.add_argument('--B', required = True, type = float)
parser.add_argument('--out_dir', required = True)
parser.add_argument('--train_mgf', required = True)
parser.add_argument('--test_mgf', required = True)
parser.add_argument('--documents_dir', required = True)
parser.add_argument('--df_substructs', type=lambda s:file_choices((".tsv"),s, '--df_substructs'), required = True)
parser.add_argument('--df_labels', type=lambda s:file_choices((".tsv"),s, '--df_labels'), required = True)
parser.add_argument('--num_iterations', required = True)
parser.add_argument('--mz_cutoff', type = float, default = 30.0)
parser.add_argument('--loss_types', choices = ['none', 'parent', 'all'], type = str, default = 'all')
args = parser.parse_args()
os.system(f"mkdir -p {args.out_dir}")
with open(os.path.join(args.out_dir, 'log.txt'), 'w') as f:
	print(f"{' '.join(sys.argv)}\n", file = f)
	print(args, file = f)
########################################################################
df_substructs = pd.read_csv(args.df_substructs, sep = '\t')

def get_mes_vec(s):
	mol = Chem.MolFromSmiles(s)
	vec = np.array([mol.HasSubstructMatch(Chem.MolFromSmarts(patt)) for patt in df_substructs['smarts']]).astype(int)
	return(vec)
def get_mes_labels(s):
	v = get_mes_vec(s)
	return(df_substructs.iloc[np.where(v > 0)]['index'].values)    
########################################################################
def glf(Q, B, x, v = 0.5, A = 0, K = 100, C = 1):
	res = (K - A) / np.power(C + Q * np.exp(-B * x), (1 / v))
	return(res)
########################################################################
doc_indices = []
doc_fnames = []
doc_smiles = []
doc_sids = []
mdl = tp.LLDAModel(seed = 2010)
print('Generating documents...')
sys.stdout.flush()

def make_doc(spec, documents_dir):
	in_fname = os.path.join(documents_dir, spec['params']['id'] + '.tsv')
	df = pd.read_csv(in_fname, sep = '\t')
	if 'mz' not in df.columns:
		df['mz'] = df['m/z']
	df['mz_rounded'] = df['mz'].round(2)
	df['rel_intensity'] = df['intensity'] / df['intensity'].max() * 100
	df['rel_intensity_rounded'] = df['rel_intensity'].round(0).astype(int)
	df = df[df['mz'] >= args.mz_cutoff] # filter by m/z
	df['intensity'] = df['intensity'].astype(int)
	df['glf_intensity'] = glf(Q = args.Q, B = args.B, x = df['rel_intensity'])

	doc_frags = np.concatenate([[w for _ in range(round(c))] for w, c in df[df['index'].str.contains('frag_')][['formula', 'glf_intensity']].values]).astype(str)

	if args.loss_types == 'parent':
		parent_index = df[df['index'].str.contains('frag_')].sort_values(by = 'mz', ascending = False).iloc[0]['index']
		df = df[(df['index'].str.contains('loss_')) & (df['from_index'] == parent_index)]

	if np.any(df['index'].str.contains('loss_')) and (not args.loss_types != 'none'):
		doc_losses = np.concatenate([[f"loss_{w}" for _ in range(round(c))] for w, c in df[df['index'].str.contains('loss_')][['formula', 'glf_intensity']].values]).astype(str)
	else:
		doc_losses = np.array([])
	doc = np.concatenate([doc_frags, doc_losses])
	return(doc)


for n_train_spectra, _ in enumerate(mgf.MGF(args.train_mgf)):
	continue
n_train_spectra += 1
f = mgf.MGF(args.train_mgf)
for spec in tqdm(f, total = n_train_spectra):
	try:
		doc = make_doc(spec, args.documents_dir)
		s = spec['params']['smiles']
		sid = spec['params']['id']		
		s_fname = os.path.join(args.documents_dir, sid + '.tsv')
		labs = get_mes_labels(s)
		di = mdl.add_doc(doc, labels = labs)
		doc_indices.append(di)
		doc_smiles.append(s)
		doc_fnames.append(s_fname)
		doc_sids.append(sid + '_train')

	except:
		s_fname = os.path.join(args.documents_dir, sid + '.tsv')
		print(s_fname)
		sys.stdout.flush()
		continue
os.system(f'mkdir -p {args.out_dir}')
df_out = pd.DataFrame({'doc_index' : doc_indices, 'fname' : doc_fnames, 'smiles' : doc_smiles})
df_out.to_csv(os.path.join(args.out_dir, 'df_train_indices.tsv'), sep = '\t')

########################################################################

for iteration in tqdm(range(0, int(args.num_iterations) + 10, 10)):
	if iteration % 500 == 0:
		out_fname = os.path.join(args.out_dir, f"iter_{iteration}.tpy")
		mdl.save(out_fname)
	mdl.train(10)

########################################################################
# Calculate predictions at last iteration
########################################################################
df_labels = pd.read_csv(args.df_labels, sep = '\t', index_col = 0)
topic_words = np.array([mdl.get_topic_word_dist(i) for i in range(len(mdl.topic_label_dict))])
dftw = pd.DataFrame(topic_words, index = list(mdl.topic_label_dict), columns = mdl.vocabs).T
# dftw.to_csv(os.path.join(args.out_dir, "dftw.tsv"), sep = '\t')
df_sims = pd.DataFrame(index = dftw.columns)
t3s = []
sids_test = []
for n_test_spectra, _ in enumerate(mgf.MGF(args.test_mgf)):
	continue
n_test_spectra += 1
for spec in tqdm(mgf.MGF(args.test_mgf), total = n_test_spectra):
	doc = make_doc(spec, args.documents_dir)
	s = spec['params']['smiles']
	sid = spec['params']['id']	
	s_fname = os.path.join(args.documents_dir, sid + '.tsv')	
	# labs = get_mes_labels(s)
	nc = np.unique(doc, return_counts=True)
	dfc = pd.DataFrame(nc[1], columns = ['current_doc'], index = nc[0]).reindex(dftw.index, fill_value = 0.0)
	if dfc['current_doc'].sum() == 0:
		print(f'Error no overlap {s_fname}')
		cos_sims = [0 for x in dftw.columns]
	else:
		cos_distances = [scipy.spatial.distance.cosine(dfc['current_doc'], dftw[c]) for c in dftw.columns]
		cos_sims = [1 - x for x in cos_distances]
	df_temp = pd.DataFrame({sid : cos_sims}, index = dftw.columns)
	df_sims = pd.concat([df_sims, df_temp], axis = 1)
	sids_test.append(sid + '_test')
	t3 = df_labels.loc[sid + '_test'][df_sims[sid].sort_values(ascending = False).iloc[0 : 3].index].sum()
	t3s.append(t3)
df_t3s = pd.DataFrame({'sid' : sids_test, 't3' : t3s})
df_sims.to_csv(os.path.join(args.out_dir, "df_preds.tsv"), sep = '\t')

with open(os.path.join(args.out_dir, 'log.txt'), 'a') as f:
	print(f'\nt3>=1 - {len([x for x in t3s if x >= 1])} / {n_test_spectra}', file = f)
	print(f'\nt3>=2 - {len([x for x in t3s if x >= 2])} / {n_test_spectra}', file = f)

