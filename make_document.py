########################################################################
# Required packages
########################################################################
import argparse
import sys
import os
import pandas as pd
import numpy as np
from tqdm.auto import tqdm
from pyteomics import mgf, auxiliary
import molmass
import subprocess
from joblib import Parallel, delayed
########################################################################
# Parse arguments
########################################################################
parser = argparse.ArgumentParser() 
def file_choices(choices,fname, argname):
	ext = os.path.splitext(fname)[1]
	if ext not in choices:
	   parser.error(f"{argname} must be one of the following filetypes ({choices})")
	return fname
parser.add_argument('--in_mgf', required = True)
parser.add_argument('--out_dir', required = True)
parser.add_argument('--eval_peak_script', required = True)
parser.add_argument('--n_jobs', default = 1, type = int)
parser.add_argument('--adduct_element', default = 'H', type = str)
parser.add_argument('--adduct_number', default = 1, type = int)
args = parser.parse_args()
########################################################################
# Auxiliary formula parsing functions
########################################################################
def add_fd(form1, fd2):
	fd1 = {x[0] : x[1] for x in molmass.Formula(form1).composition()}
	for add_element, add_num in fd2.items():
		if add_element not in fd1.keys():
			fd1[add_element] = 0
		fd1[add_element] += add_num
	s_out = ''.join([f"{k}{v}" if v > 1 else f"{k}" for k, v in fd1.items()])
	return(molmass.Formula(s_out).formula)

def add_formulas(form1, form2):
	fd1 = {x[0] : x[1] for x in molmass.Formula(form1).composition()}
	fd2 = {x[0] : x[1] for x in molmass.Formula(form2).composition()}
	for add_element, add_num in fd2.items():
		if add_element not in fd1.keys():
			fd1[add_element] = 0
		fd1[add_element] += add_num
	s_out = ''.join([f"{k}{v}" if v > 1 else f"{k}" for k, v in fd1.items()])
	return(molmass.Formula(s_out).formula)

def sub_fd(form1, fd2, mode = 'report'): # form1 - fd2
	fd1 = {x[0] : x[1] for x in molmass.Formula(form1).composition()}
	for sub_element, sub_num in fd2.items():
		if sub_element not in fd1.keys():
			fd1[sub_element] = 0
		fd1[sub_element] -= sub_num
	if mode == 'report':
		s_out = ""
		for k, v in fd1.items():
			if v != 0:
				s_out += k
				if v != 1:
					s_out += str(v)
		return(s_out)
	elif mode == 'error':
		for k, v in fd1.items():
			if v < 0:
				raise Exception(f"Error - {v} {k} in resulting formula")
	elif mode == 'ignore':
		s_out = ""
		for k, v in fd1.items():
			if v > 0:
				s_out += k
				if v > 1:
					s_out += str(v)
	return(molmass.Formula(s_out).formula)

def sub_formulas(form1, form2, mode = 'report'): # form1 - form2
	fd1 = {x[0] : x[1] for x in molmass.Formula(form1).composition()}
	fd2 = {x[0] : x[1] for x in molmass.Formula(form2).composition()}
	for sub_element, sub_num in fd2.items():
		if sub_element not in fd1.keys():
			fd1[sub_element] = 0
		fd1[sub_element] -= sub_num
	if mode == 'report':
		s_out = ""
		for k, v in fd1.items():
			if v != 0:
				s_out += k
				if v != 1:
					s_out += str(v)
		return(s_out)
	elif mode == 'error':
		for k, v in fd1.items():
			if v < 0:
				raise Exception(f"Error - {v} {k} in resulting formula")
	elif mode == 'ignore':
		s_out = ""
		for k, v in fd1.items():
			if v > 0:
				s_out += k
				if v > 1:
					s_out += str(v)
	return(molmass.Formula(s_out).formula)
########################################################################
# Create a document (in tsv format) for each spectrum in the mgf file
########################################################################
for num_spectra, s in enumerate(mgf.MGF(args.in_mgf)):
	continue
num_spectra += 1
for s in tqdm(mgf.MGF(args.in_mgf), total = num_spectra):
	df = pd.DataFrame({'mz' : s['m/z array'], 'intensity' : s['intensity array']})
	if args.adduct_number >= 0:
		form = add_fd(s['params']['formula'], {args.adduct_element : args.adduct_number})
	else:
		form = sub_fd(s['params']['formula'], {args.adduct_element : args.adduct_number})
	out_fname = os.path.join(args.out_dir, s['params']['id'] + '.tsv')
	def get_output(i_mz, mz):        
		cmd = f"Rscript {args.eval_peak_script} --formula {form} --mz_val {mz} --tolerance {0.1}"
		ret = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = ret.stdout.read().decode()
		stderr = ret.stderr.read().decode()
		if '''Error in .jcall(mfSet, "I", "size")''' in stderr:
			a = None
			tm = None
			o = None
		elif "cdkFormula:" in stdout:
			a = molmass.Formula(stdout.split(' ,')[0].split(' ')[-1].strip()).formula
			# a = molmass.Formula(stdout.split(' ,')[0].strip('cdkFormula: ')).formula
			tm = float(stdout.split(' , ')[1].strip('mass = '))
			o = stdout
		else:
			sys.exit(f'Error in rcdk run\n stdout {stdout}\nstderr {stderr}')
		return(a, tm, o)
	output_list = Parallel(n_jobs = args.n_jobs)(delayed(get_output)(i_mz, m) for i_mz, m in enumerate(df['mz']))
	rcdk_assignments = []
	rcdk_theo_masses = []
	rcdk_raw_output = []
	for a, tm, o in output_list:
		rcdk_assignments.append(a)
		rcdk_theo_masses.append(tm)
		rcdk_raw_output.append(o)
	df['rcdk_assignment_top'] = rcdk_assignments
	df['rcdk_top_mass'] = rcdk_theo_masses
	df['rcdk_all_outputs'] = rcdk_raw_output
	
	df = df.dropna(subset = ['rcdk_assignment_top'])
	df['formula'] = df['rcdk_assignment_top']
	df = df.reset_index(drop=True)
	df['index'] = [f'frag_{i}' for i in df.index]

	pis = np.where(df['mz'].values - df['mz'].values[:, None] > 0)

	a = df.values[np.array(pis[::-1]).T]
	from_indices = np.array(pis[::-1]).T[: ,0]
	to_indices = np.array(pis[::-1]).T[: ,1]
	f = lambda x : sub_formulas(*x)
	loss_forms = np.array(list(map(f, a[:, :, 2])))
	dmzs = np.subtract(a[:, :, 0][:, 0], a[:, :, 0][:, 1])
	from_intensities = a[:,:,1][:, 0]
	to_intensities = a[:,:,1][:, 1]
	mean_intensities = np.round(np.mean(a[:,:,1], axis = 1).astype(float), decimals = 5)

	dl = pd.DataFrame({'mz' : dmzs, 
				  'formula' : loss_forms, 
				  'intensity' : mean_intensities,
					   'from_index' : [f'frag_{x}' for x in from_indices], 
					   'to_index' : [f'frag_{x}' for x in to_indices]})

	dl = dl[(from_intensities - to_intensities) > 0]
	if len(dl) > 0:
		dl = dl[~dl['formula'].str.contains('-')]
		dl = dl.reset_index(drop=True)
		dl['index'] = [f'loss_{i}' for i in dl.index]
		dl = dl[dl['formula'] != '']


	dl = dl.reset_index(drop=True)
	dl['index'] = [f'loss_{i}' for i in dl.index]
	df_out = pd.concat([df, dl], sort = False)
	
	
	df_out.to_csv(out_fname, sep = '\t', index = False)
