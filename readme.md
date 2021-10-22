# README

gkreder@gmail.com

Please note that this github repo only contains the **code** directory mentioned below. The **data** directory (and a copy of the same code) is available in the Zenodo repository at https://doi.org/10.5281/zenodo.4655149

## Data

The **data** is available in the Zenodo repository https://doi.org/10.5281/zenodo.4655149

The **data** directory contains the data necessary for running the experiments described in the paper. Specifically the files included are the following:

| File / Directory      | Description |
| ----------- | ----------- |
| **df_substructs.tsv**   | The substructures being predicted in the default models (includes substructure names and SMARTS strings |
| df_labels.tsv      | The substructure labels for each spectrum in the default train/validate/test set of spectra |
| df_ids.tsv      | The spectrum ids and SMILES strings for each spectrum in the default train/validate/test set of spectra |
| **train.mgf**   | The default train spectra in .mgf format |
| **test.mgf**   | The default test spectra in .mgf format |
| **validate.mgf**   | The default validation spectra in .mgf format |
| binarized\_embedded\_spectra.npz   | The vector-embedded binarized spectra in the default train/test/validate sets |
| documents   | The geneated documents with fitted formulas and neutral losses (in tsv format) for each spectrum in the default train/validate/test set of spectra|

The data files listed in bold specify a given train/test/validate + substructure set and are used to generate those files that are not listed in bold. We include the generated files as a convenience but note that the train/test/validate set and substructures can be changed by the user as desired. The included bolded files correspond to the default train/test/validate spectra and substructures as provided in the MESSAR paper [1]. Details on how to generate the data files are shown below.

[1] Liu Y, Mrzic A, Meysman P, De Vijlder T, Romijn EP, et al. (2020) MESSAR: Automated recommendation of metabolite substructures from tandem mass spectra. PLOS ONE 15(1): e0226770. https://doi.org/10.1371/journal.pone.0226770

## Code

The **code** directory contains the scripts neccesary for preparing data, generating documents, running the LLDA and k-NN models, and calculating ROC performances from both models.

### Dependencies

The recommended software environment is provided in a Docker image containing all necessary dependencies available at [https://hub.docker.com/r/gkreder/ms2-topic-model](https://hub.docker.com/r/gkreder/ms2-topic-model). A Dockerfile used to build this image is included is the code package. If the user would like to install dependencies manually, they are the following:

#### R
- stringr 
- docopt
- rcdk

#### Python
- tqdm
- tomotopy
- rdkit
- sckit-learn
- scipy
- matplotlib
- seaborn
- pyteomics
- pandas
- xlrd
- molmass


All of which should be installable via Pip or Conda. A requirements.txt file is included in the code package that can be installed via Conda. This will install only the Python requirements, the R requiments must be installed manually.  

### Generating data and documents


To generate the df\_labels, df\_ids, and embedded\_spectra files, the user may run the following command:

```bash
python prep_data.py \
--train_mgf train.mgf \
--test_mgf test.mgf \
--validate_mgf validate.mgf \
--df_substructs df_substructs.tsv \
--embedded_spectra_out <embedded_spectra_fileName_out.npz> \
--df_labels_out <df_labels_fileName_out.tsv> \
--df_ids_out <df_ids_fileName_out.tsv>
```

Note that the extensions in the output files must be .npz and .tsv respectively. Also note that the validation set may be left empty by passing through an empty .mgf file.

To generate documents from a given mgf file, the user may run the following command:

```bash
python make_documents.py \
--in_mgf <in_mgf_file> \
--out_dir <documents_dir> \
--eval_peak_script evaluate_peak.R \
--n_jobs 1 \
--adduct_element H \
--adduct_element 1
```

This script must be pointed to the evaluate_peak.R script included with this code base. Note that this is a fairly time consuming step. In our tests generating a single spectrum using a single thread could take more than 2 minutes depending on the spectrum. The --n\_jobs flag allows for parallelization of formula finding over multiple threads and speeds up the process considerably if multiple cores are available. The --adduct\_element and --adduct\_number arguments specify if an additional element should be added to each spectrum's molecular formula before fitting molecular formulas. 

### Running LLDA and k-NN models

Once the data has been generated (or if the user is using the pre-generated data included in the code package), the LLDA model can be run using the following command:

```bash
python run_llda.py \
--Q <Q> \
--B <B> \
--out_dir <llda_out_directory> \
--train_mgf train.mgf \
--test_mgf test.mgf \
--documents_dir documents \
--df_substructs df_substructs.tsv \
--df_labels df_labels.tsv \
--num_iterations <number_iterations>
```

This produces the following files:

| File      | Description |
| ----------- | ----------- |
| log.txt      | Contains the LLDA command and parameters passed through to the model as well as the t3 >= 1 and t3 >= 2 metrics for this run |
| df_preds.tsv      | The predicted scores for each substructure in each test spectrum |
| iter\_*x*.tpy      | The LLDA model at iteration *x* (which can be loaded using Tomotopy) |
| df\_train_indices.tsv      | The document index in the Tomotopy LLDA model for each document in the train spectra along with its corresponding filename and SMILES string |

The k-NN model can be run using the following command:

```bash
python run_baseline.py \
--Q <Q> \
--B <B> \
--out_dir <llda_out_directory> \
--train_mgf train.mgf \
--test_mgf test.mgf \
--documents_dir documents \
--df_substructs df_substructs.tsv \
--df_labels df_labels.tsv \
--num_iterations <number_iterations>
```

| File      | Description |
| ----------- | ----------- |
| log.txt      | Contains the k-NN command and parameters passed through to the model as well as the t3 >= 1 and t3 >= 2 metrics for this run |
| df_preds.tsv      | The predicted scores for each substructure in each test spectrum |

### Calculating ROC performance

After running both the LLDA and k-NN models - to calculate the per-substructure area under the receiver operating characteristic curve (AUC), the user can run the following command:

```bash
python roc_calcs.py \
--df_labels df_labels.tsv \
--df_preds_knn <knn_df_preds_tsv> \
--df_preds_llda <llda_df_preds_tsv> \
--df_substructs df_substructs.tsv \
--out_dir <roc_output_directory>
```

This produce the following files in the output directory

| File      | Description |
| ----------- | ----------- |
| roc_aucs.pdf      | A plot of the per-substructure ROC AUC in both the k-NN and LLDA models, colored by train set appearance |
| precs.pdf      | A plot of the per-substructure average precision score in both the k-NN and LLDA models, colored by train set appearance |
| df_rscores.tsv      | A per-substructure breakdown of ROC AUC, average precision score, train appearance, test appearance, and number of atoms |


