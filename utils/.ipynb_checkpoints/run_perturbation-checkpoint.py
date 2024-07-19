import os
import glob
import dill
import pyranges
import warnings
warnings.filterwarnings("ignore")
import sys
import scanpy as sc
import scPipe as sp
import pickle
import pandas as pd

# scenicplus
import numpy as np
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus
import pybiomart as pbm

ATAC_Cohort='CCG1112_ATAC_QC'
GEX_Cohort='GEX_CCG1112_LowMt'
REPO = '/home/analysis/hr_brca_mo_analysis'
WORKFLOW_DATA = f'{REPO}/data/workflow'
EXTERNAL_DATA = f'{REPO}/data/external'
RESULT_DATA = f'{REPO}/data/result/cleaned_files'
RESULT_Figure = f'{REPO}/data/result/Figure'
Group_Variable = 'GEX_Celltype'


_stderr = sys.stderr
null = open(os.devnull,'wb')
import dill
scplus_obj = dill.load(open(f'{RESULT_DATA}/scenicplus/scplus_obj.pkl', 'rb'))

from scenicplus.dimensionality_reduction import run_eRegulons_pca
run_eRegulons_pca(
        scplus_obj,
        auc_key = 'eRegulon_AUC_filtered',
        reduction_name = 'eRegulons_PCA_gene_based',
        n_pcs=25,
        selected_regulons = scplus_obj.uns['selected_eRegulon']['Gene_based'])
from pycisTopic.diff_features import find_highly_variable_features
hvg = find_highly_variable_features(scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features = 200, plot = False)
from scenicplus import diff_features
diff_features.get_differential_features(scplus_obj,variable=Group_Variable)

flatten_list = lambda t: [item for sublist in t for item in sublist]
DEGs = list(set(flatten_list([list(scplus_obj.uns['DEGs'][Group_Variable][k].index) for k in scplus_obj.uns['DEGs'][Group_Variable].keys()])))
genes_to_use = list(set([*DEGs, * [x.split('_')[0] for x in scplus_obj.uns['selected_eRegulon']['Gene_based']]]))

from scenicplus.simulation import train_gene_expression_models
regressors = train_gene_expression_models(
        scplus_obj,
        eRegulon_metadata_key = 'eRegulon_metadata',
        genes = genes_to_use)

from scenicplus.simulation import simulate_perturbation
from scenicplus.simulation import _make_rankings
TFs_of_interest = list(set([x.split('_')[0] for x in scplus_obj.uns['selected_eRegulon']['Gene_based']]))
n_iter = 3

from tqdm.notebook import tqdm
import logging
logging.basicConfig(level=logging.CRITICAL)
import warnings
from scenicplus.eregulon_enrichment import score_eRegulons
from scenicplus.dimensionality_reduction import run_eRegulons_pca
for i,TF in enumerate(TFs_of_interest):
    print(f'##### Processing {i}/{len(TFs_of_interest)}#####')
    perturbed_matrices = simulate_perturbation(
        scplus_obj,
        perturbation = {TF: 0},
        regressors = regressors,
        keep_intermediate = True,
        n_iter = n_iter)
    perturbed_rankings = {k: _make_rankings(perturbed_matrices[k]) for k in perturbed_matrices.keys()}
    for k in perturbed_rankings.keys():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            score_eRegulons(
                scplus_obj,
                ranking = perturbed_rankings[k],
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = f'{TF}_KD_sim_eRegulon_AUC_iter_{k}',
                enrichment_type = 'gene',
                n_cpu = 16
            )
from scenicplus.simulation import _project_perturbation_in_embedding
shifts_PC0 = {}
shifts_PC1 = {}
import sys
sys.stderr = open(os.devnull, "w")  # silence stderr
for TF in TFs_of_interest:
    delta_embedding = _project_perturbation_in_embedding(
        scplus_obj,
        original_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_0']['Gene_based'],
        perturbed_matrix = scplus_obj.uns[f'{TF}_KD_sim_eRegulon_AUC_iter_4']['Gene_based'],
        reduction_name = f'{TF}_KD_sim_eRegulons_PCA_iter_0')
    mean_shift = pd.DataFrame(delta_embedding).groupby(scplus_obj.metadata_cell['GEX_Celltype'].to_numpy()).mean()
    shifts_PC0[TF] = mean_shift[0]
    shifts_PC1[TF] = mean_shift[1]
sys.stderr = sys.__stderr__  # unsilence stderr
shift_df = pd.DataFrame(shifts_PC0).T
factors_to_plot = [
        *shift_df.max(1).sort_values(ascending = False).head(10).index,
        *reversed(shift_df.min(1).sort_values(ascending = True).head(10).index)]
# line_order = ['MM001', 'MM011', 'MM031', 'MM087', 'MM074', 'MM057', 'MM047', 'MM029', 'MM099']
line_order = scplus_obj.metadata_cell[Group_Variable].unique().tolist()
import seaborn as sns
fig, ax = plt.subplots(figsize = (10, 5))
sns.heatmap(
    shift_df.loc[factors_to_plot, line_order].T,
    yticklabels=True,vmin = -0.3, vmax = 0.3, ax = ax, cmap = 'bwr')
for ytick in ax.get_yticklabels():
    ytick.set_color(color_dict_line[ytick.get_text()])
        
fig.savefig(f'{RESULT_Figure}/{Group_Variable}_TF_Perturbation.pdf',bbox_inches='tight')
shift_df.to_csv(f'{RESULT_DATA}/{Group_Variable}_TF_Perturbation.csv')