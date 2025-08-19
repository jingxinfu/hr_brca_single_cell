REPO = '../'
RESULT_TABLE = f'{REPO}/data/table'
RESULT_OBJ = f'{REPO}/data/object'
FIGURE_FOLDER= f'{REPO}/data/figure'
SETTING_FOLDER = f'{REPO}/data/setting'
EXTERNAL_DATA=f'{REPO}/data/external'
# load Terra API for get data table from Terra
import sys
sys.path.append(REPO)
from settings import COLOR_PAlETTE
from utils.visual import *

import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scPipe as sp
import scanpy as sc
import seaborn as sns
import warnings
import commentjson
### Additional Colormap
with open(f'{SETTING_FOLDER}/ColorMap.json','r') as f:
    colormap = commentjson.load(f)
COLOR_PAlETTE.update(colormap)

for POPULATION in [f'TumorRandom95_{i}' for i in range(10)]:

    obs_path = f'{RESULT_TABLE}/TumorRandom95_10times.csv'
    adata = sc.read(f'{RESULT_OBJ}/gex_all.h5ad')
    obs = pd.read_csv(obs_path,index_col=0)
    adata = adata[adata.obs.index.isin(obs.index[obs[POPULATION]=='selected']),:]
    print(f"Detect N={adata.shape[0]:,} tumor cells with high quality GEX profiles;")
    print(f"max(% of mitocondrial reads): {adata.obs['pct_counts_mito'].max():.0f}")
    sample_meta =  pd.read_excel(f'{RESULT_TABLE}/Supplementary Table 1.xlsx',index_col=0).replace(np.nan,'N/A')
    import pickle
    with open(f'{RESULT_TABLE}/MPs/{POPULATION}/nmf_basis.pickle', "rb") as input_file:
        programs_basis = pickle.load(input_file)
    
    with open(f'{RESULT_TABLE}/MPs/{POPULATION }/nmf_coef.pickle', "rb") as input_file:
        programs_coef = pickle.load(input_file)
        
    n_programs = programs_basis[list(programs_basis.keys())[2]].shape[1] * len(programs_basis)
    
    with open(f'{RESULT_TABLE}/MPs/{POPULATION}/robustRPH.pickle', "rb") as input_file:
        robustRPH = pickle.load(input_file)
    
    import palettable
    Cluster_Map,MP_Genesets,Programs_Order,MP_colors = sp.ext.clusterRobustRPH(
        robustRPH=robustRPH,
        programs_basis=programs_basis,
        Min_group_size=2,
        palette=palettable.tableau.Tableau_20.hex_colors + palettable.tableau.TrafficLight_9.hex_colors
    )
    
    low_quality_mps = pd.concat([
        MP_Genesets.apply(lambda c:c.str.startswith('MT-').sum(),axis=0).rename('Mitochondrial'),
        MP_Genesets.apply(lambda c:c.str.contains('^RP[LS]').sum(),axis=0).rename('Ribosomal')],axis=1
                               ).sort_values(['Mitochondrial','Ribosomal'],ascending=False).sum(axis=1).idxmax()
    
    anno = Cluster_Map.rename('MetaProgram').to_frame()
    anno['Sample'] = anno.index.map(lambda x: x.split('.')[0])
    anno = anno.reset_index().merge(sample_meta,left_on='Sample',right_on='CCG_ID',how='left').set_index('index')

    anno = anno.loc[~anno['MetaProgram'].isin([low_quality_mps]+['MP_Unknown']),:]
    New_Programs_Order = [ x for x in Programs_Order if x in anno.index]
    New_MP_Genesets = MP_Genesets[anno['MetaProgram'].unique().tolist()]
    print(f"Generated {anno['MetaProgram'].nunique()} MP clusters which covered {(anno.shape[0]/robustRPH.shape[1]):.0%} robust RPH programs.")

    mp_gmt_path = f'{RESULT_TABLE}/MPs/{POPULATION}/MP_Programs.gmt'
    with open(mp_gmt_path,'w') as f:
        for mp in New_MP_Genesets:
            f.write('\t'.join([mp,'MetaProgram_HRpos']+MP_Genesets[mp].tolist())+'\n')
    print(f'FINISH {POPULATION}')

