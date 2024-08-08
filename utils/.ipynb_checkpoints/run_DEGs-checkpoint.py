REPO = '.'
RESULT_TABLE = f'{REPO}/data/table'
RESULT_OBJ = f'{REPO}/data/object'

import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
import sys
sys.path.append(REPO)
from utils.misc import oneToRestDEGs,get_signature_matrix

# Define Paramters
#-------- Cell annotation
obs = pd.read_csv(f'{RESULT_TABLE}/GEX_OBS_Lineage.csv',index_col=0)
## add celltype 
obs['Celltype'] = obs['Lineage']
for lineage in ['T','Myeloid','Epithelial']:
    lineage_anno = pd.read_csv(f'{RESULT_TABLE}/annotation/{lineage}.csv',index_col=0)
    celltype = lineage_anno['Celltype'] if 'Celltype' in lineage_anno else lineage_anno['Lineage']
    obs.loc[celltype.index,'Celltype'] = celltype
## add cellstate
obs['Cellstate'] = obs['Celltype']
for celltype in ['CD8T','Macs','Tumor']:
    cellstate = pd.read_csv(f'{RESULT_TABLE}/MPs/{celltype}/Annotation.csv',index_col=0)['Cellstate']
    obs.loc[cellstate.index,'Cellstate'] = cellstate

# clin
sample_meta =  pd.read_excel(f'{RESULT_TABLE}/Supplementary Table 1.xlsx',index_col=0).replace(np.nan,'N/A')
sample_meta['br_short'] = sample_meta['BestResponse'].map({
    'favorable response\n(RCB 0-I)': 'R',
    'unfavorable response\n(RCB II-III)': 'NR'
})
# merge sample meta with cell meta
obs = obs.reset_index().merge(sample_meta,left_on='Sample',right_on='CCG_ID',how='left').set_index('index')
#-------- load adata
adata = sc.read(f'{RESULT_OBJ}/gex_all.h5ad')
adata = adata[obs.index,:]
adata.obs = obs

## Generate pseudo-bulk
adata.layers['counts']=adata.X
pdata = dc.get_pseudobulk(
    adata,
    sample_col='CCG_ID',
    groups_col='Cellstate',
    layer='counts',
    mode='sum',
    min_cells=10, # remove samples with some cell states less than 10 cells -> from decoupler document
    min_counts=1000 # remove samples with some cell states less than 1000 accumulated counts
)


######################## DEGs per cellstates #####################
results = oneToRestDEGs(pdata=pdata,
                        group='Cellstate',
                        confounders=['Patient'],
                        output_xlsx=f"{RESULT_TABLE}/DEGs_Cellstate.xlsx")
## CIBERSORTx signature matrix
signature_matrix = get_signature_matrix(pdata=pdata,
                                        degs_results=results,
                                        group='Cellstate',
                                        exclude_groups=['Immune','Stromal'],
                                        )
signature_matrix.to_csv(f"{RESULT_TABLE}/cellstate_signature.tsv",sep='\t')
################################################################################################
######################## DEGs per celltype #####################
results = oneToRestDEGs(pdata=pdata,
                        group='Celltype',
                        confounders=['Patient'],
                        output_xlsx=f"{RESULT_TABLE}/DEGs_Celltype.xlsx")
## CIBERSORTx signature matrix
signature_matrix = get_signature_matrix(pdata=pdata,
                                        degs_results=results,
                                        group='Cellstate',
                                        exclude_groups=['Immune','Stromal'],
                                        )
signature_matrix.to_csv(f"{RESULT_TABLE}/celltype_signature.tsv",sep='\t')
################################################################################################
######################## DEGs per cellstates per celltype #####################
for celltype in ['Tumor']:
    adata_sub = pdata[pdata.obs['Celltype'] == celltype].copy()
    results = oneToRestDEGs(pdata=adata_sub,
                        group='Cellstate',
                        confounders=['Patient'],
                        output_xlsx=f"{RESULT_TABLE}/DEGs_{celltype}.xlsx")