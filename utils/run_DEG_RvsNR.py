REPO = '.'
RESULT_TABLE = f'{REPO}/data/table'
RESULT_OBJ = f'{REPO}/data/object'

import scanpy as sc
import sys
sys.path.append(REPO)
from utils.misc import oneToRestDEGs
import numpy as np
import pandas as pd
import decoupler as dc

# Define Paramters
obs_path = f"{RESULT_TABLE}/GEX_OBS_Lineage.csv"
obs = pd.read_csv(obs_path,index_col=0)
adata = sc.read(f'{RESULT_OBJ}/gex_all.h5ad')
# keep cells with annotation
adata = adata[obs.index,:]
# add clinical information
sample_meta =  pd.read_excel(f'{RESULT_TABLE}/Supplementary Table 1.xlsx',index_col=0).replace(np.nan,'N/A')
sample_meta['br_short'] = sample_meta['BestResponse'].map({
    'favorable response\n(RCB 0-I)': 'R',
    'unfavorable response\n(RCB II-III)': 'NR'
})
# append meta information to the adata.obs data frame
adata.obs = adata.obs.reset_index().merge(sample_meta,left_on='Sample',right_on='CCG_ID',how='left').set_index('index')

adata.obs['Celltype'] = 'Misc.'
## add celltype 
for lineage in ['T','Myeloid','Epithelial']:
    lineage_anno = pd.read_csv(f'{RESULT_TABLE}/annotation/{lineage}.csv',index_col=0)
    celltype = lineage_anno['Celltype'] if 'Celltype' in lineage_anno else lineage_anno['Lineage']
    adata.obs.loc[celltype.index,'Celltype'] = celltype

## Generate pseudo-bulk
adata.layers['counts']=adata.X
pdata = dc.get_pseudobulk(
    adata,
    sample_col='CCG_ID',
    groups_col='Celltype',
    layer='counts',
    mode='sum',
    min_cells=10, # remove samples with some cell states less than 10 cells -> from decoupler document
    min_counts=1000 # remove samples with some cell states less than 1000 accumulated counts
)

pdata = pdata[pdata.obs.Timepoint=='Baseline',:]
######################## DEGs (NR vs R) per celltype #####################
for celltype in ['CD8T','Macs','Tumor']:
    adata_sub = pdata[pdata.obs['Celltype'] == celltype].copy()
    results = oneToRestDEGs(pdata=adata_sub,
                        group='br_short',
                        output_xlsx=f"{RESULT_TABLE}/DEGs_Response/{celltype}.xlsx")