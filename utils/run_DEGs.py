REPO = '..'
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
obs_path = f"{RESULT_TABLE}/GEX_OBS_Lineage.csv"
obs = pd.read_csv(obs_path,index_col=0)
adata = sc.read(f'{RESULT_OBJ}/gex_all.h5ad')
# keep cells with annotation
adata = adata[obs.index,:]
# add clinical information
sample_meta =  pd.read_excel(f'{RESULT_TABLE}/Supplementary Table 1.xlsx',index_col=0).replace(np.nan,'N/A')

# append meta information to the adata.obs data frame
adata.obs = adata.obs.reset_index().merge(sample_meta,left_on='Sample',right_on='CCG_ID',how='left').set_index('index')

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
pdata.write(f"{RESULT_OBJ}/pseudobulk.h5ad")

######################## DEGs per cellstates #####################
results = oneToRestDEGs(pdata=pdata,
                        group='Cellstate',
                        confounders=['Patient'],
                        output_xlsx=f"{Manuscript_RESULT}/DEGs_Cellstate.xlsx")
## CIBERSORTx signature matrix
signature_matrix = get_signature_matrix(pdata=pdata,
                                        degs_results=results,
                                        group='Cellstate',
                                        exclude_groups=['Immune','Stromal'],
                                        )
signature_matrix.to_csv(f"{Manuscript_RESULT}/cellstate_signature.tsv",sep='\t')
################################################################################################
######################## DEGs per celltype #####################
results = oneToRestDEGs(pdata=pdata,
                        group='Celltype',
                        confounders=['Patient'],
                        output_xlsx=f"{Manuscript_RESULT}/DEGs_Celltype.xlsx")
## CIBERSORTx signature matrix
signature_matrix = get_signature_matrix(pdata=pdata,
                                        degs_results=results,
                                        group='Cellstate',
                                        exclude_groups=['Immune','Stromal'],
                                        )
signature_matrix.to_csv(f"{Manuscript_RESULT}/celltype_signature.tsv",sep='\t')
################################################################################################
######################## DEGs per cellstates per celltype #####################
for celltype in ['CD8T','Macs','Tumor','Endothelial','CAF']:
    adata_sub = pdata[pdata.obs['Celltype'] == celltype].copy()
    results = oneToRestDEGs(pdata=adata_sub,
                        group='Cellstate',
                        confounders=['Patient'],
                        output_xlsx=f"{Manuscript_RESULT}/DEGs_{celltype}.xlsx")