GEX_Cohort='GEX_CCG1112_LowMt'
REPO = './'
WORKFLOW_DATA = f'{REPO}/data/workflow'
Manuscript_RESULT = f'{REPO}/data/result/manuscript_table/'


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
obs_path = f"{Manuscript_RESULT}/GEX_OBS.csv"
obs = pd.read_csv(obs_path,index_col=0)
adata = sc.read(f'{WORKFLOW_DATA}/{GEX_Cohort}/gex_qc.h5ad')
# keep cells with annotation
adata = adata[obs.index,:]
# append meta information to the adata.obs data frame
for c in ['Cellstate','Celltype','BestResponse','Patient','Timepoint','Sample_Short','Treatment_Arm','RCB']:
    adata.obs[c] = obs.loc[adata.obs.index,c]

## Generate pseudo-bulk
adata.layers['counts']=adata.X
pdata = dc.get_pseudobulk(
    adata,
    sample_col='Sample_Short',
    groups_col='Cellstate',
    layer='counts',
    mode='sum',
    min_cells=10, # remove samples with some cell states less than 10 cells -> from decoupler document
    min_counts=1000 # remove samples with some cell states less than 1000 accumulated counts
)
pdata.write(f"{Manuscript_RESULT}/pseudobulk.h5ad")

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