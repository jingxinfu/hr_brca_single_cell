GEX_Cohort='GEX_CCG1112_LowMt'
REPO = './'
WORKFLOW_DATA = f'{REPO}/data/workflow'
Manuscript_RESULT = f'{REPO}/data/result/manuscript_table/'

import scanpy as sc
import sys
sys.path.append(REPO)
from utils.misc import oneToRestDEGs

# Define Paramters
pdata= sc.read(f"{Manuscript_RESULT}/pseudobulk.h5ad")
pdata = pdata[pdata.obs.Timepoint=='Baseline',:]
######################## DEGs (NR vs R) per celltype #####################
for celltype in ['CD8T','Macs','Tumor','Endothelial','CAF']:
    adata_sub = pdata[pdata.obs['Celltype'] == celltype].copy()
    results = oneToRestDEGs(pdata=adata_sub,
                        group='BestResponse',
                        output_xlsx=f"{Manuscript_RESULT}/DEGs_Response/{celltype}.xlsx")