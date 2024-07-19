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


pdata =sc.read(f"{Manuscript_RESULT}/pseudobulk.h5ad")


adata_sub = pdata[pdata.obs['Celltype'] == 'Tumor'].copy()
adata_sub.obs['aggr_cellstate'] = adata_sub.obs.Cellstate.map({
    'Tumor.EMT-III':'EMT',
    'Tumor.EMT-II':'EMT',
    'Tumor.ER-II':'ER',
    'Tumor.ER-I':'ER',
    'Tumor.Cell_Cycle':'Cell_Cycle',
    'Tumor.Interferon/MHCII(I)':'Interferon'
})

results = oneToRestDEGs(pdata=adata_sub,
                    group='aggr_cellstate',
                    confounders=['Patient'],
                    output_xlsx=f"{Manuscript_RESULT}/DEGs_Tumor_Aggr.xlsx")