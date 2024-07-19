GEX_Cohort='GEX_CCG1112_LowMt'
REPO = './'
WORKFLOW_DATA = f'{REPO}/data/workflow'
Manuscript_RESULT = f'{REPO}/data/result/manuscript_table/'

import pandas as pd
import scanpy as sc
import sys
sys.path.append(REPO)


# Define Paramters
obs_path = f"{Manuscript_RESULT}/GEX_OBS.csv"
degs_path = f"{Manuscript_RESULT}/DEGs_Cellstate.xlsx"
degs_celltype_path = f"{Manuscript_RESULT}/DEGs_Celltype.xlsx"
obs = pd.read_csv(obs_path,index_col=0)
adata = sc.read(f'{WORKFLOW_DATA}/{GEX_Cohort}/gex_qc.h5ad')
downsampling_h5ad_path = f"{Manuscript_RESULT}/downsampling.h5ad"
degs_h5ad_path = f"{Manuscript_RESULT}/gex_qc_markers.h5ad"
# keep cells with annotation
adata = adata[obs.index,:]
# append meta information to the adata.obs data frame
for c in ['Cellstate','Celltype','BestResponse','Patient','Timepoint','Sample_Short','Treatment_Arm','RCB']:
    adata.obs[c] = obs.loc[adata.obs.index,c]
# remove uncharacterized cells
adata = adata[~adata.obs.Cellstate.isin(['Stromal','Immune']),:]

def downsampling(adata,celltype_col,sample_col,total_n=10000):
    prop_sample = adata.obs[sample_col].value_counts(normalize=True).to_dict()
    selected_barcodes = []
    for sample_name,obs in adata.obs.groupby(sample_col):
        prop_celltype = obs[celltype_col].value_counts(normalize=True).to_dict()
        for celltype,df in obs.groupby(celltype_col):
            n_cells = int(prop_celltype[celltype] * prop_sample[sample_name] * total_n)
            if n_cells ==0:
                print(prop_celltype[celltype])
                print(f'{sample_name} has 0 {celltype} if downsampling to {total_n} with proportional method.')
            barcodes = df.sample(n=n_cells, random_state=1,replace=False).index.tolist()
            selected_barcodes = selected_barcodes + barcodes
    print(f'Proportionally selected {len(selected_barcodes):,} from {adata.shape[0]:,} cells.')
    return adata[selected_barcodes,:].copy()

sub_adata = downsampling(adata,celltype_col='Cellstate',sample_col='Sample_Short',total_n=5000)
degs_dict = {}
for cellstate in sub_adata.obs.Cellstate.unique():
    if cellstate in ['Immune','Stromal']:
        continue
    sheet_name = cellstate.replace('/','_')  if '/' in cellstate else cellstate
    degs= pd.read_excel(degs_path,sheet_name=sheet_name,index_col=0)
    degs = degs.index[(degs.padj<0.05) & (degs.log2FoldChange>1)].tolist()
    degs_dict[cellstate] = degs

sub_adata.uns['Cellstate_Markers'] = degs_dict 
adata.uns['Cellstate_Markers'] = degs_dict 

degs_celltype_dict = {}
for celltype in sub_adata.obs.Celltype.unique():
    if celltype in ['Immune','Stromal']:
        continue
    degs= pd.read_excel(degs_celltype_path,sheet_name=celltype,index_col=0)
    degs = degs.index[(degs.padj<0.05) & (degs.log2FoldChange>1)].tolist()
    degs_celltype_dict[celltype] = degs

sub_adata.uns['Celltype_Markers'] = degs_celltype_dict
adata.uns['Celltype_Markers'] = degs_celltype_dict

# store result
sub_adata.write(downsampling_h5ad_path)
adata.write(degs_h5ad_path)