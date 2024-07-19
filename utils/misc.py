from typing import Tuple,List,Union,Dict
from pathlib import Path
import glob
import pandas as pd
import numpy as np
import os
def mergeRsemGeneOutput(folder_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Merge single sample's rsem gene quantification files into cohort level gene by sample data frame

    Parameters
    ----------
    folder_path : str
        Path to the cohort folder, where all salmon quantification output files locate

    Returns
    -------
    (pd.DataFrame,pd.DataFrame)
        1. gene (Hugo Symbol) by sample Raw Count data frame
        2. gene (Hugo Symbol) by sample TPM data frame

    """
    file_suffix = ".genes.results.gz"
    tpms = []
    counts = []
    if isinstance(folder_path,str):
        folder_path = Path(folder_path)
    for file_path in glob.glob(str(folder_path / "*")):
        if not file_path.endswith(file_suffix):
            continue
        sample_name = os.path.basename(file_path).replace(file_suffix,'')
        
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        tpm = df['TPM']
        tpm.name = sample_name
        tpms.append(tpm)

        count = df["expected_count"]
        count.name = sample_name
        counts.append(count.astype(int))
    return (
        pd.concat(counts, axis=1),
        pd.concat(tpms, axis=1),
    )
### FOR DEGs
import anndata
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
import scanpy as sc
def oneToRestDEGs(pdata:anndata.AnnData,
                  group:str,
                  confounders:Union[List[str],None]=None,
                  output_xlsx:Union[Path,None]=None)->Dict[str,pd.DataFrame]:
    
    pdata_test = pdata.copy()
    if pdata_test.obs[group].nunique() <2:
        return pd.DataFrame()
    pdata_test.obs = pdata_test.obs.merge(pd.get_dummies(pdata_test.obs[group]).replace({1:'this',0:'rest'}),
                            left_index=True,right_index=True)
    results = {}
    for group in pdata_test.obs[group].unique():
        # Obtain genes that pass the thresholds
        genes = dc.filter_by_expr(pdata_test, group=group, min_count=10, min_total_count=15)
        # Filter by these genes
        pdata_test_filtered = pdata_test[:, genes].copy()
        # rename the group to avoid "-" in the name
        pdata_test_filtered.obs.rename(columns={group:'group'},inplace=True)
        # Build DESeq2 object
        inference = DefaultInference(n_cpus=8)
        dds = DeseqDataSet(
            adata=pdata_test_filtered,
            design_factors=['group'] +confounders if isinstance(confounders,list) else ['group'],
            refit_cooks=True,
            inference=inference,
        )
        # Compute LFCs
        dds.deseq2()
        
        # Running Wald test
        stat_res = DeseqStats(
            dds,
            contrast=['group', 'this', 'rest'],
            inference=inference,
        )
        stat_res.summary()
        results[group] = stat_res.results_df
        
    ## Store results if path is provided
    if output_xlsx is not None:
        with pd.ExcelWriter(output_xlsx) as f:
            for k,v in results.items():
                v.to_excel(f,sheet_name=k.replace('/','_'),index=True)
    return results
    
def get_signature_matrix(pdata:anndata.AnnData,
                         degs_results:Dict[str,pd.DataFrame],
                         group:str,
                         min_genes:int=300,
                         max_genes:int=500,
                         n_gene_steps:int=25,
                         exclude_groups:Union[List[str],None]=None)->pd.DataFrame:
    cond_number = np.inf
    optimal_sig_matrix = None
    # using CPM
    pdata_test = pdata.copy()
    sc.pp.normalize_total(pdata_test,target_sum=1e6)
    # pick the signature matrix that are most stable defined by lowest condition number.
    for n in np.arange(min_genes,max_genes,n_gene_steps):
        selected_genes = []
        for k,v in degs_results.items():
            if k in  exclude_groups:
                continue
            genes = v.sort_values(by='stat',ascending=False).index[(v.padj<0.05)&(v.log2FoldChange>0)].tolist()[:n]
            selected_genes = selected_genes + genes
        selected_genes = list(set(selected_genes))
        sig_matrix = sc.get.obs_df(pdata_test[~pdata_test.obs.Celltype.isin(exclude_groups),:],selected_genes+[group])
        sig_matrix = sig_matrix.groupby([group]).mean().T
        curr_cond_number = np.linalg.cond(sig_matrix.values,p=None) # 2-norm, computed directly using the SVD
        if curr_cond_number < cond_number:
            cond_number = curr_cond_number
            optimal_sig_matrix  = sig_matrix
    print(f'Find optimal signature matrix with {optimal_sig_matrix.shape[0]} genes, condition number = {cond_number}')
    return optimal_sig_matrix