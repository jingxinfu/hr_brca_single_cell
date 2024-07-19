import subprocess
import scPipe as sp
import scanpy as sc
import pandas as pd
GEX_Cohort='GEX_CCG1112_LowMt'
REPO = '.'
WORKFLOW_DATA = f'{REPO}/data/workflow'
EXTERNAL_DATA = f'{REPO}/data/external'
RESULT_DATA = f'{REPO}/data/result/cleaned_files'


def scSHC_clustering(group,subset_name,subset_var):
    """
    Perform scSHC on `subset_var` of the `group` 
    by selected `subset_var` on the `subset_name` in the `group` annotation
    """
    ## Load Data
    adata = sc.read(f'{WORKFLOW_DATA}/{GEX_Cohort}/gex_qc.h5ad')

    adata.obs.drop(adata.obs.columns,axis=1,inplace=True)
    ## Attach selected group annotation
    obs = pd.read_csv(f'{RESULT_DATA}/Annotation_{group}.csv',
                     index_col=0,low_memory=False)
    ## subset immune cells
    adata = adata[obs.index,:]
    adata.obs = obs
    assert adata.obs[subset_name].isna().sum()==0,'Duplicated cell name in obs.[annotation]'

    print(adata.obs[subset_name].value_counts())
    adata = adata[adata.obs[subset_name]==subset_var,:]

    print(f"Detect N={adata.shape[0]:,} {subset_var} cells with high quality GEX profiles;")
    print(f"max(% of mitocondrial reads): {adata.obs['pct_counts_mito'].max():.0f}")

    sp.utils.write_10X_h5(filename=f"{RESULT_DATA}/RawCount/{subset_var}.h5",
                         matrix=adata.X.toarray().T,
                         features=adata.var.index,
                         barcodes=adata.obs.index,
                         datatype='GEX')
    adata.obs['Sample'].to_frame().to_csv(f"{RESULT_DATA}/scSHC/{subset_var}_batch.csv")

    with open(f'{RESULT_DATA}/scSHC/{subset_var}.R','w') as f:
        f.write(
    f"""library(Seurat)
library(scSHC)
data <- Read10X_h5('{RESULT_DATA}/RawCount/{subset_var}.h5')
batch <- read.csv('{RESULT_DATA}/scSHC/{subset_var}_batch.csv')[,'Sample']
clusters <- scSHC(data,cores=20,alpha=0.05,batch=batch)

# alpha, which controls the family-wise error rate (default 0.05). 
# If the goal is discovery, consider setting a more lenient alpha, such as 0.25.
write.csv(clusters[[1]],'{RESULT_DATA}/scSHC/{subset_var}_clusters.csv')

pdf('{RESULT_DATA}/scSHC/{subset_var}.pdf',width=6, height=6)
par(bty='l',las=1,lwd=2,cex=0.7,oma = c(5, 1, 1, 1), mar = c(5, 5, 5,5))
plot(as.dendrogram(clusters[[2]]),center=T)
title('{subset_var}')
dev.off()

    """)
    with open(f'{RESULT_DATA}/scSHC/{subset_var}.R','r') as f:
        print(f.read())
    subprocess.run(f"Rscript {RESULT_DATA}/scSHC/{subset_var}.R",shell=True)

# Step1:
#scSHC_clustering(group='Compartment',subset_name='Compartment',subset_var='Immune')
# Step2:
# After manually annotate the scSHC clusters of the Immune compartment
scSHC_clustering(group='Immune',subset_name='Lineage',subset_var='T')
scSHC_clustering(group='Immune',subset_name='Lineage',subset_var='Myeloid')
# scSHC_clustering(group='Compartment',subset_name='Compartment',subset_var='Stromal')


