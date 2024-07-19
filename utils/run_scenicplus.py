import os
import glob
import dill
import pyranges
import warnings
warnings.filterwarnings("ignore")
import sys
import scanpy as sc
import scPipe as sp
import pickle
import pandas as pd

# scenicplus
import numpy as np
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus
import pybiomart as pbm

ATAC_Cohort='CCG1112_ATAC_QC'
GEX_Cohort='GEX_CCG1112_LowMt'
REPO = '/home/analysis/hr_brca_mo_analysis'
WORKFLOW_DATA = f'{REPO}/data/workflow'
EXTERNAL_DATA = f'{REPO}/data/external'
RESULT_DATA = f'{REPO}/data/result/cleaned_files'
Group_Variable = 'Celltype'
# Set stderr to null to avoid strange messages from ray
_stderr = sys.stderr
null = open(os.devnull,'wb')
tmp_dir = '/mnt/disks/brca_16_466_data/Test/'

adata = sc.read(f'{WORKFLOW_DATA}/{GEX_Cohort}/gex_qc.h5ad')
sc.pp.filter_genes(adata,min_cells=3) # consistent with the clustering method
for pickle_file in glob.glob(f'{WORKFLOW_DATA}/{GEX_Cohort}/{GEX_Cohort}/leiden_2000variable_genes/*.pickle'):
    if 'colors' in pickle_file:
        continue
    sp.utils.attach_attr_from_pickle(data=adata,pickle_file=pickle_file)
    
cell_data = pd.read_csv(f'{RESULT_DATA}/Celltype.csv',index_col=0)
cell_data['atac_index']= cell_data.index.map(lambda x:x.split('_')[0]+'___') + cell_data['Sample']
adata = adata[cell_data.index,:]
adata.obs = adata.obs.merge(cell_data.drop(['Sample'],axis=1),left_index=True,right_index=True)
adata.obs.set_index('atac_index',inplace=True)

cistopic_obj = dill.load(open(os.path.join(RESULT_DATA, 'mo_cistopic_obj_with_model.pickle'), 'rb'))


menr = dill.load(open(os.path.join(RESULT_DATA, 'motifs/menr.pkl'), 'rb'))


######################################################## Create  SCENIC object ########################################################


scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
######################################################## Select correct enseml version ########################################################
ensembl_version_dict = {'105': 'http://www.ensembl.org',
                        '104': 'http://may2021.archive.ensembl.org/',
                        '103': 'http://feb2021.archive.ensembl.org/',
                        '102': 'http://nov2020.archive.ensembl.org/',
                        '101': 'http://aug2020.archive.ensembl.org/',
                        '100': 'http://apr2020.archive.ensembl.org/',
                        '99': 'http://jan2020.archive.ensembl.org/',
                        '98': 'http://sep2019.archive.ensembl.org/',
                        '97': 'http://jul2019.archive.ensembl.org/',
                        '96': 'http://apr2019.archive.ensembl.org/',
                        '95': 'http://jan2019.archive.ensembl.org/',
                        '94': 'http://oct2018.archive.ensembl.org/',
                        '93': 'http://jul2018.archive.ensembl.org/',
                        '92': 'http://apr2018.archive.ensembl.org/',
                        '91': 'http://dec2017.archive.ensembl.org/',
                        '90': 'http://aug2017.archive.ensembl.org/',
                        '89': 'http://may2017.archive.ensembl.org/',
                        '88': 'http://mar2017.archive.ensembl.org/',
                        '87': 'http://dec2016.archive.ensembl.org/',
                        '86': 'http://oct2016.archive.ensembl.org/',
                        '80': 'http://may2015.archive.ensembl.org/',
                        '77': 'http://oct2014.archive.ensembl.org/',
                        '75': 'http://feb2014.archive.ensembl.org/',
                        '54': 'http://may2009.archive.ensembl.org/'}


def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
    except:
        print('Host not reachable')
v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")
biomart_host = ensembl_version_dict[v]

######################################################## RUN SCENICPLUS ########################################################
scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]
outfolder = os.path.join(RESULT_DATA, 'scenicplus')
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = [f'GEX_{Group_Variable}'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = f'{EXTERNAL_DATA}/TF_db/TF_names_v_1.01.txt',
        save_path = outfolder,
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = f'{EXTERNAL_DATA}/TF_db',
        n_cpu = 12,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(f'{outfolder}/scplus_obj.pkl', 'wb'), protocol=-1)
    raise(e)
