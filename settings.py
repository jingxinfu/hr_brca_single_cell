import logging
import colorlog
from pathlib import Path

# ---------------------------------------- Output Setting----------------------------------------------------------------------
HERE = Path(__file__).parent.absolute()
PlotStyle = HERE / 'paper.mplstyle'

def init_logger(dunder_name, testing_mode=False) -> logging.Logger:
    log_format = (
        '%(asctime)s - '
        '%(name)s - '
        '%(funcName)s - '
        '%(levelname)s - '
        '%(message)s'
    )
    bold_seq = '\033[1m'
    colorlog_format = f'{bold_seq} ' '%(log_color)s ' f'{log_format}'
    colorlog.basicConfig(format=colorlog_format)
    logger = logging.getLogger(dunder_name)

    if testing_mode:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    # # Output full log
    # fh = logging.FileHandler('app.log')
    # fh.setLevel(logging.DEBUG)
    # formatter = logging.Formatter(log_format)
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)

    # # Output warning log
    # fh = logging.FileHandler('app.warning.log')
    # fh.setLevel(logging.WARNING)
    # formatter = logging.Formatter(log_format)
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)

    # # Output error log
    # fh = logging.FileHandler('app.error.log')
    # fh.setLevel(logging.ERROR)
    # formatter = logging.Formatter(log_format)
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)

    return logger


# ---------------------------------------- Palette Setting----------------------------------------------------------------------
import palettable

VIVID_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
COLOR_PAlETTE = {
    'Treatment_Arm': {"Chemo->Combo": '#8d99ae', "ICI->Combo": "#fed9b7"},
    'RCB': {
        0: '#006d77',
        '0': '#006d77',
        'I': '#83c5be',
        'II': '#ffddd2',
        'III': '#e29578',
    },
    'stage':{
        'II':'pink',
        'III':'plum'
    },
    'BestResponse': {'favorable response\n(RCB 0-I)': '#006d77',
                     'unfavorable response\n(RCB II-III)': '#e29578'},
    'Timepoint': {
        'Baseline': "#352208",
        "W3D1": "#E1BB80",
        "W7D1": "#6c8152",
        "Surg+AC": "#685634",
        "AfterSurg": "#806443",
        "EOT?": "gray",
    },
    'Tech':{"5' v2 gene expression":'skyblue','Multiome ATAC + Gene Expression':'steelblue'},
    'Copy number alteration':{
        'CN=0':'gray',
        'Low-level amplification':'salmon',
        'High-level amplification':'darkred',
        'Low-level deletion':'lightsteelblue',
        'High-level deletion':'darkblue'
    },
    'Mutation type': {
        'Missense': VIVID_10[5],
        'Nonsense': VIVID_10[0],
        'In frame indel': VIVID_10[1],
        'Frameshift indel': VIVID_10[4],
        'Splice site': VIVID_10[9],
    },
    'PAM50':{'LumB':'goldenrod', 'LumA':'orange', 'Normal':'green', 'Basal':'salmon'},
    'WES_Profile':{'Y': 'black', 'N': 'white'},
    'BulkRNA_Profile':{'Y': 'black', 'N': 'white'},
    'WES_APOBEC_Enriched': {'yes': 'salmon', 'no': 'skyblue'},
    'WES_absolute_purity':{'color':'Purples'},
    'TMB': {'WES_Nonsynonymous_TMB': 'purple', 'WES_Synonymous_TMB': 'pink'},
    'Mutation clonality': {'WES_Clonal_TMB': 'purple', 'WES_Subclonal_TMB': 'pink'},
    "Mutation_Signature": {
        'WES_defective DNA mismatch repair': '#d8e2dc',
        'WES_APOBEC Cytidine Deaminase (C>T)': '#f0a6ca',
        'WES_Defects in DNA-DSB repair by HR': '#b8bedd',
    },
    "WGD": {"Yes": 'black', "No": 'white'},
    "bulkRNA_PAM50": {
        'Normal': 'steelblue',
        'Basal': 'orange',
        'LumB': 'yellow',
        'LumA': 'lightgreen',
        'Her2': 'purple',
    },
    "scSubtype": {
        'Basal_SC': 'orange',
        'LumB_SC': 'yellow',
        'LumA_SC': 'lightgreen',
        'Her2E_SC': 'purple',
    },
    "Compartment":{
      'Immune':'#DBAD6A',
      'Stromal':'#D0CE7C',
      'Epithelial':'#628395'
    },
    "Lineage":{
        'T':'orange',
        'B':'steelblue',
        'Plasma':'skyblue',
        'Myeloid':'purple',
        'Immune':'coral',
        'CAF':'darkgreen',
        'Pericyte':'lightgreen',
        'SMC':'tan',
        'Endothelial':'wheat',
        'Adipocytes':'lightcyan',
        'Stromal':'teal',
        'Epithelial':'rosybrown',
        'Tumor':'darkred',

    },
    "Celltype":{
        'CD8T':'gold',
        'CD4T':'khaki',
        'Macs':'plum',
        'Myeloid':'purple',
        'CAF':'darkgreen',
        'Endothelial':'wheat',
        'Tumor':'darkred',
        'B':'steelblue',
        'Plasma':'skyblue',
        'Adipocytes':'gray',
        'Pericyte':'green',
        'Epithelial':'salmon',
    },
    "Cellstate":{
        'NS':'white',
        'Tumor.Interferon/MHCII(I)':'green',
        'Tumor.EMT-I':'lightblue','Tumor.EMT-II':'skyblue','Tumor.EMT-III':'steelblue',
        'Tumor.ER-I':'orange','Tumor.ER-II':'red', 
        'Tumor.Cell_Cycle':'purple',
        'Tumor.Stress':'peru','Tumor.Apelin':'slategray',
        
        ## Leiden
        'Macs.APOE High':"#7209b7",'Macs.HDAC9 High':"#4a4e69",'Macs.PLCG2 High':"#9a8c98",'Macs.ZNF331 High':"#c9ada7",
        'CD8T.em.GZMK+':"#faf0ca", 'CD8T.ex.TCF7+':'#f4d35e','CD8T.m.IL7R High':"#ee964b", 'CD8T.m.FTL High':"#f95738",
        'CAF.Inf':'#678d58','CAF.MHCI High':'#74d3ae','CAF.Myo':'#a6c48a',
        'Endo.ACKR1+':'#a68a64','Endo.CXCL12+':'#d6ccc2','Endo.MHCII High':'#7f4f24', 'Endo.PRKG1+':'#582f0e',
        ## MetaPrograms
        'Macs.Lipid-IGF1':'purple',
        'Macs.Cholinergic':'red',
        'Macs.Adipogenesis':'green',
        'Macs.MHCII':'steelblue',
        'Macs.Lipid-APOE':'black',
        'Macs.Interferon':'darkorange',
        
        'Macs.Lipid':'purple',
        'Macs.Adhesion':'gray',
        'Macs.Endocytosis':'green',
        'Macs.Presentation':'steelblue',
        'Macs.Secretion':'gold',
        
        'CD8T.Cytotoxic':'#faf0ca',
        'CD8T.Dysfunction':'#f4d35e',
        'CD8T.Naive':'steelblue',
        
        'Endo.Allograft_Rejection':'#7f4f24',
        'Endo.Endo1':'#a68a64',
        'Endo.Endo2':'#d6ccc2',
        'Endo.Endo5':'#582f0e',
        'Endo.HEV1':'salmon',
        'Endo.Notch-signaling':'black',
        
        'Endothelial.Allograft_Rejection':'#7f4f24',
        'Endothelial.Endo1':'#a68a64',
        'Endothelial.Endo2':'#d6ccc2',
        'Endothelial.Endo5':'#582f0e',
        'Endothelial.HEV1':'salmon',
        'Endothelial.Notch-signaling':'black',
        'Endothelial.Fatty_Acid':'olive',
        'Endothelial.Collagen':'darkgreen',
        
        'CAF.Allograft_Rejection':'#7f4f24',
        'CAF.CAF1':'#a68a64',
        'CAF.CAF2':'#d6ccc2',
        'CAF.CAF3':'#582f0e',
        'CAF.Complement':'gray',
        'CAF.Interferon':'salmon',
        'CAF.KRAS_Signaling_up':'steelblue',
        'CAF.Pericyte-like':'green',
        'CAF.Stress':'black'
    }
}

## Update N/A category
for _, v in COLOR_PAlETTE.items():
    v.update({'N/A': '#e5e5e5'})
