Manuscript_RESULT = f'data/result/manuscript_table/'
import pandas as pd
import os

## Get Fastq Path
fastq_path = pd.read_csv(f"{Manuscript_RESULT}/dbGaP/sample.tsv",sep='\t',index_col=0)[['Sample','Fastq_Folder']]
dbGaP_fastq = []
sn_S_map = {
    'CCG1112_03_T1_A1_SN_5GEX':'S1',
    'CCG1112_03_T2_A1_SN_5GEX':'S2',
    'CCG1112_08_T1_A1_SN_5GEX':'S3',
    'CCG1112_08_T2_A1_SN_5GEX':'S4',
    'CCG1112_08_T3_A1_SN_5GEX':'S5',
    'CCG1112_13_T1_A1_SN_5GEX':'S6',
    'CCG1112_13_T2_A1_SN_5GEX':'S7',
    'CCG1112_13_T3_A1_SN_5GEX':'S8'
}
for _,row in fastq_path.iterrows():
    if '[' in row['Fastq_Folder']:
        fastq_list = row['Fastq_Folder'].replace('["','').replace('"]','').split('","')
        for fastq in fastq_list:
            t_row = pd.Series({
                'SAMPLE_ID':'P'+'.'.join(row['Sample'].split('_')[1:3]),
                'ANALYTE_TYPE':'DNA' if 'atac' in fastq else 'RNA',
                'FASTQ_PATH':fastq
            })
            dbGaP_fastq.append(t_row)
    else:
        sample_name = os.path.basename(row['Fastq_Folder'])
        s_suffix = sn_S_map[sample_name]
        for suffix in ['L001_R1_001.fastq.gz','L001_R2_001.fastq.gz',
                      'L002_R1_001.fastq.gz','L002_R2_001.fastq.gz']:
            t_row = pd.Series({
                'SAMPLE_ID':'P'+'.'.join(row['Sample'].split('_')[1:3]),
                'ANALYTE_TYPE':'RNA',
                'FASTQ_PATH':row['Fastq_Folder']+'/'+sample_name+f"_{s_suffix}_{suffix}"
            })
            dbGaP_fastq.append(t_row)
dbGaP_fastq = pd.concat(dbGaP_fastq,axis=1).T.sort_values('SAMPLE_ID')

## Output files
obs = pd.read_csv(f"{Manuscript_RESULT}/GEX_OBS.csv",index_col=0)
meta = (obs.
        drop_duplicates(['Sample_Short']).
        sort_values(['Patient']).
        rename(columns={'Patient':'SUBJECT_ID','Sample_Short':'SAMPLE_ID',
                        'age':'AGE_ONSET','race':'RACE','sex':'SEX','ethnic':'ETHNIC','stage':'STAGE',
                        'histologic_grade':'TUMOR_GRADE','histology':'HISTOLOGY',
                        'Timepoint':'TREATMENT_TIMEPOINT','Treatment_Arm':'TREATMENT_ARM'})
       )
meta = dbGaP_fastq.merge(meta,on='SAMPLE_ID',how='left')
## 2ds
dbGaP_2ds = meta[['SUBJECT_ID']].copy()
dbGaP_2ds['CONSENT']=2
dbGaP_2ds['SEX'] = 2
dbGaP_2ds.to_csv(f"{Manuscript_RESULT}/dbGaP/dbGaP_2ds.txt",sep='\t',index=False)
## 3ds
meta[['SUBJECT_ID','SAMPLE_ID']].to_csv(f"{Manuscript_RESULT}/dbGaP/dbGaP_3ds.txt",sep='\t',index=False)
## 5ds
meta[['SUBJECT_ID','AGE_ONSET','RACE','SEX','ETHNIC','STAGE']].to_csv(f"{Manuscript_RESULT}/dbGaP/dbGaP_5ds.txt",sep='\t',index=False)

## 6ds
dbGaP_6ds=meta[['SUBJECT_ID','SAMPLE_ID','TUMOR_GRADE','STAGE','HISTOLOGY','TREATMENT_TIMEPOINT','TREATMENT_ARM','ANALYTE_TYPE']].copy()
dbGaP_6ds['PRIMARY_TUMOR_LOCATION']='Breast'
dbGaP_6ds['PRIMARY_METASTATIC_TUMOR']='Primary'
dbGaP_6ds['BODY_SITE']='Breast'
dbGaP_6ds['IS_TUMOR']='Y'
dbGaP_6ds.rename(columns={'STAGE':'TUMOR_STAGE'},inplace=True)
dbGaP_6ds[['SUBJECT_ID',
           'SAMPLE_ID',
           'ANALYTE_TYPE',
           'PRIMARY_TUMOR_LOCATION',
           'PRIMARY_METASTATIC_TUMOR',
           'BODY_SITE',
           'IS_TUMOR',
           'TUMOR_GRADE',
           'TUMOR_STAGE',
           'HISTOLOGY',
           'TREATMENT_TIMEPOINT',
           'TREATMENT_ARM']].to_csv(f"{Manuscript_RESULT}/dbGaP/dbGaP_6ds.txt",sep='\t',index=False)
## fastq path
meta[['SUBJECT_ID','SAMPLE_ID','FASTQ_PATH']].to_csv(f"{Manuscript_RESULT}/dbGaP/dbGaP_FASTQ.txt",sep='\t',index=False)