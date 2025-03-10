#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 15:06:07 2023

@author: aaronatkinson
"""

#import packages
import pandas as pd 
import numpy as np
import glob, os, shutil, argparse
#import re


print('\n\n|||||||||||||||||||||||||||||||||||||||||||||| SUPLEMENTING META DATA XLSX FILE FOR cBIO ||||||||||||||||||||||||||||||||||||||||||||||\n\n')
# Find all Excel files in the current directory
excel_files = glob.glob("*.xlsx")

# Check if any files were found
if not excel_files:
    raise FileNotFoundError("No .xlsx files found in the directory.")

# Load the first (or latest) Excel file found
excel_file = excel_files[0]  
print(f"Loading: {excel_file}")

# Read the Excel file with engine specified and build temp txt file for further parsing. Found missing diagnosis among data creating the temp txt file resolves this. 
meta1 = pd.read_excel(excel_file, sheet_name=0, engine="openpyxl")  # Use 'openpyxl' for .xlsx files

meta1 = meta1.loc[~meta1['Load'].isin(['no'])]
meta1.drop(['Comment (do not load to cBio_', 'Study ID', 'Load'], axis = 1, inplace=True)
meta1.to_csv('Temp.txt', sep='\t', index=False)
meta = pd.read_csv('Temp.txt', sep='\t')
#Change Case and fill spaces
meta.columns = meta.columns.astype(str).str.replace(' ', '_')
meta.rename(columns = str.upper , inplace = True)
# strip white space from name columns
meta['PATIENT_ID'] = meta['PATIENT_ID'].str.strip()
#Build Patient list from lab meta data
patient = meta[['PATIENT_ID', 'GENDER']]
patient_list = patient['PATIENT_ID'].tolist()

def reorder_columns(df, columns):
    # Filter the columns that are present in the DataFrame
    present_columns = [col for col in columns if col in df.columns]
    # Add any additional columns that are not in the specified list
    additional_columns = [col for col in df.columns if col not in present_columns]
    # Reorder the DataFrame columns
    return df[present_columns + additional_columns]

# Define the columns to be used
columns = ['PATIENT_ID', 'SAMPLE_ID', 'SAMPLE_CLASS', 'PRIMARY_OR_RECURRENCE', 'WGS_ID', 'WGS_PT_NORMAL_ID', 'WES_ID', 'WES_PT_NORMAL_ID', 'RNASEQ', 'RNASEQ_META', 'RNASEQ_LIBRARY', 'CNV_(SEGFILESAMPLEIDS)', 'CNV_PT_NORMAL_(SEGFILESAMPLEIDS)2', 'CNV_PLATFORM/CHIP', 'AGE_AT_TISSUE_COLLECTION', 'ALT_ID', 'TOWARDS_ID', 'CC', 'CORE_SAMPLE_ID', 'COLLECTION_SITE', 'DISEASE', 'ER_PATIENT', 'PR_PATIENT', 'HER2_PATIENT', 'ER_PDX', 'PR_PDX', 'HER2_PDX', 'HCI_ID', 'INJECTION_SIDE', 'ORGANOID_PREP', 'DNA_EXTRACTION_', 'MOUSE_PASSAGE', 'ORGANOID_PASSAGE', 'CNV_PASSAGE', 'PASSAGE_DAY', 'CNV_DAY', 'SAMPLE_TYPE', 'SHADOW_ID', 'CORE_PATIENT_ID', 'TX_NAIVE_OR_PRETREATED', 'STUDY', 'EXPERIMENT', 'PDX_GROW_NOGROW', 'DATE_OF_COLLECTION_(ONLY_SAMPLES_W-O_SHADOWID)']

sample = meta
del meta
sample.drop(['GENDER'], axis = 1, inplace=True)
print(*sample)
# strip white space from name columns
sample['SAMPLE_ID'] = sample['SAMPLE_ID'].str.strip()
sample.rename(columns={'CC_NUMBER': 'CC', 'PRIMARY_-_RECURRENCE': 'PRIMARY_OR_RECURRENCE'}, inplace=True)
sample['CC'] = sample['CC'].str.strip()
#Build supplemental sample file
sample['SAMPLE_TYPE'] = sample['SAMPLE_TYPE'].str.upper()
sample['SAMPLE_CLASS'] = sample['SAMPLE_TYPE']

sample = sample.replace({'SAMPLE_CLASS' : {'PDXEV': 'xenograft','PATIENT_TUMOR': 'Biopsy', 'PDXOX': 'xenograft', 'PDXO': 'Organoid','PDOX': 'xenograft','CELLS': 'CellLine'}})
sample = sample.replace({'SAMPLE_TYPE' : {'PDXEV': 'PDXEv','PATIENT_TUMOR': 'Patient_tumor', 'PDXOX': 'PDXoX', 'PDXO': 'PDxO','PDOX': 'PDoX','CELLS': 'Cells'}})
sample['DATE_OF_COLLECTION_(ONLY_SAMPLES_W-O_SHADOWID)'] = ''
# Reorder columns for sample and sample2
sample = reorder_columns(sample, columns)
#create universal sample ID between CC and CoreID
sample['CC'] = sample['CC'].fillna(sample['CORE_SAMPLE_ID'])
sample2 = reorder_columns(sample, columns)
#Found blank lines
sample2['PATIENT_ID'] = sample2['PATIENT_ID'].fillna('REMOVE')
sample2 = sample2.loc[~sample2['PATIENT_ID'].isin(['REMOVE'])]

#Build HCI_ID heatmap legend from Sample info to group samples by patient in custom assay
HCI = sample2[['SAMPLE_ID']]
HCI['HCI_ID'] = HCI['SAMPLE_ID'].str.split('_', expand=True)[0]
HCI = HCI[HCI['HCI_ID'].astype(str).str.contains('HCI', case=False)]
# Transpose the dataframe
HCI_transposed = HCI.set_index('SAMPLE_ID').T
HCI_merged = pd.concat([HCI, HCI_transposed], axis=1)

# Fill down the values for the HCI_transposed columns
HCI_merged.fillna(method='ffill', inplace=True)
HCI_merged.fillna(method='bfill', inplace=True)

HCI_merged.insert(0, column='ENTITY_STABLE_ID', value=HCI['HCI_ID'])
HCI_merged.insert(3, column='DESCRIPTION', value='Samples grouped by patient/HCI_ID')
HCI_merged.drop(['SAMPLE_ID'], axis = 1, inplace=True)
HCI_merged.rename(columns={'HCI_ID': 'NAME'}, inplace=True)
for column in HCI_merged.columns.difference(['ENTITY_STABLE_ID', 'NAME', 'DESCRIPTION']):
   HCI_merged[column] = HCI_merged[column] == HCI_merged['NAME']
   HCI_merged[column] = HCI_merged[column].replace({True: '1', False: '0'})
HCI_merged.drop_duplicates(subset=['ENTITY_STABLE_ID'], keep='first', inplace=True)
HCI_merged.to_csv('HCI.txt', sep='\t', index=False)
currentDateTime = pd.Timestamp.now().strftime("%d%b%Y")
str_currentDateTime = str(currentDateTime)
file_name = str_currentDateTime+"-PATIENT_ONCOPRINT_LEGEND.txt"
HCI_merged.to_csv(file_name, sep='\t', index=False)
del file_name

#Build samplesheet for molecular data
MolSample = sample2[['SAMPLE_ID', 'WGS_ID', 'WGS_PT_NORMAL_ID', 'WES_ID', 'WES_PT_NORMAL_ID', 'RNASEQ', 'RNASEQ_META', 'RNASEQ_LIBRARY', 'CNV_(SEGFILESAMPLEIDS)']]
file_name = str_currentDateTime+"-MolecularSampleSheet.txt"
MolSample.to_csv(file_name, sep='\t', index=False)
del file_name

#Take samples for case lists now before headers change etc. 
cases_all_list = sample2['SAMPLE_ID'].tolist()
# Count the number of strings in the list
cases_all = len(cases_all_list)
# Convert the list to a tab-separated string
cases_all_delim = '\t'.join(map(str, cases_all_list))

#Sequenced cases
sequenced = sample2[['SAMPLE_ID', 'WES_ID', 'WES_PT_NORMAL_ID']]
# Drop rows without a value in either 'WES_ID' or 'WES_PT_NORMAL_ID'
sequenced = sequenced.dropna(subset=['WES_ID', 'WES_PT_NORMAL_ID'], how='all')
cases_sequenced_list = sequenced['SAMPLE_ID'].tolist()
# Count the number of strings in the list
cases_seq = len(cases_sequenced_list)
# Convert the list to a tab-separated string
cases_seq_delim = '\t'.join(map(str, cases_sequenced_list))

#CNV cases
CNV = sample2[['SAMPLE_ID', 'WES_ID', 'WES_PT_NORMAL_ID', 'CNV_(SEGFILESAMPLEIDS)']]
# Drop rows without a value in either 'WES_ID' or 'WES_PT_NORMAL_ID', or 'CNV_(SEGFILESAMPLEIDS)
CNV = CNV.dropna(subset=['CNV_(SEGFILESAMPLEIDS)'], how='all')
cases_cna_list = CNV['SAMPLE_ID'].tolist()
# Count the number of strings in the list
cases_cna = len(cases_cna_list)
# Convert the list to a tab-separated string
cases_cna_delim = '\t'.join(map(str, cases_cna_list))

#CNV samples with sequencing
CNV = CNV.dropna(subset=['WES_ID', 'WES_PT_NORMAL_ID'], how='all')
cases_cnaseq_list = CNV['SAMPLE_ID'].tolist()
# Count the number of strings in the list
cases_cnaseq = len(cases_cnaseq_list)
# Convert the list to a tab-separated string
cases_cnaseq_delim = '\t'.join(map(str, cases_cnaseq_list))

#RNAseq cases
RNAseq = sample2[['SAMPLE_ID', 'RNASEQ']]
# Drop rows without a value in 'RNASEQ'
RNAseq = RNAseq.dropna(subset=['RNASEQ'], how='all')
cases_rna_seq_list = RNAseq['SAMPLE_ID'].tolist()
# Count the number of strings in the list
cases_rna_seq = len(cases_rna_seq_list)
# Convert the list to a tab-separated string
cases_rna_seq_delim = '\t'.join(map(str, cases_rna_seq_list))


#Build corrected sample file
# Determine data types and replace with STRING/NUMBER
dtype_mapping = sample2.dtypes.replace({
    'object': 'STRING',
    'float64': 'NUMBER',
    'int64': 'NUMBER'
})

# Convert dtype_mapping to a single-row DataFrame
dtype_df = pd.DataFrame([dtype_mapping], columns=sample2.columns)
# Add '#' prefix to the first column in dtype_df (Row 2)
dtype_df.iloc[:, 0] = "#" + dtype_df.iloc[:, 0].astype(str)

# Create other required rows
header_dup_1 = pd.DataFrame([sample2.columns], columns=sample2.columns)  # Duplicate header (Row 1)
header_dup_2 = pd.DataFrame([sample2.columns], columns=sample2.columns)  # Duplicate header (Row 4)
ones_row = pd.DataFrame([["1"] * len(sample2.columns)], columns=sample2.columns)  # Row of "1" (Row 3)
# Add '#' prefix to the first column in header_dup_1, header_dup_2, and ones_row
header_dup_1.iloc[:, 0] = "#" + header_dup_1.iloc[:, 0].astype(str)
ones_row.iloc[:, 0] = "#" + ones_row.iloc[:, 0].astype(str)

# Combine all components
sample2 = pd.concat([header_dup_1, dtype_df, ones_row, header_dup_2, sample2], ignore_index=True)
sample2.rename(columns={'PATIENT_ID': '#PATIENT_ID'}, inplace=True)
file_name = str_currentDateTime+"-SAMPLE_ATTRIBUTES.txt"
sample2.to_csv(file_name, sep='\t', index=False)
del file_name
del header_dup_1, dtype_df, ones_row, header_dup_2

sample['CCtemp'] = sample['CC']
sample['CC'] = sample['CC'].astype(str).str.replace(r'\D+', '')
print(*sample)


print('\n\n|||||||||||||||||||||||||||||||||||||||||||||| BUILDING Treatment FILE FOR cBIO ||||||||||||||||||||||||||||||||||||||||||||||\n\n')
df = pd.read_csv('KTVPatients.txt', sep='|', encoding= 'utf-8')
print(*df)
df = df.drop([0])
df.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df = df.loc[df['PATIENT_ID'].isin(patient_list)]

df.rename(columns = str.upper , inplace = True)
df.rename(columns={'VITAL_STATUS': 'OS_STATUS'}, inplace=True)
#keep DOB if needed later
dob = df[['PATIENT_ID', 'DATE OF BIRTH']]
dob['DATE OF BIRTH'] = pd.to_datetime(dob['DATE OF BIRTH']).dt.date
dob['DATE OF BIRTH'] = pd.to_datetime(dob['DATE OF BIRTH'], format='%Y/%m/%d')
df.drop(['ERRORS', 'GENDER'], axis = 1, inplace=True)
patient = patient.merge(df, on = 'PATIENT_ID', how='left')
patient.rename(columns={'IS DECEASED?': 'OS_STATUS'}, inplace=True)
patient['OS_STATUS'] = patient['OS_STATUS'].str.replace('Y', '1:DECEASED')
patient['OS_STATUS'] = patient['OS_STATUS'].str.replace('N', '0:LIVING')
patient.to_csv("patientTemp.txt", sep='\t', index=False)


RFS = pd.read_csv('KTVRecur.txt', sep='|', encoding= 'utf-8')
print(*RFS)
RFS = RFS.drop([0])
RFS.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
RFS = RFS.loc[RFS['PATIENT_ID'].isin(patient_list)]
RFS.rename(columns = str.upper , inplace = True)
RFS = RFS.sort_values(by=['PATIENT_ID', 'RECUR DAYS AFTER DIAG'])
RFS.rename(columns={'DIAGNOSTICSER': 'RECURRENCE_CPT_CODE', 'RECUR TYPE': 'RECURRENCE_TYPE'}, inplace=True)
RFS['COURSE NUM'] = RFS['COURSE NUM'].astype(str)
##Build Rucerrence timeline file separately. 
RecurDF = RFS
RecurDF.drop(['HCI PERSON ID', 'FIRST RECUR?'], axis = 1, inplace=True)
#There are only 4 secondary recurrences. For RFS take time to first occurence
RFS = RFS.loc[RFS['COURSE NUM'].isin(['1'])]

##Need diagnosis date save RecurDF for patient merge later also follow this with OS.
RecurDF['START_DATE'] = ''
RecurDF['STOP_DATE'] = ''
RecurDF['EVENT_TYPE'] = 'STATUS'
RecurDF['RECURRENCE'] = 'RECURRENCE'
RecurDF.rename(columns={'RECURRENCE_TYPE': 'SOURCE'}, inplace=True)
RecurDF = RecurDF.merge(dob, on = 'PATIENT_ID', how = 'left')
#RecurDF['RECUR DATE'] = pd.to_datetime(RecurDF['RECUR DATE']).dt.date
RecurDF['RECUR DATE'] = pd.to_datetime(RecurDF['RECUR DATE'], format='%m/%d/%Y')

RecurDF['START_DATE'] = RecurDF['RECUR DATE'] - RecurDF['DATE OF BIRTH']
RecurDF['START_DATE'] = RecurDF['START_DATE'].astype(str)
RecurDF['START_DATE'] = RecurDF['START_DATE'].str.replace(' days', '')
RecurDF = RecurDF.loc[~RecurDF['START_DATE'].isin(['NaT'])]
RecurDF.drop(['RECUR DATE'], axis = 1, inplace=True)
RecurDF.drop(['DATE OF BIRTH'], axis = 1, inplace=True)
RecurDF = RecurDF[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SOURCE']]
RecurDF.to_csv("STATUSrecurTemp.txt", sep='\t', index=False)
###save dataframe for later append of death and diagnosis

#####
df1 = pd.read_csv('KTVTx.txt', sep='|', encoding= 'utf-8')
print(*df1)
df1 = df1.drop([0])
df1.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df1 = df1.loc[df1['PATIENT_ID'].isin(patient_list)]
df1.rename(columns = str.upper , inplace = True)
df1.drop(['HCI PERSON ID'], axis = 1, inplace=True)
df1.rename(columns={'DIAGNOSTICSER': 'CPT CODE'}, inplace=True)
df1.insert(1, column='START_DATE', value='')
df1.insert(2, column='STOP_DATE', value='')
df1.insert(3, column='EVENT_TYPE', value='TREATMENT')

df1 = df1.merge(dob, on = 'PATIENT_ID', how='left')
df1['CHEMO START DATE'] = pd.to_datetime(df1['CHEMO START DATE']).dt.date
df1['CHEMO START DATE'] = pd.to_datetime(df1['CHEMO START DATE'], format='%Y/%m/%d')
df1['CHEMO END DATE'] = pd.to_datetime(df1['CHEMO END DATE']).dt.date
df1['CHEMO END DATE'] = pd.to_datetime(df1['CHEMO END DATE'], format='%Y/%m/%d')
df1['RAD THERAPY START DATE'] = pd.to_datetime(df1['RAD THERAPY START DATE']).dt.date
df1['RAD THERAPY START DATE'] = pd.to_datetime(df1['RAD THERAPY START DATE'], format='%Y/%m/%d')
df1['RAD THERAPY END DATE'] = pd.to_datetime(df1['RAD THERAPY END DATE']).dt.date
df1['RAD THERAPY END DATE'] = pd.to_datetime(df1['RAD THERAPY END DATE'], format='%Y/%m/%d')
df1['IMMUNO START DATE'] = pd.to_datetime(df1['IMMUNO START DATE']).dt.date
df1['IMMUNO START DATE'] = pd.to_datetime(df1['IMMUNO START DATE'], format='%Y/%m/%d')
df1['HORMONE START DATE'] = pd.to_datetime(df1['HORMONE START DATE']).dt.date
df1['HORMONE START DATE'] = pd.to_datetime(df1['HORMONE START DATE'], format='%Y/%m/%d')
df1['SYSTEMIC THERAPY START DATE'] = pd.to_datetime(df1['SYSTEMIC THERAPY START DATE']).dt.date
df1['SYSTEMIC THERAPY START DATE'] = pd.to_datetime(df1['SYSTEMIC THERAPY START DATE'], format='%Y/%m/%d')
df1['FIRST SURG DATE'] = pd.to_datetime(df1['FIRST SURG DATE']).dt.date
df1['FIRST SURG DATE'] = pd.to_datetime(df1['FIRST SURG DATE'], format='%Y/%m/%d')
#chemo
chemo = df1[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'CHEMOTHERAPY TYPE', 'COURSE NUM','CPT CODE','CHEMO START DATE', 'CHEMO END DATE', 'DATE OF BIRTH']]
chemo.rename(columns={'CHEMOTHERAPY TYPE': 'SUBTYPE_AGENT'}, inplace=True)
chemo['SUBTYPE_AGENT'].fillna('remove', inplace=True)
chemo = chemo.loc[~chemo['SUBTYPE_AGENT'].isin(['remove'])]
chemo.insert(4, column='TREATMENT_TYPE', value='CHEMOTHERAPY')
chemo['START_DATE'] = (chemo['CHEMO START DATE'] - chemo['DATE OF BIRTH']).dt.days
chemo['STOP_DATE'] = chemo['CHEMO END DATE'] - chemo['DATE OF BIRTH']
chemo['START_DATE'] = chemo['START_DATE'].astype(str)
chemo['START_DATE'] = chemo['START_DATE'].str.replace(' days', '')
chemo = chemo.loc[~chemo['START_DATE'].isin(['NaT'])]
chemo['START_DATE'] = chemo['START_DATE'].str.replace('nan', '0')
chemo['STOP_DATE'] = chemo['STOP_DATE'].astype(str)
chemo['STOP_DATE'] = chemo['STOP_DATE'].str.replace(' days', '')
chemo['STOP_DATE'] = chemo['STOP_DATE'].str.replace('NaT', '')

chemo.drop(['CHEMO START DATE'], axis=1, inplace = True)
chemo.drop(['CHEMO END DATE'], axis=1, inplace = True)
chemo.to_csv("chemoTemp.txt", sep='\t', index=False)
#rad
rad = df1[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'RAD THERAPY TYPE', 'COURSE NUM','CPT CODE','RAD NUM TX', 'RAD THERAPY START DATE', 'RAD THERAPY END DATE', 'DATE OF BIRTH']]
rad.rename(columns={'RAD THERAPY TYPE': 'SUBTYPE_AGENT'}, inplace=True)
rad['SUBTYPE_AGENT'].fillna('remove', inplace=True)
rad = rad.loc[~rad['SUBTYPE_AGENT'].isin(['remove'])]
rad.insert(4, column='TREATMENT_TYPE', value='RADIATION')
rad['START_DATE'] = rad['RAD THERAPY START DATE'] - rad['DATE OF BIRTH']
rad['STOP_DATE'] = rad['RAD THERAPY END DATE'] - rad['DATE OF BIRTH']
rad['START_DATE'] = rad['START_DATE'].astype(str)
rad['START_DATE'] = rad['START_DATE'].str.replace(' days', '')
rad = rad.loc[~rad['START_DATE'].isin(['NaT'])]
rad['STOP_DATE'] = rad['STOP_DATE'].astype(str)
rad['STOP_DATE'] = rad['STOP_DATE'].str.replace(' days', '')
rad['STOP_DATE'] = rad['STOP_DATE'].str.replace('NaT', '')
rad.drop(['RAD THERAPY START DATE'], axis=1, inplace = True)
rad.drop(['RAD THERAPY END DATE'], axis=1, inplace = True)
rad.to_csv("radTemp.txt", sep='\t', index=False)
#immuno
immuno = df1[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'IMMUNO TYPE', 'COURSE NUM','CPT CODE', 'IMMUNO START DATE', 'DATE OF BIRTH']]
immuno.rename(columns={'IMMUNO TYPE': 'SUBTYPE_AGENT'}, inplace=True)
immuno['SUBTYPE_AGENT'].fillna('remove', inplace=True)
immuno = immuno.loc[~immuno['SUBTYPE_AGENT'].isin(['remove'])]
immuno.insert(4, column='TREATMENT_TYPE', value='Immunoptherapy')
immuno['START_DATE'] = immuno['IMMUNO START DATE'] - immuno['DATE OF BIRTH']
immuno['START_DATE'] = immuno['START_DATE'].astype(str)
immuno['START_DATE'] = immuno['START_DATE'].str.replace(' days', '')
immuno = immuno.loc[~immuno['START_DATE'].isin(['NaT'])]
immuno.drop(['IMMUNO START DATE'], axis=1, inplace = True)
immuno.to_csv("immunoTemp.txt", sep='\t', index=False)
#hormone
hormone = df1[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'HORMONE THERAPY TYPE', 'COURSE NUM','CPT CODE', 'HORMONE START DATE', 'DATE OF BIRTH']]
hormone.rename(columns={'HORMONE THERAPY TYPE': 'SUBTYPE_AGENT'}, inplace=True)
hormone['SUBTYPE_AGENT'].fillna('remove', inplace=True)
hormone = hormone.loc[~hormone['SUBTYPE_AGENT'].isin(['remove'])]
hormone.insert(4, column='TREATMENT_TYPE', value='HORMONE THERAPY')
hormone['START_DATE'] = hormone['HORMONE START DATE'] - hormone['DATE OF BIRTH']
hormone['START_DATE'] = hormone['START_DATE'].astype(str)
hormone['START_DATE'] = hormone['START_DATE'].str.replace(' days', '')
hormone = hormone.loc[~hormone['START_DATE'].isin(['NaT'])]
hormone.drop(['HORMONE START DATE'], axis=1, inplace = True)
hormone.to_csv("hormoneTemp.txt", sep='\t', index=False)
#systemic therapy
SYSTEMIC = df1[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SYSTEMIC THERAPY PROCEDURE', 'COURSE NUM','CPT CODE', 'SYSTEMIC THERAPY START DATE', 'DATE OF BIRTH']]
SYSTEMIC.rename(columns={'SYSTEMIC THERAPY PROCEDURE': 'SUBTYPE_AGENT'}, inplace=True)
SYSTEMIC.dropna(subset=['SYSTEMIC THERAPY START DATE'])
SYSTEMIC.insert(4, column='TREATMENT_TYPE', value='SYSTEMIC THERAPY')
SYSTEMIC['START_DATE'] = SYSTEMIC['SYSTEMIC THERAPY START DATE'] - SYSTEMIC['DATE OF BIRTH']
SYSTEMIC['START_DATE'] = SYSTEMIC['START_DATE'].astype(str)
SYSTEMIC['START_DATE'] = SYSTEMIC['START_DATE'].str.replace(' days', '')
SYSTEMIC = SYSTEMIC.loc[~SYSTEMIC['START_DATE'].isin(['NaT'])]
SYSTEMIC.drop(['SYSTEMIC THERAPY START DATE'], axis=1, inplace = True)
SYSTEMIC.to_csv("systemicTemp.txt", sep='\t', index=False)
#surgery/sammples
SURGERY = df1[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SURG PRIM SITE PROCEDURE', 'COURSE NUM','CPT CODE', 'FIRST SURG DATE', 'DATE OF BIRTH']]
#SURGERY.rename(columns={'SURGERY TYPE': 'SUBTYPE_AGENT'}, inplace=True)
SURGERY['START_DATE'] = SURGERY['FIRST SURG DATE'] - SURGERY['DATE OF BIRTH']
SURGERY['START_DATE'] = SURGERY['START_DATE'].astype(str)
SURGERY['START_DATE'] = SURGERY['START_DATE'].str.replace(' days', '')
SURGERY = SURGERY.loc[~SURGERY['START_DATE'].isin(['NaT'])]
SURGERY.drop(['FIRST SURG DATE', 'COURSE NUM', 'DATE OF BIRTH'], axis=1, inplace = True)
SURGERY['EVENT_TYPE'] = 'Surgery'

df2 = pd.read_csv('KTVDx.txt', sep='|', encoding= 'utf-8')
df2 = df2.drop([0])
df2.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df2.rename(columns = str.upper , inplace = True)

df2 = df2.loc[df2['PATIENT_ID'].isin(patient_list)]
df2 = df2.merge(dob, on = 'PATIENT_ID', how='left')
df2.insert(1, column='START_DATE', value='')
df2.insert(2, column='STOP_DATE', value='')
df2.insert(3, column='EVENT_TYPE', value='Diagnosis')
df2.drop(['HCI PERSON ID'], axis = 1, inplace=True)
df2.drop(['AGE AT DX'], axis = 1, inplace=True)
df2['DATE INITIAL DX'] = pd.to_datetime(df2['DATE INITIAL DX'])

#save secondary non-BC for later to incorporate into sample data
df2['ICD-O SITE'] = df2['ICD-O SITE'].str.replace('Nipple', 'Nipple/Breast')
df2['ICD-O SITE'] = df2['ICD-O SITE'].str.strip('  ')
df2noBC = df2[~df2['ICD-O SITE'].astype(str).str.contains('breast', case=False)]
df2 = df2[df2['ICD-O SITE'].astype(str).str.contains('breast', case=False)]
df2['COUNT'] = df2.groupby(['PATIENT_ID']).cumcount()+1
#keep first cancers
df2secondary = df2.loc[~df2['COUNT'].isin([1])]
df2secondary = pd.concat([df2noBC])
df2 = df2.loc[df2['COUNT'].isin([1])]


##separate diagnosis date dataframe for DFS and OS and later patient merge.
diagBC = df2[['PATIENT_ID', 'DATE INITIAL DX']]
df2['START_DATE'] = df2['DATE INITIAL DX'] - df2['DATE OF BIRTH']
df2['START_DATE'] = df2['START_DATE'].astype(str)
df2['START_DATE'] = df2['START_DATE'].str.replace(' days', '')
df2['START_DATE'] = df2['START_DATE'].str.replace('NaT', '0')
df2.drop(['DATE INITIAL DX'], axis = 1, inplace=True)
df2.drop(['DATE OF BIRTH'], axis = 1, inplace=True)
df2['START_DATE'].fillna(0)
df2 = df2[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'ICD-O SITE CODE', 'ICD-O SITE', 'ICD-O HISTO CODE', 'ICD-O HISTOLOGY', 'GENERAL STAGE', 'T', 'N', 'M', 'STAGE', 'STAGED FROM', 'STAGE SEER SYSTEM', 'TUMOR SIZE', '# REG NODES EXAMINED', '# POS REG NODES', 'METS AT DX?', 'METS AT DX DESC (IF KNOWN)', 'MIN RECUR DAYS AFTER DX', 'RECUR TYPE', 'HEIGHT AT DX (CM)', 'DAYS FROM DX TO HEIGHT', 'WEIGHT AT DX (KG)', 'DAYS FROM DX TO WEIGHT', 'BMI AT DX']]
df2recur = df2[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'MIN RECUR DAYS AFTER DX', 'RECUR TYPE']]

df2recur.rename(columns={'MIN RECUR DAYS AFTER DX': 'DAYS'}, inplace=True)
df2recur['DAYS'] = pd.to_numeric(df2recur['DAYS'], errors='coerce')
df2recur['DAYS'].fillna(0, inplace=True)
df2recur = df2recur[df2recur['DAYS'] > 0]

# Extract the numeric part of START_DATE and convert to int
df2recur['START_DATE'] = df2recur['START_DATE'].str.extract(r'(\d+)').astype(int)
df2recur = df2recur[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'RECUR TYPE']]
df2 = df2[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'ICD-O SITE CODE', 'ICD-O SITE', 'ICD-O HISTO CODE', 'ICD-O HISTOLOGY', 'GENERAL STAGE', 'T', 'N', 'M', 'STAGE', 'STAGED FROM', 'STAGE SEER SYSTEM', 'TUMOR SIZE', '# REG NODES EXAMINED', '# POS REG NODES', 'METS AT DX?', 'METS AT DX DESC (IF KNOWN)', 'HEIGHT AT DX (CM)', 'DAYS FROM DX TO HEIGHT', 'WEIGHT AT DX (KG)', 'DAYS FROM DX TO WEIGHT', 'BMI AT DX']]
df2 =  pd.concat([df2, df2recur], ignore_index=True )
df2 = df2.sort_values(by=['PATIENT_ID', 'START_DATE'], ascending=[True, False])
#Low/No clinical data patient filter
df2 = df2.loc[~df2['PATIENT_ID'].str.contains('BCM')]
df2 = df2.loc[~df2['PATIENT_ID'].str.contains('CW01')]
#Keep Diagnoses as separate timeline file
file_name = str_currentDateTime+"-DIAGNOSIS-TIMELINE.txt"
df2.to_csv(file_name, sep='\t', index=False)
del file_name

#SURGERY = pd.concat([df2])
SURGERY = SURGERY.loc[~SURGERY['START_DATE'].isin(['NaT'])]

#merge in pathology report on sample to verify surgical path accession and diagnosis. 
path = pd.read_csv('KTVBreastPath.txt', sep='|', encoding= 'utf-8')
path = path.drop([0])
path.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
path.rename(columns = str.upper , inplace = True)

path = path.loc[path['PATIENT_ID'].isin(patient_list)]
path = path.merge(dob, on = 'PATIENT_ID', how='left')
path.insert(1, column='START_DATE', value='')
path.insert(2, column='STOP_DATE', value='')
path.insert(3, column='EVENT_TYPE', value='Pathology Report')
path['OBSERVATION DATE'] = pd.to_datetime(path['OBSERVATION DATE'], format='%m/%d/%Y')
path['START_DATE'] = path['OBSERVATION DATE'] - path['DATE OF BIRTH']
path['START_DATE'] = path['START_DATE'].astype(str)
path['START_DATE'] = path['START_DATE'].str.replace(' days', '')
path.drop(['HCI PERSON ID', 'ORDERING PHYSICIAN', 'REPORT', 'PATH TEXT', '# MATCHED BREAST PATH REPORTS', 'OBSERVATION DATE', 'DATE OF BIRTH'], axis = 1, inplace=True)
path.to_csv("PATH_REPORT_Temp.txt", sep='\t', index=False)

SURGERY = pd.concat([SURGERY, path], ignore_index=True )
SURGERY = SURGERY.sort_values(by=['PATIENT_ID', 'START_DATE'], ascending=True)
SURGERY['STYLE_SHAPE'] = SURGERY['EVENT_TYPE']
SURGERY['STYLE_SHAPE'] = SURGERY['STYLE_SHAPE'].str.replace('Surgery', 'circle')
SURGERY['STYLE_SHAPE'] = SURGERY['STYLE_SHAPE'].str.replace('Pathology Report', 'square')
SURGERY.drop_duplicates(subset=['PATIENT_ID', 'EVENT_TYPE', 'SP#'], keep='first', inplace=True)
#Found blank lines
SURGERY['PATIENT_ID'] = SURGERY['PATIENT_ID'].fillna('REMOVE')
SURGERY = SURGERY.loc[~SURGERY['PATIENT_ID'].isin(['REMOVE'])]
SURGERY = SURGERY.loc[~SURGERY['START_DATE'].str.contains('NaT')]

#Low/No clinical data patient filter
SURGERY = SURGERY.loc[~SURGERY['PATIENT_ID'].str.contains('BCM')]
SURGERY = SURGERY.loc[~SURGERY['PATIENT_ID'].str.contains('CW01')]

file_name = str_currentDateTime+"-SURGERY-TIMELINE.txt"
SURGERY.to_csv(file_name, sep='\t', index=False)
del file_name

### Pul RFS to merge with diagBC - breast diagnosis dates so all RFS is from date of diagnosis

RFS = RFS.merge(diagBC, on = 'PATIENT_ID', how='right')
RFS['RECURRENCE_TYPE'] = RFS[['RECURRENCE_TYPE']].fillna('No Noted Recurrence')
RFS['RECUR DATE'] = pd.to_datetime(RFS['RECUR DATE'])
RFS['DATE INITIAL DX'] = pd.to_datetime(RFS['DATE INITIAL DX'], format='%m/%d/%Y')
RFS['RFS_STATUS'] = '1:Progressed'
#Build data for those disease free patients
RFS.loc[RFS['RECURRENCE_TYPE'].str.contains("No Noted Recurrence"),'RFS_STATUS'] = '0:DiseaseFree'
RFS['TODAY'] = pd.to_datetime('today').normalize()
RFS['TEMP'] = RFS['TODAY'] - RFS['DATE INITIAL DX']
RFS['TEMP'] = RFS['TEMP'].astype(str).str.strip(' days')
RFS['TEMP'] = RFS['TEMP'].astype(float)
RFS.loc[RFS['RECUR DAYS AFTER DIAG'].isnull(),'RECUR DAYS AFTER DIAG'] = RFS['TEMP']
RFS['RFS_MONTHS'] = RFS['RECUR DAYS AFTER DIAG'].astype(float).div(30.436875)
#delete and reorder columns - ie remove temps, rebuild today post merge
RFS = RFS[['PATIENT_ID', 'RFS_MONTHS', 'RFS_STATUS', 'DATE INITIAL DX', 'RECURRENCE_TYPE']]
RFS.to_csv("RFSTemp.txt", sep='\t', index=False)
patient = patient.merge(RFS, on = 'PATIENT_ID', how='left')
patient = patient.sort_values(by=['PATIENT_ID', 'DATE INITIAL DX'], ascending=False)
patient.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)

#living OS
patient['TODAY'] = pd.to_datetime('today').normalize()
patient['TEMP'] = patient['TODAY'] - patient['DATE INITIAL DX']
patient['TEMP'] = patient['TEMP'].astype(str).str.strip(' days')
patient['TEMP'] = pd.to_numeric(patient.TEMP, errors='coerce')
patient['TEMP'] = patient[['TEMP']].fillna(0)
patient['TEMP'] = patient['TEMP'].astype(float).div(30.436875)
patient['DATE OF DEATH'] = pd.to_datetime(patient['DATE OF DEATH'], format='%m/%d/%Y')
#deceased OS
patient['TEMP2'] = patient['DATE OF DEATH'] - patient['DATE INITIAL DX']
patient['TEMP2'] = patient['TEMP2'].astype(str).str.strip(' days')
patient['TEMP2'] = pd.to_numeric(patient.TEMP2, errors='coerce')
patient['TEMP2'] = patient['TEMP2'].fillna(0)
patient['TEMP2'] = patient['TEMP2'].astype(float).div(30.436875)
#add a blank float column to accept numbers from temp columns above
patient = patient.assign(OS_MONTHS=pd.Series(dtype='float'))
patient['VITAL STATUS'] = patient['VITAL STATUS'].fillna('remove')
patient.loc[patient['VITAL STATUS'].str.contains("Dead"),'OS_MONTHS'] = patient['TEMP2']
patient.loc[patient['VITAL STATUS'].str.contains("Alive"),'OS_MONTHS'] = patient['TEMP']
#patient['OS_MONTHS'] = patient['OS_MONTHS'].astype(float, errors='ignore')
patient['VITAL STATUS'] = patient['VITAL STATUS'].str.replace('remove', '')

#calculate age at death and AGE_AT_DIAGNOSIS.
patient.rename(columns={'DATE OF BIRTH': 'DATE_OF_BIRTH', 'AGE AT DIAGNOSIS': 'AGE_AT_DIAGNOSIS'}, inplace=True)
patient['DATE_OF_BIRTH'].fillna(0, inplace=True)
patient['DATE_OF_BIRTH'] = pd.to_datetime(patient['DATE_OF_BIRTH'], errors='coerce')
patient['AGE_AT_DIAGNOSIS'] = patient['DATE INITIAL DX'] - patient['DATE_OF_BIRTH']
patient['AGE_AT_DIAGNOSIS'] = patient['AGE_AT_DIAGNOSIS'].astype(str).str.replace(' days', '')
#patient['AGE_AT_DIAGNOSIS'].fillna(0, inplace=True)
#patient['AGE_AT_DIAGNOSIS'] = patient['AGE_AT_DIAGNOSIS'].astype(float).div(365.25)
patient['AGE_AT_DIAGNOSIS'] = pd.to_numeric(patient.AGE_AT_DIAGNOSIS, errors='coerce').div(365.25)

patient['AGE_AT_DEATH'] = patient['DATE OF DEATH'] - patient['DATE_OF_BIRTH']
patient['AGE_AT_DEATH'] = patient['AGE_AT_DEATH'].astype(str).str.replace(' days', '')
#grab status events before month calculation
status2 = patient[['PATIENT_ID', 'DATE_OF_BIRTH', 'DATE OF DEATH', 'VITAL STATUS', 'AGE_AT_DIAGNOSIS']]
patient['DATE_OF_BIRTH'].fillna(0, inplace=True)
patient['AGE_AT_DEATH'] = pd.to_numeric(patient.AGE_AT_DEATH, errors='coerce').div(365.25)
#cleanup columns
patient.rename(columns={'# BODY FLUID AVAIL': 'BODY FLUID AVAIL', '# SOLID TISSUE AVAIL': 'SOLID TISSUE AVAIL', '# POSSIBLE PATH BLOCKS': 'POSSIBLE PATH BLOCKS', '# EXTRACTED BREAST PATH REPORTS': 'EXTRACTED BREAST PATH REPORTS'}, inplace=True)
#patient = patient[['PATIENT_ID', 'HCI PERSON ID', 'GENDER', 'RACE', 'ETHNICITY', 'OS_STATUS', 'OS_MONTHS', 'RFS_MONTHS', 'RFS_STATUS', 'TUMOR REG DX', 'BODY FLUID AVAIL', 'TOTAL BODY FLUID', 'SOLID TISSUE AVAIL', 'TOTAL SOLID TISSUE', 'POSSIBLE PATH BLOCKS', 'EXTRACTED BREAST PATH REPORTS']]


##Build remaing status start dates
StatusDeath = status2.loc[status2['VITAL STATUS'].str.contains('Dead')]
StatusDeath['STOP_DATE'] = ''
StatusDeath['EVENT_TYPE'] = 'STATUS'
StatusDeath['SOURCE'] = ''
#split diagnosis out of deaths for separate lines
DeathDiag = StatusDeath
DeathDiag['SOURCE'] = 'AGE_AT_DIAGNOSIS'
DeathDiag.rename(columns={'AGE_AT_DIAGNOSIS': 'START_DATE'}, inplace=True)
DeathDiag['START_DATE'] = DeathDiag['START_DATE'].astype(float)*365.25
DeathDiag = DeathDiag[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SOURCE']]
DeathDiag = DeathDiag.reset_index()
DeathDiag.to_csv("DeathDiagTemp.txt", sep='\t', index=False)


StatusDeath['SOURCE'] = 'AGE_AT_DEATH'
#StatusDeath.drop(['AGE_AT_DIAGNOSIS'], axis = 1, inplace=True)
StatusDeath = StatusDeath[['PATIENT_ID', 'STOP_DATE', 'EVENT_TYPE', 'SOURCE', 'DATE_OF_BIRTH', 'DATE OF DEATH']]
StatusDeath['AGE_AT_DEATH'] = StatusDeath['DATE OF DEATH'] - StatusDeath['DATE_OF_BIRTH']
StatusDeath['AGE_AT_DEATH'] = StatusDeath['AGE_AT_DEATH'].astype(str).str.replace(' days', '')
#StatusDeath[['AGE_AT_DEATH']] = StatusDeath[['AGE_AT_DEATH']].fillna('REMOVE')
StatusDeath = StatusDeath.loc[~StatusDeath['AGE_AT_DEATH'].str.contains('NaT')]
StatusDeath.rename(columns={'AGE_AT_DEATH': 'START_DATE'}, inplace=True)

StatusDeath = StatusDeath[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SOURCE']]
StatusDeath = StatusDeath.reset_index()
StatusDeath.to_csv("StatusDeathTemp.txt", sep='\t', index=False)
StatusDeath = pd.read_csv('StatusDeathTemp.txt', sep='\t')

StatusAlive = status2.loc[status2['VITAL STATUS'].str.contains('Alive')]
StatusAlive['STOP_DATE'] = ''
StatusAlive['EVENT_TYPE'] = 'STATUS'
StatusAlive['SOURCE'] = 'AGE_AT_DIAGNOSIS'
StatusAlive.rename(columns={'AGE_AT_DIAGNOSIS': 'START_DATE'}, inplace=True)
StatusAlive['START_DATE'] = StatusAlive['START_DATE'].astype(float)*365.25
StatusAlive = StatusAlive[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SOURCE']]
StatusAlive = StatusAlive.reset_index()
StatusAlive.to_csv("StatusAliveTemp.txt", sep='\t', index=False)

RecurDF = RecurDF.reset_index()
RecurDF = pd.concat([RecurDF, DeathDiag])
RecurDF = pd.concat([RecurDF, StatusDeath]) 
RecurDF = pd.concat([RecurDF, StatusAlive]) 
RecurDF = RecurDF.sort_values(by=['PATIENT_ID', 'START_DATE'], ascending=True)
RecurDF.drop(['index'], axis = 1, inplace=True)
RecurDF.rename(columns={'SOURCE': 'STATUS_EVENT'}, inplace=True)

#default recurrence to purple circles and change deaths and diagnosis
RecurDF['STYLE_COLOR'] = '#AF7AC5'
RecurDF['STYLE_SHAPE'] = 'circle'
#black square death
RecurDF.loc[RecurDF['STATUS_EVENT'].str.contains('AGE_AT_DEATH'), 'STYLE_COLOR'] = '#040404'
RecurDF.loc[RecurDF['STATUS_EVENT'].str.contains('AGE_AT_DEATH'), 'STYLE_SHAPE'] = 'square'
#red
RecurDF.loc[RecurDF['STATUS_EVENT'].str.contains('AGE_AT_DIAGNOSIS'), 'STYLE_COLOR'] = '#E74C3C'
RecurDF.loc[RecurDF['STATUS_EVENT'].str.contains('AGE_AT_DIAGNOSIS'), 'STYLE_SHAPE'] = 'diamond'
RecurDF['START_DATE'] = RecurDF['START_DATE'].astype(int)

#Low/No clinical data patient filter
RecurDF = RecurDF.loc[~RecurDF['PATIENT_ID'].str.contains('BCM')]
RecurDF = RecurDF.loc[~RecurDF['PATIENT_ID'].str.contains('CW01')]

file_name = str_currentDateTime+"-STATUS-TIMELINE.txt"
RecurDF.to_csv(file_name, sep='\t', index=False)
del file_name

# Read and preprocess the first file
df3a = pd.read_csv('KTVMed1.txt', sep='|', encoding='utf-8')
df3a = df3a.drop([0])
df3a.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df3a = df3a.loc[df3a['PATIENT_ID'].isin(patient_list)]
df3a.rename(columns=str.upper, inplace=True)

# Read and preprocess the second file
df3b = pd.read_csv('KTVMed2.txt', sep='|', encoding='utf-8')
df3b = df3b.drop([0])
df3b.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df3b = df3b.loc[df3b['PATIENT_ID'].isin(patient_list)]
df3b.rename(columns=str.upper, inplace=True)

# Concatenate the dataframes
df3 = pd.concat([df3a, df3b], ignore_index=True)

# Filter the dataframe
df3 = df3.loc[df3['MEDICATION CLASS'].isin(["ESTROGENS", "ENDOCRINE AND METABOLIC AGENTS - MISC.", "CHEMO", "IMMUNO", "HORMONE", "ANTINEOPLASTICS AND ADJUNCTIVE THERAPIES"])]

# Merge with date of birth dataframe
df3 = df3.merge(dob, on='PATIENT_ID', how='left')

# Insert new columns
df3['START_DATE'] = ''
df3['STOP_DATE'] = ''
df3['EVENT_TYPE'] = 'TREATMENT'
df3['TREATMENT_TYPE'] = 'MEDICATION'

# Drop unnecessary columns
df3.drop(['HCI PERSON ID', 'AGE AT START', 'AGE AT END'], axis=1, inplace=True)

# Process date columns
df3['START DATE'] = df3['START DATE'].astype(str).str.split(' ', expand=True)[0]
df3['END DATE'] = df3['END DATE'].astype(str).str.split(' ', expand=True)[0]
df3['START DATE'] = pd.to_datetime(df3['START DATE'])
df3['END DATE'] = pd.to_datetime(df3['END DATE'])
df3['START_DATE'] = (df3['START DATE'] - df3['DATE OF BIRTH']).astype(str).str.replace(' days', '')
df3['STOP_DATE'] = (df3['END DATE'] - df3['DATE OF BIRTH']).astype(str).str.replace(' days', '').str.replace('NaT', '')

# Drop original date columns
df3.drop(['END DATE', 'START DATE', 'DATE OF BIRTH', 'CURRENT?', 'MEDICATION CLASS', 'DOSE MIN', 'DOSE MAX', 'DOSE UNIT'], axis=1, inplace=True)

# Rename columns
df3.rename(columns={'MEDICATION': 'SUBTYPE_AGENT', 'GPI CODE': 'CPT CODE'}, inplace=True)

# Drop duplicates
df3.drop_duplicates(subset=['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'TREATMENT_TYPE', 'SUBTYPE_AGENT'], inplace=True)

## Therapy total
df3_tx = pd.concat([chemo, rad, immuno, hormone, SYSTEMIC], ignore_index=True)

df3_tx.rename(columns={'SUBTYPE_AGENT': 'AGENT'}, inplace=True)
df3_tx.drop(['DATE OF BIRTH'], axis=1, inplace=True)
df3_tx.drop(['RAD NUM TX'], axis=1, inplace=True)
df3_tx.rename(columns={'CPT CODE': 'CPT/GPI CODE'}, inplace=True)

df3_tx['START_DATE'] = df3_tx['START_DATE'].str.replace('NaT', '0')

# Found excess duplicates where one had stop_date and the duplicate entry did not.
df3_tx = df3_tx.sort_values(by=['PATIENT_ID', 'START_DATE', 'STOP_DATE'], ascending=False)
df3_tx.drop_duplicates(subset=['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'AGENT'], keep='first', inplace=True)
df3_tx['START_DATE'] = df3_tx['START_DATE'].astype(str)
df3_tx['START_DATE'] = df3_tx['START_DATE'].str.replace('nan', '0')

df3_tx['START_DATE'] = df3_tx['START_DATE'].str.extract(r'(\d+)').astype(int)
df3_tx = df3_tx[(df3_tx['START_DATE'] > 0)]
# Found duplicate meds with the same start date but absent stop dates ending in duplicates - sort to eliminate those dups without stop dates
df3_tx = df3_tx.sort_values(by=['PATIENT_ID', 'AGENT', 'START_DATE', 'STOP_DATE'], ascending=False)
df3_tx.drop_duplicates(subset=['PATIENT_ID', 'START_DATE', 'AGENT'], keep='first', inplace=True)

#Low/No clinical data patient filter
df3_tx = df3_tx.loc[~df3_tx['PATIENT_ID'].str.contains('BCM')]
df3_tx = df3_tx.loc[~df3_tx['PATIENT_ID'].str.contains('CW01')]

file_name = str_currentDateTime + "-TREATMENTS-TIMELINE.txt"
df3_tx.to_csv(file_name, sep='\t', index=False)
del file_name

print('\n\n\***************** LABS **********************\n\n')
df4 = pd.read_csv('KTVArupPath.txt', sep='|', encoding= 'utf-8')
df4 = df4.drop([0])
df4.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df4 = df4.loc[df4['PATIENT_ID'].isin(patient_list)]
df4.rename(columns = str.upper , inplace = True)
df4 = df4.merge(dob, on = 'PATIENT_ID', how='left')
df4.insert(1, column='START_DATE', value='')
df4.insert(2, column='STOP_DATE', value='')
df4.insert(3, column='EVENT_TYPE', value='LAB_TEST')
df4.insert(4, column='TEST', value='')
df4.drop(['HCI PERSON ID'], axis = 1, inplace=True)
df4.drop(['PROCEDURE DATE'], axis = 1, inplace=True)
df4.drop(['SP#'], axis = 1, inplace=True)
df4['REPORT DATE'] = pd.to_datetime(df4['REPORT DATE'])
df4['START_DATE'] = df4['REPORT DATE'] - df4['DATE OF BIRTH']
df4['START_DATE'] = df4['START_DATE'].astype(str)
df4['START_DATE'] = df4['START_DATE'].str.replace(' days', '')
df4.drop(['REPORT DATE'], axis = 1, inplace=True)
df4.drop(['DATE OF BIRTH'], axis = 1, inplace=True)
df4.drop_duplicates(subset=['PATIENT_ID', 'START_DATE'], keep='first', inplace=True)


#split into different tests to make flat file with one row for each result for display purposes
ER = df4[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'TEST', 'ER']]
ER['TEST'] = "ER"
ER.rename(columns={'ER': 'RESULT'}, inplace=True)

PR = df4[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'TEST', 'PR']]
PR['TEST'] = "PR"
PR.rename(columns={'PR': 'RESULT'}, inplace=True)

HER2 = df4[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'TEST', 'HER2']]
HER2['TEST'] = "HER2_IHC"
HER2.rename(columns={'HER2': 'RESULT'}, inplace=True)

HER2fish = df4[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'TEST', 'HER2 FISH', 'HER2 SCORE']]
HER2fish[['HER2 FISH']] = HER2fish[['HER2 FISH']].fillna('BLANK')
HER2fish['HER2 FISH'] = HER2fish['HER2 FISH'] + ': ' + HER2fish['HER2 SCORE']
HER2fish['TEST'] = "HER2_FISH"
HER2fish.rename(columns={'HER2 FISH': 'RESULT'}, inplace=True)
HER2fish.drop(['HER2 SCORE'], axis = 1, inplace=True)

df4b = pd.concat([ER, PR, HER2, HER2fish], ignore_index=True)
df4b = df4b.sort_values(by=['PATIENT_ID', 'START_DATE'])
df4b[['RESULT']] = df4b[['RESULT']].fillna('REMOVE')
df4b = df4b.loc[~df4b['RESULT'].isin(['REMOVE'])]


print('\n\n\***************** LIQUID SPECIMENS **********************\n\n')
df5 = pd.read_csv('KTVBodyFluid.txt', sep='|', encoding= 'utf-8')
df5 = df5.drop([0])
df5.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df5['PATIENT_ID'] = df5['PATIENT_ID'].str.strip()
df5 = df5.loc[df5['PATIENT_ID'].isin(patient_list)]
df5.rename(columns = str.upper , inplace = True)
df5.rename(columns={'DEPLETED?': 'DEPLETED'}, inplace=True)
df5 = df5.merge(dob, on = 'PATIENT_ID', how='left')
df5.insert(1, column='START_DATE', value='')
df5.insert(2, column='STOP_DATE', value='')
df5.insert(3, column='EVENT_TYPE', value='SPECIMEN')
df5.insert(4, column='SPECIMEN_TYPE', value='Liquid Samples')
#df5.insert(5, column='TISSUE TYPE', value='')
#df5.to_csv("df5.txt", sep='\t', index=False)
df5.drop(['HCI PERSON ID'], axis = 1, inplace=True)
df5['AVAIL AMOUNT'] = df5['AVAIL AMOUNT'].astype(float)
df5['AVAIL AMOUNT'] = df5['AVAIL AMOUNT'].round(decimals=2)
df5['AVAIL AMOUNT'] = df5['AVAIL AMOUNT'].astype(str)
df5['AVAIL AMOUNT'] = df5['AVAIL AMOUNT'].str.replace('0.0', '')
df5['AVAIL AMOUNT'] = df5['AVAIL AMOUNT'] + df5['AMOUNT UNIT'].astype(str)
df5['CONCENTRATION'] = df5['CONCENTRATION'].astype(float)
df5['CONCENTRATION'] = df5['CONCENTRATION'].round(decimals=2)
df5['CONCENTRATION'] = df5['CONCENTRATION'].astype(str) + df5['CONCENTRATION UNIT'].astype(str)
df5['CONCENTRATION'] = df5['CONCENTRATION'].str.replace('nannan', '')
df5['CONCENTRATION'] = df5['CONCENTRATION'].str.replace('nan', '')
df5.drop(['AMOUNT UNIT', 'CONCENTRATION UNIT'], axis = 1, inplace=True)
#df5['DATE OF BIRTH'] = pd.to_datetime(df5['DATE OF BIRTH'])
#df5['DATE OF BIRTH'] = pd.to_datetime(df5['DATE OF BIRTH'], format='%Y/%m/%d')
df5['COLLECT DATE'] = pd.to_datetime(df5['COLLECT DATE'])
df5['COLLECT DATE'] = pd.to_datetime(df5['COLLECT DATE'], format='%Y/%m/%d')
df5['START_DATE'] = df5['COLLECT DATE'] - df5['DATE OF BIRTH']
df5['START_DATE'] = df5['START_DATE'].astype(str)
df5['START_DATE'] = df5['START_DATE'].str.replace(' days', '')
df5.drop(['COLLECT DATE', 'DATE OF BIRTH', 'CONCENTRATION SOURCE', 'CELL COUNT', 'CELL COUNT UNIT'], axis = 1, inplace=True)
df5 = df5 [['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SPECIMEN_TYPE', 'TISSUE TYPE', 'CC#', 'SAMPLE TYPE', 'PREP METHOD', 'COLLECTION PREP TYPE', 'STORAGE PREP TYPE', 'AVAIL AMOUNT', 'CONCENTRATION', 'COLLECTION PROTOCOL', 'PROJECT(S)', 'DEPLETED']]
df5['SOURCE'] = 'BodyFluid'
df5.to_csv("LiquidSpecimensTemp.txt", sep='\t', index=False)


print('\n\n\***************** SOLID TISSUE SPECIMENS **********************\n\n')
df6 = pd.read_csv('KTVSolidTissue.txt', sep='|', encoding= 'utf-8')
df6 = df6.drop([0])
df6.rename(columns={'Shadow ID': 'PATIENT_ID'}, inplace=True)
df6['PATIENT_ID'] = df6['PATIENT_ID'].str.strip()
df6 = df6.loc[df6['PATIENT_ID'].isin(patient_list)]
df6.rename(columns = str.upper , inplace = True)
df6.rename(columns={'DEPLETED?': 'DEPLETED' }, inplace=True)
df6 = df6.merge(dob, on = 'PATIENT_ID', how='left')

df6.insert(1, column='START_DATE', value='')
df6.insert(2, column='STOP_DATE', value='')
df6.insert(3, column='EVENT_TYPE', value='SPECIMEN')
df6.insert(4, column='SPECIMEN_TYPE', value='Solid Tissue Samples')
#df5.insert(5, column='SAMPLE_TYPE', value='')
df6.drop(['HCI PERSON ID'], axis = 1, inplace=True)
df6['AVAIL AMOUNT'] = df6['AVAIL AMOUNT'].astype(float)
df6['AVAIL AMOUNT'] = df6['AVAIL AMOUNT'].round(decimals=3)
df6['AVAIL AMOUNT'] = df6['AVAIL AMOUNT'].astype(str) + df6['AMOUNT UNIT'].astype(str)
df6['AVAIL AMOUNT'] = df6['AVAIL AMOUNT'].str.replace('nannan', '')
df6['AVAIL AMOUNT'] = df6['AVAIL AMOUNT'].str.replace('nan', '')
df6.drop(['AMOUNT UNIT'], axis = 1, inplace=True)
df6['COLLECT DATE'] = pd.to_datetime(df6['COLLECT DATE']).dt.date
df6['COLLECT DATE'] = pd.to_datetime(df6['COLLECT DATE'], format='%Y/%m/%d')
df6['START_DATE'] = df6['COLLECT DATE'] - df6['DATE OF BIRTH']
df6['START_DATE'] = df6['START_DATE'].astype(str)
df6['START_DATE'] = df6['START_DATE'].str.replace(' days', '')
df6.drop(['COLLECT DATE', 'DATE OF BIRTH'], axis = 1, inplace=True)
df6 = df6 [['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SPECIMEN_TYPE', 'TISSUE TYPE', 'CC#', 'MATCHED ARUP PATHOLOGY', 'SAMPLE TYPE', 'PREP METHOD', 'COLLECTION PREP TYPE', 'STORAGE PREP TYPE', 'AVAIL AMOUNT', 'COLLECTION PROTOCOL', 'PROJECT(S)', 'DEPLETED']]
df6.drop_duplicates(subset=['PATIENT_ID', 'START_DATE', 'CC#'], keep='first', inplace=True)
df6['SOURCE'] = 'SolidTissue'
df6.to_csv("SolidSpecimenTemp.txt", sep='\t', index=False)

##merge in liquid samples and labs
df5 = pd.concat([df5, df6], ignore_index=True)
LifeSample = df5.loc[df5['SPECIMEN_TYPE'].isin(['Solid Tissue Samples', 'Liquid Samples'])]
LifeSample['CC#'] = LifeSample['CC#'].astype(str)
LifeSample['CC#'] = LifeSample['CC#'].str.split('.',expand=True)[0]
#remove alpabetic characters
LifeSample['CCtemp2'] = LifeSample['CC#']
LifeSample['CC#'] = LifeSample['CC#'].astype(str).str.replace(r'\D+', '')
LifeSample.drop_duplicates(subset=['START_DATE', 'CC#'], keep='first', inplace=True)
LifeSample.rename(columns={'CC#': 'CC'}, inplace=True)
LifeSample['CC'] = LifeSample['CC'].astype(str)
LifeSample['CC'] = LifeSample['CC'].str.strip()


print('\n\n\n\*********************************** dytpes here:\n\n')
print(pd.api.types.infer_dtype(LifeSample['CC']))
##CC trim
#LifeSample['CC'] = LifeSample['CC'].str[:-1]
sample['CC'] = sample['CC'].astype(str)
print(pd.api.types.infer_dtype(sample['CC']))
#use to match sans last character but have now taken that out. 
#sample['CC'] = sample['CC'].str[:-1]

#merge wtih Sample data
LifeSample.drop(['SAMPLE TYPE', 'TISSUE TYPE', 'PREP METHOD', 'STORAGE PREP TYPE', 'AVAIL AMOUNT', 'COLLECTION PROTOCOL', 'PROJECT(S)', 'DEPLETED'], axis = 1, inplace=True)
LifeSample2 = pd.merge(LifeSample, sample,  on =['PATIENT_ID', 'CC'], how='outer')

#Consider taking new sample information 

#LifeSample2.drop(['TISSUE TYPE', 'PREP METHOD', 'STORAGE PREP TYPE', 'AVAIL AMOUNT', 'COLLECTION PROTOCOL', 'PROJECT(S)', 'DEPLETED'], axis = 1, inplace=True)
###LifeSample2 = LifeSample2[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SAMPLE_ID', 'SAMPLE_TYPE', 'SPECIMEN_TYPE', 'PRIMARY_OR_RECURRENCE', 'CCtemp', 'COLLECTION_SITE', 'COLLECTION PREP TYPE', 'MATCHED ARUP PATHOLOGY', 'WES_ID', 'RNASEQ', 'ER_PDX', 'PR_PDX', 'HER2_PDX', 'PDX_GROW_NOGROW', 'SOURCE']]
LifeSample2['CCtemp'] = LifeSample2['CCtemp'].fillna(LifeSample2['CCtemp2'])
LifeSample2.to_csv('CCmergeTemp.txt', sep='\t', index=False)
LifeSample2.rename(columns={'CCtemp': 'CC#'}, inplace=True)
#Duplicates 
LifeSample2 = LifeSample2.sort_values(by=['PATIENT_ID', 'SAMPLE_ID', 'SPECIMEN_TYPE'], ascending=False)
LifeSample2.drop_duplicates(subset=['PATIENT_ID', 'START_DATE', 'SAMPLE_ID', 'CC#'], keep='first', inplace=True)
#take dates from provided sample data where CC# not available. 
LifeSample2['START_DATE'] = LifeSample2['START_DATE'].fillna(LifeSample2['DATE_OF_COLLECTION_(ONLY_SAMPLES_W-O_SHADOWID)'])

#split specimens from sequencing
LifeSample2['SAMPLE_ID'].fillna('UNUSED SPECIMEN', inplace=True)
#Build expanded sample file
del sample
sample = LifeSample2.loc[~LifeSample2['SAMPLE_ID'].isin(['UNUSED SPECIMEN'])]
SampleCDEmeta = patient.dtypes
sample = pd.concat([sample, SampleCDEmeta], ignore_index=True)
cols = sample.columns.tolist()
sample.loc[-1] = cols
sample = sample.reindex(np.roll(sample.index, shift=1))
sample.loc[-1] = ''  # adding a row
sample = sample.reindex(np.roll(sample.index, shift=1))
#sample.to_csv("SAMPLE_AUDIT.txt", sep='\t', index=False)
#Splitting specimens may look odd on timeline. Take subset for sequencing but keep all for specimens. 
#SPECIMENS = LifeSample2.loc[LifeSample2['SAMPLE_ID'].isin(['SPECIMEN'])]
SPECIMENS = LifeSample2
SPECIMENS = SPECIMENS[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SAMPLE_ID', 'SAMPLE_TYPE', 'SPECIMEN_TYPE', 'PRIMARY_OR_RECURRENCE', 'CC#', 'COLLECTION_SITE', 'COLLECTION PREP TYPE', 'MATCHED ARUP PATHOLOGY']]

SPECIMENS.drop_duplicates(subset=['START_DATE', 'SAMPLE_ID'], keep='first', inplace=True)
SPECIMENS.rename(columns={'CC#': 'SPECIMEN_ID'}, inplace=True)
SPECIMENS['EVENT_TYPE'] = 'SPECIMEN'
SPECIMENS.rename(columns={'MATCHED ARUP PATHOLOGY': 'MATCHED PATHOLOGY ACCESSION'}, inplace=True)
SPECIMENS['START_DATE'] = SPECIMENS['START_DATE'].fillna(0).replace('', 0).astype(int)
SPECIMENS['START_DATE'] = SPECIMENS['START_DATE'].astype(int)
SPECIMENS = SPECIMENS.sort_values(by=['PATIENT_ID', 'START_DATE'])
#Drop blanks as some came through parser using xlsx starter
SPECIMENS[['PATIENT_ID']] = SPECIMENS[['PATIENT_ID']].fillna('REMOVE')
SPECIMENS = SPECIMENS.loc[~SPECIMENS['PATIENT_ID'].isin(['REMOVE'])]
#tumor on top
SPECIMENS = SPECIMENS.iloc[(~SPECIMENS['SAMPLE_ID'].str.contains('tumor')).argsort()]
#patient specimen counts
SpecimenCount = SPECIMENS['PATIENT_ID'].value_counts(dropna=False).to_frame()

#Low/No clinical data patient filter
SPECIMENS = SPECIMENS.loc[~SPECIMENS['PATIENT_ID'].str.contains('BCM')]
SPECIMENS = SPECIMENS.loc[~SPECIMENS['PATIENT_ID'].str.contains('CW01')]

file_name = str_currentDateTime+"-SPECIMENS-TIMELINE.txt"
SPECIMENS.to_csv(file_name, sep='\t', index=False)
del file_name


#Build model timeline event
LifeSample2 = LifeSample2.loc[~LifeSample2['SAMPLE_ID'].isin(['UNUSED SPECIMEN'])]
LifeSample2.drop(['SOURCE'], axis = 1, inplace=True)
LifeSample2['EVENT_TYPE'] = 'Sequencing'

#format timeline display
df4b['STYLE_COLOR'] = df4b['TEST']
df4b['STYLE_SHAPE'] = df4b['TEST']
LifeSample2['STYLE_COLOR'] = LifeSample2['SAMPLE_TYPE']
LifeSample2['STYLE_SHAPE'] = LifeSample2['PDX_GROW_NOGROW']

df4b['STYLE_COLOR'] = df4b['STYLE_COLOR'].str.replace('PR', '#AF7AC5')
df4b['STYLE_COLOR'] = df4b['STYLE_COLOR'].str.replace('HER2_IHC', '#2ECC71')
df4b['STYLE_COLOR'] = df4b['STYLE_COLOR'].str.replace('HER2_FISH', '#2ECC71')
df4b['STYLE_COLOR'] = df4b['STYLE_COLOR'].str.replace('ER', '#E74C3C')
df4b['STYLE_COLOR'].fillna('#3498DB', inplace=True)

df4b['STYLE_SHAPE'] = df4b['STYLE_SHAPE'].str.replace('PR', 'diamond')
df4b['STYLE_SHAPE'] = df4b['STYLE_SHAPE'].str.replace('HER2_IHC', 'square')
df4b['STYLE_SHAPE'] = df4b['STYLE_SHAPE'].str.replace('HER2_FISH', 'square')
df4b['STYLE_SHAPE'] = df4b['STYLE_SHAPE'].str.replace('ER', 'triangle')

#Grab model information for lab_test
Model = LifeSample2[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SAMPLE_ID', 'SAMPLE_TYPE' , 'PDX_GROW_NOGROW']]
#black
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('patient_tumor', '#040404')
#red
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('PDXEv', '#E74C3C')
#purple
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('pdxox', '#AF7AC5')
#bluegreen
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('pdxo', '#29b9b4')
#green
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('pdx', '#2ECC71')
#yellow
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('pdox', '#b9b929')
#orange
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('pdo', '#B99729')

#grey
LifeSample2['STYLE_COLOR'] = LifeSample2['STYLE_COLOR'].str.replace('cells', '#6b6669')
LifeSample2['STYLE_COLOR'].fillna('#3498DB', inplace=True)

LifeSample2['STYLE_SHAPE'] = LifeSample2['STYLE_SHAPE'].str.replace('NoGrow', 'square')
LifeSample2['STYLE_SHAPE'] = LifeSample2['STYLE_SHAPE'].str.replace('Grow', 'circle')
LifeSample2['STYLE_SHAPE'].fillna('triangle', inplace=True)
#drop hard coded columns that come over with sample data and duplicate. 
LifeSample2.drop(['CC', 'CCtemp2', 'STYLE_SHAPE', 'STYLE_COLOR'], axis = 1, inplace=True)
LifeSample2.rename(columns={'CC#': 'SPECIMEN_ID'}, inplace=True)
#tumor on top
LifeSample2 = LifeSample2.iloc[(~LifeSample2['SAMPLE_ID'].str.contains('tumor')).argsort()]
LifeSample2['START_DATE'] = LifeSample2['START_DATE'].fillna(0).replace('', 0).astype(int)
LifeSample2['START_DATE'] = LifeSample2['START_DATE'].astype(int)

#Low/No clinical data patient filter
LifeSample2 = LifeSample2.loc[~LifeSample2['PATIENT_ID'].str.contains('BCM')]
LifeSample2 = LifeSample2.loc[~LifeSample2['PATIENT_ID'].str.contains('CW01')]

file_name = str_currentDateTime+"-SEQUENCING-TIMELINE.txt"
LifeSample2.to_csv(file_name, sep='\t', index=False)
del file_name



Model['EVENT_TYPE'] = 'LAB_TEST'
Model['STYLE_COLOR'] = Model['SAMPLE_TYPE']
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.upper()
Model.rename(columns={'SAMPLE_TYPE': 'TEST'}, inplace=True)
#Add space to move models to top
Model['TEST'] = ' ' + Model['TEST']
Model['TEST'] = Model['TEST'].str.replace('Patient_tumor', 'Tumor')
Model.insert(4, column='SUBTYPE', value='MODEL')


#black
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('PATIENT_TUMOR', '#040404')
#red
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('PDXEV', '#E74C3C')
#purple
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('PDXOX', '#AF7AC5')
#bluegreen
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('PDXO', '#29b9b4')
#green
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('PDX', '#2ECC71')
#yellow
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('PDOX', '#b9b929')
#orange
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('PDO', '#B99729')

#grey
Model['STYLE_COLOR'] = Model['STYLE_COLOR'].str.replace('CELLS', '#6b6669')
Model['STYLE_COLOR'].fillna('#3498DB', inplace=True)
#Sort tumor first
Model = Model.sort_values(by=['TEST'], ascending=[False])
Model = Model.iloc[(~Model['TEST'].str.contains('Tumor')).argsort()]

#Model.to_csv("Model.txt", sep='\t', index=False)

df4b.insert(4, column='SUBTYPE', value='HORMONE RECEPTOR')
df4b.insert(7, column='SAMPLE_ID', value='')
#invert append to have models first
Model = pd.concat([Model, df4b], ignore_index=True)
Model = Model[['PATIENT_ID', 'START_DATE', 'STOP_DATE', 'EVENT_TYPE', 'SUBTYPE', 'TEST', 'RESULT', 'SAMPLE_ID', 'PDX_GROW_NOGROW', 'STYLE_COLOR', 'STYLE_SHAPE']]

Model['TEST'] = Model['TEST'].str.strip()
sort_order = ['Tumor', 'PDX', 'PDXEv', 'PDxO', 'PDXoX', 'PDO', 'PDoX', 'Cells', 'HER2_IHC', 'HER2_FISH', 'ER', 'PR']
Model['TEST'] = pd.Categorical(Model['TEST'], categories=sort_order, ordered=True)
Model = Model.sort_values(by=['TEST'])

#Model = Model.iloc[(~Model['TEST'].str.contains('Tumor')).argsort()]
Model['STYLE_SHAPE'].fillna('circle', inplace=True)
Model['START_DATE'] = Model['START_DATE'].fillna(0).replace('', 0).astype(int)

Model['START_DATE'] = Model['START_DATE'].astype(int)

#Low/No clinical data patient filter
Model = Model.loc[~Model['PATIENT_ID'].str.contains('BCM')]
Model = Model.loc[~Model['PATIENT_ID'].str.contains('CW01')]

file_name = str_currentDateTime+"-LAB_RESULTS-TIMELINE.txt"
Model.to_csv(file_name, sep='\t', index=False)
del file_name

DIAGNOSTICS = Model.loc[Model['SUBTYPE'].isin(['MODEL'])]
DIAGNOSTICS['EVENT_TYPE'] = 'DIAGNOSTICS'

#Low/No clinical data patient filter
DIAGNOSTICS = DIAGNOSTICS.loc[~DIAGNOSTICS['PATIENT_ID'].str.contains('BCM')]
DIAGNOSTICS = DIAGNOSTICS.loc[~DIAGNOSTICS['PATIENT_ID'].str.contains('CW01')]

file_name = str_currentDateTime+"-DIAGNOSTICS_MODEL-TIMELINE.txt"
DIAGNOSTICS.to_csv(file_name, sep='\t', index=False)
del file_name


######################  FINAL HEADER CORRECTION #################
#add cBio specific header meta data
patient.drop_duplicates(subset=['PATIENT_ID'], keep='first', inplace=True)

#Found blank lines
patient['PATIENT_ID'] = patient['PATIENT_ID'].fillna('REMOVE')
patient = patient.loc[~patient['PATIENT_ID'].isin(['REMOVE'])]
patient = patient[['PATIENT_ID', 'HCI PERSON ID', 'GENDER', 'RACE', 'ETHNICITY', 'OS_STATUS', 'VITAL STATUS', 'AGE_AT_DIAGNOSIS', 'AGE_AT_DEATH', 'OS_MONTHS', 'RFS_MONTHS', 'RFS_STATUS', 'RECURRENCE_TYPE', 'TUMOR REG DX', 'BODY FLUID AVAIL', 'TOTAL BODY FLUID', 'SOLID TISSUE AVAIL', 'TOTAL SOLID TISSUE', 'POSSIBLE PATH BLOCKS', 'EXTRACTED BREAST PATH REPORTS']]
# Determine data types and replace with STRING/NUMBER
dtype_mapping = patient.dtypes.replace({
    'object': 'STRING',
    'float64': 'NUMBER',
    'int64': 'NUMBER'
})

# Convert dtype_mapping to a single-row DataFrame
dtype_df = pd.DataFrame([dtype_mapping], columns=patient.columns)
# Add '#' prefix to the first column in dtype_df (Row 2)
dtype_df.iloc[:, 0] = "#" + dtype_df.iloc[:, 0].astype(str)

# Create other required rows
header_dup_1 = pd.DataFrame([patient.columns], columns=patient.columns)  # Duplicate header (Row 1)
header_dup_2 = pd.DataFrame([patient.columns], columns=patient.columns)  # Duplicate header (Row 4)
ones_row = pd.DataFrame([["1"] * len(patient.columns)], columns=patient.columns)  # Row of "1" (Row 3)
# Add '#' prefix to the first column in header_dup_1, header_dup_2, and ones_row
header_dup_1.iloc[:, 0] = "#" + header_dup_1.iloc[:, 0].astype(str)
ones_row.iloc[:, 0] = "#" + ones_row.iloc[:, 0].astype(str)

# Combine all components
patient = pd.concat([header_dup_1, dtype_df, ones_row, header_dup_2, patient], ignore_index=True)
patient.rename(columns={'PATIENT_ID': '#PATIENT_ID'}, inplace=True)

file_name = str_currentDateTime+"-PATIENT_ATTRIBUTES.txt"
patient.to_csv(file_name, sep='\t', index=False)
del file_name

# Clean up temp files
dir_name = os.getcwd()
test = os.listdir(dir_name)

for item in test:
    if item.endswith("Temp.txt"):
        os.remove(os.path.join(dir_name, item))


# Set up argument parser
parser = argparse.ArgumentParser(description="Generate meta files for text files.")
parser.add_argument('-s', '--study', type=str, default="XXXXX", help="Cancer study identifier")
parser.add_argument('-c', '--cancer', type=str, default="breast", help="Type of cancer")
args = parser.parse_args()

# Define the meta prefix
meta_prefix = "meta"

# Function to generate meta file for each data file
def generate_meta_files(cancer_study_identifier, str_currentDateTime):
    # List all files in the current directory
    files = os.listdir('.')
    
    # Iterate over each file
    for file_name in files:
        # Check if the file is a text file
        if file_name.startswith(str_currentDateTime):
            # Generate the meta file name
            meta_file_name = meta_prefix + file_name.replace(str_currentDateTime, "")
            
            # Extract the last _ delimited string prior to .txt
            datatype = file_name.rsplit('-', 1)[-1].replace('.txt', '')
            
            # Write meta information to the meta file
            with open(meta_file_name, 'w') as meta_file:
                meta_file.write(f"cancer_study_identifier: {cancer_study_identifier}\n")
                meta_file.write("genetic_alteration_type: CLINICAL\n")
                meta_file.write(f"datatype: {datatype}\n")
                meta_file.write(f"data_filename: {file_name}\n")
# Function to generate the meta-study file
def generate_meta_study_file(cancer_study_identifier, type_of_cancer, str_currentDateTime):
    meta_study_content = f"""type_of_cancer: {type_of_cancer}
cancer_study_identifier: {cancer_study_identifier}
name: HCI - {cancer_study_identifier} {str_currentDateTime}: Patient, Germline, and PDX/PDO/PDOX samples - Welm/Varley
description: Whole exome of tumor and tumor-normal samples with SNVs, INDELs, CNVs, treatment responses, RNAseq
add_global_case_list: true
reference_genome: hg38"""
    
    with open("meta-study.txt", 'w') as meta_study_file:
        meta_study_file.write(meta_study_content)

def generate_meta_generic_assay(cancer_study_identifier, type_of_cancer, str_currentDateTime):
    meta_generic_assay_content = f"""cancer_study_identifier: {cancer_study_identifier}
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: HCI_ID_LEGEND
datatype: LIMIT-VALUE
stable_id: PATIENT_SAMPLES
profile_name: PATIENT SAMPLES LEGEND
profile_description: HCI_ID grouped samples for OncoPrint Legend. 
data_filename: {str_currentDateTime}-PATIENT_ONCOPRINT_LEGEND.txt
show_profile_in_analysis_tab: true
pivot_threshold_value: 0
value_sort_order: ASC
generic_entity_meta_properties: NAME,DESCRIPTION"""

    meta_file_name = f"meta-{str_currentDateTime}-PATIENT_ONCOPRINT_LEGEND.txt"
    with open(meta_file_name, 'w') as meta_file:
        meta_file.write(meta_generic_assay_content)

# Generate meta files
generate_meta_files(args.study, str_currentDateTime)

# Generate the meta-study file
generate_meta_study_file(args.study, args.cancer, str_currentDateTime)

print("Meta files generated successfully.")

# Remove the specific file if it exists
specific_file = "meta-MolecularSampleSheet.txt"
specific_file_path = os.path.join(dir_name, specific_file)
if os.path.exists(specific_file_path):
    os.remove(specific_file_path)
else:
    print(f"{specific_file} does not exist in the directory.")
    
#Return to case lists now that all files are present and arguments are in place
#other case lists should be made for treatment data.
#write case_all files:
data = f"""cancer_study_identifier: {args.study}
stable_id: {args.study}_
case_list_name: All samples ({cases_all})
case_list_description: This is this case list that contains all samples that are profiled.
case_list_ids:\t{cases_all_delim}"""

# Directory where the file will be saved
directory = "case_lists"

# Create the directory if it doesn't exist
if not os.path.exists(directory):
    os.makedirs(directory)

# File path
file_path = os.path.join(directory, "cases_all.txt")

# Write data to the file
with open(file_path, 'w') as file:
    file.write(data)
del data

#write cases_sequenced files:
data = f"""cancer_study_identifier: {args.study}
stable_id: {args.study}_sequenced
case_list_name: Samples with mutation data
case_list_description: Samples with mutation data ({cases_seq}).
case_list_ids:\t{cases_seq_delim}"""

# File path
file_path = os.path.join(directory, "cases_sequenced.txt")

# Write data to the file
with open(file_path, 'w') as file:
    file.write(data)
del data

#write cases_cna files:
data = f"""cancer_study_identifier: {args.study}
stable_id: {args.study}_cna
case_list_name: Samples with CNA data
case_list_description: Samples with CNA data ({cases_cna}).
case_list_ids:\t{cases_cna_delim}"""

# File path
file_path = os.path.join(directory, "cases_cna.txt")

# Write data to the file
with open(file_path, 'w') as file:
    file.write(data)
del data

#write seq files:
data = f"""cancer_study_identifier: {args.study}
stable_id: {args.study}_cnaseq
case_list_name: Samples with mutation and CNA data
case_list_description: Samples with mutation and CNA data ({cases_cnaseq}).
case_list_ids:\t{cases_cnaseq_delim}"""

# File path
file_path = os.path.join(directory, "cases_cnaseq.txt")

# Write data to the file
with open(file_path, 'w') as file:
    file.write(data)
del data   

#write cases_rna_seq files:
data = f"""cancer_study_identifier: {args.study}
stable_id: {args.study}_rna_seq
case_list_name: Samples with RNAseq data
case_list_description: Samples with RNAseq data ({cases_rna_seq}).
case_list_ids:\t{cases_rna_seq_delim}"""

# File path
file_path = os.path.join(directory, "cases_rnaseq.txt")

# Write data to the file
with open(file_path, 'w') as file:
    file.write(data)
del data  

# Directory name with the current date and time
directory_name = f"Upload_{str_currentDateTime}"

# Create the directory if it doesn't exist
if not os.path.exists(directory_name):
    os.makedirs(directory_name)

# Move all files with the str_currentDateTime prefix into the directory
for file_name in os.listdir('.'):
    if file_name.startswith(str_currentDateTime) or file_name.startswith("meta"):
        shutil.move(file_name, os.path.join(directory_name, file_name))
        
# Move the "case_lists" directory into the new directory
if os.path.exists("case_lists"):
    shutil.move("case_lists", os.path.join(directory_name, "case_lists"))

# Copy all .xlsx files into the directory for records
shutil.move(excel_file, os.path.join(directory_name, os.path.basename(excel_file)))
print(f"Moved {excel_file} to {directory_name}")

print(f"Meta files generated and moved to {directory_name} successfully.")
print(f"Cancer study identifier: {args.study}")
print(f"Type of cancer: {args.cancer}")