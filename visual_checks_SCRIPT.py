#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 08:39:34 2020

@author: ibarlow

Script to find the videos to check that drugs effects can be see by eye for
Spirindolines and levamisol

nAChR agonist	lev-1 receptor	CSAA020364
vAChT inhibitor	spiroindoline type	CSAC270548
vAChT inhibitor	spiroindoline type	CSCC202642
vAChT inhibitor	spiroindoline type	CSCC230937

Output a spreadsheet with imgstore name, prestim, wells and drug type, 
drug_concentration

Only check N2s

"""

from pathlib import Path
import pandas as pd

metadata_fname = Path('/Users/ibarlow/Imperial College London/Minga, Eleni - SyngentaStrainScreen/AuxiliaryFiles/metadata.csv')
out_fname = metadata_fname.parent / 'drug_effects_to_check.csv'

if out_fname.exists():
    print ('Warning: overwriting ')
    
worm_strain = 'N2'

drug_dict = {'CSAA020364': 'nAChR agonist; paralysing',
             'CSAC270548': 'vAChT inhibitor; coiling',
             'CSCC202642': 'vAChT inhibitor; coiling',
             'CSCC230937': 'vAChT inhibitor; coiling'}

if __name__ == '__main__':
    metadata = pd.read_csv(metadata_fname)
        
    columns_to_keep = [i for i in metadata.columns if 'bad' in i]
    columns_to_keep.extend([i for i in metadata.columns if 'drug' in i])
    columns_to_keep.extend(['date_yyyymmdd',
                            'well_name',
                           'imgstore_name',
                           'worm_strain',
                           'imaging_run_number',
                           'camera_serial',
                           'drug_mode',
                           ])
    columns_to_keep = columns_to_keep[::-1]
    assert len(set(columns_to_keep).intersection(metadata.columns)) == len(columns_to_keep)-1
    
    #now extract out only the interesting drugs
    files_to_check = []
    for dkey in drug_dict:
        files_to_check.append(metadata[metadata['drug_type'] == dkey])
            
    files_to_check = pd.concat(files_to_check)
    files_to_check['drug_mode'] = files_to_check['drug_type'].map(drug_dict)
    
    # filter these files, first only keep the prestim files
    files_to_check['stim_type'] = [r['imgstore_name'].split('_')[-3]
                                   for i, r in files_to_check.iterrows()]
    files_to_check = files_to_check[files_to_check['stim_type']=='prestim']
    files_to_check = files_to_check[files_to_check['worm_strain'] == worm_strain]
    
    unique_files = list(files_to_check.imgstore_name.unique())
    
    files_to_check_grouped = files_to_check.groupby('imgstore_name')
    
    #pull into a dataframe
    files_to_check_concise = []
    for fname in unique_files:
        _temp = pd.DataFrame()
        for col in columns_to_keep:
            if files_to_check_grouped.get_group(fname)[col].unique().shape[0]>1:
                try:
                    _temp[col] = files_to_check_grouped.get_group(fname)[col].str.cat(
                        sep='; ')
                except AttributeError:
                    _temp[col] = [files_to_check_grouped.get_group(fname)[col].to_list()]
            else:
                _temp[col] = files_to_check_grouped.get_group(fname)[col].unique()
        _temp
        files_to_check_concise.append(_temp)
        del _temp
        
    files_to_check_concise = pd.concat(files_to_check_concise)
    
    files_to_check_concise.to_csv(out_fname, index=False)
    
    
    
    
    
    
    
