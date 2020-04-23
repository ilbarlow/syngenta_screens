#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:09:23 2020

@author: ibarlow

Script for making the day metadata for the January Syngenta screen of 100
drugs at 3 doses on 12 strains

Parent directory (Auxiliary Files) contains individual folders for each
day of experiments, and contains the following files:
    YYYYMMDD_manual_metadata.csv (DD = day of tracking (d))
    YYYYMMDD_wormsorter.csv (DD = day of tracking)
    YYYMMDD_robot_bad_imaging_wells.csv (DD = d-1)
    
    A folder called sourceplates containing the shuffled sourceplates


The following files will be generated: Nb. may be a good idea to put robot,
plate and merged metadata into a separate folder within AuxiliaryFiles
    YYYYMMDD_plate_metadata.csv:
        this is generated from the wormsorter files
    YYYYMMDD_day_metadata.csv:
        compiled metadata from all the UI .csv files
    
    metadata.csv - final concatenated metadata of all the days combined


"""

import pandas as pd
from pathlib import Path
from tierpsytools.hydra.compile_metadata import populate_96WPs,\
    convert_bad_wells_lut, get_day_metadata, day_metadata_check,\
    concatenate_days_metadata
import re
,
date_regex = r"\d{8}"
plate_id_regex= r""
PROJECT_DIRECTORY = Path('/Users/ibarlow/' +
                         'OneDrive - Imperial College London/' +
                         'Documents/' +
                         'SyngentaScreens_WFH/' +
                         'SyngentaStrainScreen/' +
                         'AuxiliaryFiles')

# import the shuffled plates as a reference
sourceplates = list(PROJECT_DIRECTORY.rglob('*shuffled.csv'))

drug_plates = []
for file in sourceplates:
    drug_plates.append(pd.read_csv(file))
drug_plates = pd.concat(drug_plates)

# %% 
if __name__=='__main__':
    day_root_dirs = [d for d in PROJECT_DIRECTORY.glob("*")
                     if d.is_dir() and re.search(date_regex, str(d))
                     is not None]
    
    print('Calculating metadata for {} days of experiments'.format(
            len(day_root_dirs)))
    
    for count, day in enumerate(day_root_dirs):
        exp_date = re.findall(date_regex, str(day))[0]
        manualmetadata_file = list(day.rglob('*_manual_metadata.csv'))[0]
        assert (exp_date in str(manualmetadata_file))
        wormsorter_file = list(day.rglob('*_wormsorter.csv'))[0]
        assert (exp_date in str(wormsorter_file))
        bad_wells_file = list(day.rglob('*_robot_bad_imaging_wells.csv'))[0]
        assert (exp_date in str(bad_wells_file))
        
        print('Collating wormsorter files: {}'.format(wormsorter_file))
        
        plate_metadata = populate_96WPs(wormsorter_file,
                                        del_if_exists=True)
        
        bad_wells_df = convert_bad_wells_lut(bad_wells_file)
        
        plate_metadata = pd.merge(plate_metadata,
                                  bad_wells_df,
                                  on=['imaging_plate_id', 'well_name'],
                                  how='outer')
        
        metadata_file = day / '{}_day_metadata.csv'.format(exp_date)
    
        print('Generating day metadata: {}'.format(
                metadata_file))
        
        
        try:
            day_metadata = get_day_metadata(plate_metadata,
                                            manualmetadata_file,
                                            saveto=metadata_file,
                                            del_if_exists=True,
                                            include_imgstore_name=True)
        except ValueError:
            print ('imgstore error')
            day_metadata = get_day_metadata(plate_metadata,
                                            manualmetadata_file,
                                            saveto=metadata_file,
                                            del_if_exists=True,
                                            include_imgstore_name=False)
        
        day_metadata = pd.merge(day_metadata,
                               drug_plates,
                               left_on=['source_plate_id',
                                        'well_name'],
                               right_on=['shuffled_plate_id',
                                         'destination_well'],
                               suffixes=('_day', '_robot'),
                               how='outer')
        day_metadata.drop(
            day_metadata[day_metadata['imaging_plate_id'].isna()].index,
                        inplace=True)
        
        files_to_check = day_metadata_check(day_metadata, day, plate_size=48)
 
        day_metadata.to_csv(metadata_file, index=False)
    
    # %%concatenate day metadata
    
    # combine all the metadata files
    concatenate_days_metadata(PROJECT_DIRECTORY,
                              list_days=None,
                              saveto=None)
    
