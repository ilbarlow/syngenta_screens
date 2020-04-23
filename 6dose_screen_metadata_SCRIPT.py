#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 10:17:46 2020

@author: ibarlow

Script for compiling the metadata from the Syngenta screen done in December
2019

Parent directory (Auxiliary Files) contains individual folders for each
day of experiments, and contains the following files:
    YYYYMMDD_sourceplates.csv (DD = d-1 for day of tracking)
    YYYYMMDD_manual_metadata.csv (DD = day of tracking (d))
    YYYYMMDD_wormsorter.csv (DD = day of tracking)
    YYYMMDD_robot_bad_imaging_wells.csv (DD = d-1)

The following files will be generated: Nb. may be a good idea to put robot,
plate and merged metadata into a separate folder within AuxiliaryFiles
    YYYYMMDD_robot_metadata.csv (DD = d-1 for day of tracking):
        this is generated from the robot runlogs that are in a subdirectory
    YYYYMMDD_plate_metadata.csv:
        this is generated from the wormsorter files
    YYYYMMDD_merged_metadata.csv:
        this is a merge of the robot and wormsorter metadata
    YYYYMMDD_day_metadata.csv:
        compiled metadata from all the UI .csv files

"""

from pathlib import Path
import re
import pandas as pd
import numpy as np

from tierpsytools.hydra.compile_metadata import populate_96WPs,\
    merge_robot_metadata, merge_robot_wormsorter, get_day_metadata,\
    concatenate_days_metadata

date_regex = r"\d{8}"
PROJECT_DIRECTORY = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/SyngentaScreens_WFH/SyngentaScreen')

# %%
if __name__ == '__main__':

    day_root_dirs = [d for d in (PROJECT_DIRECTORY /
                                 'AuxiliaryFiles').glob("*")
                     if d.is_dir() and re.search(date_regex, str(d))
                     is not None]

    print('Calculating metadata for {} days of experiments'.format(
            len(day_root_dirs)))

    for count, day in enumerate(day_root_dirs):
        exp_date = re.findall(date_regex, str(day))[0]
        sourceplates_file = list(day.rglob("*_sourceplates.csv"))[0]
        assert (exp_date in str(sourceplates_file))
        manualmetadata_file = list(day.rglob('*_manual_metadata.csv'))[0]
        assert (exp_date in str(manualmetadata_file))
        wormsorter_file = list(day.rglob('*_wormsorter.csv'))[0]
        assert (exp_date in str(wormsorter_file))
        bad_wells_file = list(day.rglob('*_robot_bad_imaging_wells.csv'))[0]
        assert (exp_date in str(bad_wells_file))

        print('Collating manual metadata files: {},'.format(sourceplates_file)
              + '\n' '{}'.format(wormsorter_file))

#        plate_metadata = populate_96WPs(wormsorter_file,
#                                        del_if_exists=False)
#        robot_metadata = merge_robot_metadata(sourceplates_file,
#                                              del_if_exists=False)
        plate_file = list(day.rglob("*plate_metadata.csv"))[0]
        plate_metadata = pd.read_csv(plate_file)
        robot_file = list(day.rglob("*robot_metadata.csv"))[0]
        robot_metadata = pd.read_csv(robot_file)

        concat_metadata = merge_robot_wormsorter(day,
                                                 robot_metadata,
                                                 plate_metadata,
                                                 bad_wells_file,
                                                 del_if_exists=True)

        metadata_file = day / '{}_day_metadata.csv'.format(exp_date)

        print('Generating day metadata: {}'.format(
                metadata_file))

        day_metadata = get_day_metadata(concat_metadata,
                                        manualmetadata_file,
                                        saveto=metadata_file,
                                        del_if_exists=True)

        assert (day_metadata.shape[0] % 96 == 0)

        if AssertionError:
            # if there number of rows of day metadata is not divisible by
            # 96 means there has been an issue with propagating
            files_to_check = []
            plate_list = list(day_metadata['imaging_plate_id'].unique())
            day_metadata_grouped = day_metadata.groupby('imaging_plate_id')
            for plate in plate_list:
                _checking = day_metadata_grouped.get_group(plate)
                if _checking.shape[0] % 96 != 0:
                    wells = _checking['well_name'].unique()
                    print(wells)
                    for well in wells:
                        if (_checking['well_name'] == well).sum() != 3:
                            print(well)
                            files_to_check.append(
                                    _checking[_checking['well_name'] == well]
                                    ['imgstore_name'].to_list())
            files_to_check = [i for sublist in files_to_check for i in sublist]
            files_to_check = list(np.unique(files_to_check))
            print('These files need to be checked {}'.format(files_to_check))

# %%
    # combine all the metadata files
    concatenate_days_metadata(PROJECT_DIRECTORY / 'AuxiliaryFiles',
                              list_days=None,
                              saveto=None)

