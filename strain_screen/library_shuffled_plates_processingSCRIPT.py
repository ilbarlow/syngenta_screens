#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:45:50 2020

@author: ibarlow

Compile SyngentaStrainScreen shuffled library plates.

Make sourceplates for each robot run:
    - include robot runlog
    - source slot
    - drug type

Use robotlogs and template library plates that contained 3 doses of 100 drugs
across 3.5 x 96WP. No Compound and DMSO controls included on each plate and the
opentrons robot was used ot randomly shuffle each plate three times into
shuffled library plates
"""

import pandas as pd
from pathlib import Path
from tierpsytools.hydra.compile_metadata import merge_robot_metadata
import warnings

AUXILIARY_FILES = Path('/Users/ibarlow/' +
                         'OneDrive - Imperial College London/' +
                         'Documents/SyngentaScreens_WFH/' +
                         'SyngentaStrainScreen/' +
                         'AuxiliaryFiles')
PREPROCESSING_REQUIRED = False
sourceplate_files = list(
            AUXILIARY_FILES.rglob('2020SygentaLibrary3doses*'))
sourceplate_files = [i for i in sourceplate_files if 'shuffled' not in str(i)]

if __name__ == '__main__':
    if PREPROCESSING_REQUIRED == True:
        robot_logs = list(AUXILIARY_FILES.rglob('*runlog.csv'))

        # robot runlogs have extra rows for when the robot was mixing the drug
        # by pipetting. Need to remove these
        for count, file in enumerate(robot_logs):
            rlog = pd.read_csv(file)
            rlog = rlog.drop(rlog[
                    rlog['source_slot'] == rlog['destination_slot']
                    ].index)
            rlog.to_csv(file.parent / (file.stem + '_clean.csv'), index=False)

        # update the sourceplates file
        for file in sourceplate_files:
            splate = pd.read_csv(file)
            splate.loc[:, 'robot_runlog_filename'] = [r[
                                                'robot_runlog_filename'].replace(
                                                'runlog.csv', 'runlog_clean.csv')
                                                for i, r in
                                                splate.iterrows()
                                                ]
            warnings.warn('{} file being edited to update robot log'.format(file))
            splate.to_csv(file, index=False)

    for file in sourceplate_files:
        robot_metadata = merge_robot_metadata(file,
                                              saveto=None,
                                              del_if_exists=True,
                                              compact_drug_plate=True,
                                              drug_by_column=False)
        robot_metadata.sort_values(by=['source_plate_number',
                                       'destination_well'],
                                   ignore_index=True,
                                   inplace=True)

        robot_metadata['shuffled_plate_id'] = [r.source_plate_id +
                                               '_sh%02d' %(r.robot_run_number)
                                               for i, r in
                                               robot_metadata.iterrows()]
        robot_metadata['is_bad_well'].fillna(False, inplace=True)
        robot_metadata.to_csv(str(file).replace('.csv', '_shuffled.csv'),
                              index=False)
        del robot_metadata