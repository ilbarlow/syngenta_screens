#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:20:42 2020

@author: ibarlow

Script to combine the bad well information with the sourceplates files

"""

import pandas as pd
from pathlib import Path

ROBOT_STOCKPLATES = Path('/Volumes/behavgenom$/Ida/Data/Hydra/SyngentaScreen/'
                         + 'AuxiliaryFiles/'
                         + 'Syngenta_robot_library_compiled.xlsx')
MANUAL_STOCKPLATES = Path('/Volumes/behavgenom$/Ida/Data/Hydra/SyngentaScreen/'
                          + 'AuxiliaryFiles/'
                          + 'Syngenta_manual_library_compiled.xlsx')

SOURCEPLATE_FILES = list(ROBOT_STOCKPLATES.parent.rglob("*_sourceplates.csv"))

if __name__ == '__main__':
    robot_df = pd.read_excel(ROBOT_STOCKPLATES)
    manual_df = pd.read_excel(MANUAL_STOCKPLATES)

    stockplate_df = pd.concat([robot_df, manual_df], axis=0, sort=True)

    for c, file in enumerate(SOURCEPLATE_FILES):
        print('Updating sourceplate file: {}'.format(file))
        _sp = pd.read_csv(file)
        _sp = pd.merge(_sp,
                       stockplate_df[['stock_plate_id',
                                      'column',
                                      'bad_wells'
                                      ]],
                       on=['stock_plate_id', 'column'])
        _sp.to_csv(file, index=False)
