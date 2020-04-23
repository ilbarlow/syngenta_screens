#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 14:03:57 2020

@author: ibarlow

Compile the syngenta library into new stock plates:
    1. Every drug in serial dilution from A1 - H12
    2. Import library and assign to a stock plate well

"""

import pandas as pd
from pathlib import Path
import numpy as np
import itertools
import math
import warnings

SYNGENTA_LIBRARY = Path('/Users/ibarlow/OneDrive - Imperial College London/'
                        + 'Documents/DrugScreening/DrugLibraries/Syngenta/'
                        + '/LibraryInformation/'
                        + 'SyngentaLibrary_strainscreen2020-MRC-10638.xlsx')

CONTROLS = {'DMSO': 1,
            'NoCompound': 1}
NO_CONTROLS = CONTROLS['DMSO'] + CONTROLS['NoCompound']

# 96well plate format
COLUMNS96WP = np.arange(1, 13)
ROWS96WP = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
WELLS96WP = [''.join([i[0], str(i[1])]) for i in
             list(itertools.product(ROWS96WP, COLUMNS96WP))]
save_to = SYNGENTA_LIBRARY.parent / '2020SygentaLibraryPlates3doses.csv'


if __name__ == '__main__':
    syngenta_drugs = pd.read_excel(SYNGENTA_LIBRARY,
                                   sheet_name='powder_and_liquid')
    syngenta_drugs['number_concentrations'] = [len(
                                        r.final_concentrations.split(','))
                                        for i, r in syngenta_drugs.iterrows()
                                        ]

    number_conditions = sum(syngenta_drugs['number_concentrations'])
    number_plates = math.ceil(number_conditions /
                              (len(WELLS96WP) - NO_CONTROLS))

    libraryDF = pd.DataFrame(columns=['drug_type',
                                      'drug_code',
                                      'drug_concentration'
                                      ])
    libraryDF['library_plate_number'] = sum([len(WELLS96WP)*[p] for p in
                                             range(1, number_plates+1)], [])
    libraryDF['well_name'] = WELLS96WP * number_plates

    # loop through the drugs and assign to the wells
    well_counter = 0
    for i, r in syngenta_drugs.iterrows():
        if libraryDF.loc[well_counter].well_name >= 'H7':
            libraryDF.loc[well_counter:well_counter+CONTROLS['DMSO']-1,
                          ['drug_type', 'drug_code']] = 'DMSO', 'DMSO'
            well_counter += CONTROLS['DMSO']
            libraryDF.loc[well_counter:well_counter+CONTROLS['NoCompound']-1,
                          ['drug_type', 'drug_code']] =\
                'NoCompound', 'NoCompound'

            well_counter += CONTROLS['NoCompound']

            libraryDF.loc[well_counter:well_counter+r.number_concentrations-1,
                          ['drug_type', 'drug_code']] = r.CSN, r.Code
            libraryDF.loc[well_counter:well_counter+r.number_concentrations-1,
                          'drug_concentration'] =\
                r.final_concentrations.split(',')[::-1]
            well_counter += r.number_concentrations

        else:
            libraryDF.loc[well_counter:well_counter+r.number_concentrations-1,
                          ['drug_type', 'drug_code']] = r.CSN, r.Code
            libraryDF.loc[well_counter:well_counter+r.number_concentrations-1,
                          'drug_concentration'] =\
                r.final_concentrations.split(',')[::-1]
            well_counter += r.number_concentrations

    if save_to.exists():
        warnings.warn('Sygenta 3 dose .csv file already exists')
    else:
        libraryDF.to_csv(save_to)
