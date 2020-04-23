#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:22:27 2019

@author: ibarlow
"""

""" create a template for expanding out the syngenta library into 96 well plates
suitable for drug experiments

Requirements:
    1. Each compounds made into 6 concentrations from stock (0.1, 1, 10, 50, 100, 150 mM)
    2. Compounds stored at stock concentration in a column of a 96 WP
    3. STOCK plates can be used to make SOURCE plates
    4. SOURCE plates are used to make DESTINATION plates via a random shuffle across the plate

    Terminology:
        LIBRARY plate = 96WP containing compounds at highest stock concentration
        STOCK plate = 96 well plate with columns of different concentration of compounds stored in DMSO
        SOURCE plate = 96 well plate with columns of different concentrations of compounds in DMSO + water
        DESTINATION plate = square well 96 well plate used for tracking worms

        """

import pandas as pd
import numpy as np
from pathlib import Path
import itertools
import math
import re

#specify drug library
SyngentaLibrary = Path('/Users/ibarlow/OneDrive - Imperial College London/Documents/DrugScreening/DrugLibraries/Syngenta/SygentaLibrary_stocks2019.xlsx')

#other parameters
CONTROLS = {'DMSO':1,
            'NoCompound':12}
CONCENTRATIONS = { '6doses' : [150,100,30,10,3,1],
                   '5doses' : [100,30,10,3,1],
                   '4doses' : [30,10,3,1],
                   '3doses' : [100, 10, 1]}
#96well plate format
COLUMNS96WP = np.arange(1,13)
ROWS96WP = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
WELLS96WP = [''.join([i[0], str(i[1])]) for i in list(itertools.product(ROWS96WP, COLUMNS96WP))]

n_controls=len(CONTROLS)
#n_concentrations = len(concentrations)

#load drug library
libraryDF = pd.read_excel(SyngentaLibrary, sheet_name=list(CONCENTRATIONS.keys()))
libraryDF= pd.concat(libraryDF, sort=True)
libraryDF.reset_index(level=0, inplace=True)
libraryDF.rename(columns = {'level_0':'no_doses'}, inplace=True)
libraryDF.reset_index(drop=True, inplace=True)

#group the dataframe by the number of doses and then loop through to find the
#ones that need to made manually
libraryDF_grouped = libraryDF.groupby('no_doses')
dose_groups = list(libraryDF_grouped.groups.keys())
for d in dose_groups:
    if d in list(CONCENTRATIONS.keys()):
        _drugDF = libraryDF_grouped.get_group(d)
        drugs_to_sort = _drugDF[_drugDF['Prediluted']==1]
        drugs_for_opentrons = _drugDF[_drugDF['Prediluted']==0]

        n_drugs = drugs_to_sort.shape[0]
        n_concentrations = len(CONCENTRATIONS[d])

        print ('{} drugs being sorted into {}: {} mM'.format(n_drugs, d, CONCENTRATIONS[d]))

        n_plates_required = math.ceil((n_drugs*n_concentrations)/(COLUMNS96WP.shape[0]-n_controls))

        print ('{} plates required'.format(n_plates_required))

        WP96dict = {p:COLUMNS96WP for p in range(1,n_plates_required+1)}

        plate_iterator = [(k,t) for k,v in WP96dict.items() for t in v]

        #now fill a dataframe with plates containing drugs
        stockDF = pd.DataFrame(columns = ['drug_type',
                                          'drug_concentration'])#columns = ['plate',
#                                   'well',
#                                   'source_column',
#                                   'drug_type',
#                                   'drug_concentration'])
        stockDF['plate'] = sum([96*[p] for p in range(1,n_plates_required+1)], []) #slow if too long
        stockDF['well'] = n_plates_required*WELLS96WP
        stockDF['source_column'] = [int(re.findall(r'\d+', r.well)[0]) for i,r in stockDF.iterrows()]

        #assign the control columns
        stockDF['drug_type'].iloc[stockDF[stockDF['source_column']==CONTROLS['DMSO']].index] = 'DMSO'
        stockDF['drug_type'].iloc[stockDF[stockDF['source_column']==CONTROLS['NoCompound']].index] = 'NoCompound'

        stockDF = stockDF.sort_values(by=['plate', 'source_column']).reset_index(drop=True)
        stockDF_grouped = stockDF.groupby(['plate', 'source_column'])

        #iterate over the plates to assign the stock plates
        SourcePlatesDF = pd.DataFrame()
        plate_counter=0
        for i,r in drugs_to_sort.iterrows():
            print (r.CSN)
            plate = plate_iterator[plate_counter]
            concentration_counter = 0
            while concentration_counter<n_concentrations:
                c= CONCENTRATIONS[d][concentration_counter]
                plate = plate_iterator[plate_counter]
                #        print (plate, r.CSN, c)
                if stockDF_grouped.get_group(plate).drug_type.isna().sum()>0:
                    SourcePlatesDF = SourcePlatesDF.append(stockDF_grouped.get_group(plate).assign(
                                drug_type=r.CSN,
                                drug_concentration=c,
                                MOA_general = r['MOA general'],
                                MOA_specific =r['MOA specific'],
                                MW = r.MW,
                                CODE = r.Code), sort=True)
                    print(plate, r.CSN, c)
                    concentration_counter += 1
                    plate_counter+=1

                else:
                    SourcePlatesDF = SourcePlatesDF.append(stockDF_grouped.get_group(plate), sort=True)
                    print (plate, 'control')
                    plate_counter +=1

        SourcePlatesDF.to_csv(SyngentaLibrary.parent / '{}_Syngenta_manualsourceplates.csv'.format(d), index=False)
        SourcePlatesDF.drop_duplicates(subset=['plate', 'source_column']).to_csv(SyngentaLibrary.parent / '{}_Syngenta_manual_sourceplates_abridged_summary.csv'.format(d), index=False)
        drugs_for_opentrons.to_csv(SyngentaLibrary.parent / '{}_Syngenta_forOpentrons.csv'.format(d), index=False)

#%%
