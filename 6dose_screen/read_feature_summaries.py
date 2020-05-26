#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a script for the initial inspection of the syngentascreen data.

Created on Tue Apr  7 15:02:35 2020

@author: em812
"""

import pandas as pd
from pathlib import Path
from tierpsytools.read_data.hydra_metadata import read_hydra_metadata, align_bluelight_conditions
from tierpsytools.analysis import significant_features
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

control_drug = ['DMSO', 'NoCompound']
control_strain = 'N2'
bad_feat_thresh = 3 # 3 standard deviations away from the mean

# %% Input
saveto = Path('.').resolve() / 'Documents/GitHub/pythonScripts/SyngentaScreens/6dose_screen' / 'visual_checks'
saveto.mkdir(parents=True, exist_ok=True)

root = Path('/Users/ibarlow/Imperial College London/Minga, Eleni - SyngentaScreen')
# root = Path('/Users/em812/OneDrive - Imperial College London/share/SyngentaScreen')
feat_file = root / 'SummaryFiles' / 'features_summaries_compiled.csv'
fname_file = root / 'SummaryFiles' / 'filenames_summaries_compiled.csv'

metadata_file = root / 'AuxiliaryFiles' / 'full_metadata.csv'
# moa_file = './AllCompoundsMoA.csv'
moa_file = root / 'analysis' / 'AllCompoundsMoA.csv'

bad_well_cols = ['is_bad_well_from_gui', 'is_bad_well_misplaced_plate', 'is_bad_well_ledfailure']

#%% Read data
feat = pd.read_csv(feat_file, comment='#')
fname = pd.read_csv(fname_file, comment='#')

meta = pd.read_csv(metadata_file, index_col=None)

# %% Match metaddata and features
feat, meta = read_hydra_metadata(feat, fname, meta)
feat, meta = align_bluelight_conditions(feat, meta)

# %% Add moa info
moa = pd.read_csv(moa_file, index_col=None)
moa = moa.rename(columns={"CSN": "drug_type"})
meta = pd.merge(
    meta, moa[['MOA_general', 'MOA_specific', 'drug_type', 'MOA_group']],
    on='drug_type', how='left')

# %% Preprocess (optional)
# Remove bad wells
bad = meta[bad_well_cols].any(axis=1)

feat = feat.loc[~bad,:]
meta = meta.loc[~bad,:]

# fill in nans for DMSO and NoCompound where concentrations are missing
DMSO_locs = meta[meta['drug_type'] == 'DMSO'].index
meta.loc[DMSO_locs, 'imaging_plate_drug_concentration'] = 0.1

NC_locs = meta[meta['drug_type'] == 'NoCompound'].index
meta.loc[NC_locs, 'imaging_plate_drug_concentration'] = 0

# Remove wells missing bluelight conditions
imgst_cols = [col for col in meta.columns if 'imgstore_name' in col]
miss = meta[imgst_cols].isna().any(axis=1)

feat = feat.loc[~miss,:]
meta = meta.loc[~miss,:]

# Remove features with too many nans
feat = feat.loc[:, feat.isna().sum(axis=0)/feat.shape[0]<0.10]
# Impute remaining nans
feat = feat.fillna(feat.mean())

prestim_feats = [x for x in feat.columns if 'prestim' in x]

# %% look at coiling and paralysed worms
coiled = 'spiroindoline-type'
paralysed = 'lev-1 receptor'

# select only coiled worms
coiled_meta = [meta[meta['MOA_specific']==coiled]]
for control in control_drug:
    coiled_meta.append(meta[meta['drug_type'] == control])
coiled_meta = pd.concat(coiled_meta)

coiled_meta = coiled_meta[coiled_meta['worm_strain'] == control_strain]

# filter down the feat mat 
coiled_feat = feat.loc[coiled_meta.index, prestim_feats]

# now do statistical tests 
y_classes = list(zip(coiled_meta.drug_type,
                     coiled_meta.imaging_plate_drug_concentration))
y_classes = np.array(['{}, {}'.format(x[0], x[1]) for x in y_classes])
sig_feats = significant_features.k_significant_feat(coiled_feat,
                                                    y_classes,
                                                    k=20)

speed_feats = ['speed_10th_prestim',
                'speed_50th_prestim',
                'speed_90th_prestim',
                'speed_IQR_prestim']
curvature_feats = ['curvature_midbody_abs_10th_prestim',
                    'curvature_midbody_abs_50th_prestim',
                    'curvature_midbody_abs_90th_prestim',
                    'curvature_midbody_abs_IQR_prestim',
                    'curvature_midbody_norm_abs_10th_prestim',
                    'curvature_midbody_norm_abs_50th_prestim',
                    'curvature_midbody_norm_abs_90th_prestim',
                    'curvature_midbody_norm_abs_IQR_prestim',
                    'curvature_midbody_w_backward_abs_10th_prestim',
                    'curvature_midbody_w_backward_abs_50th_prestim',
                    'curvature_midbody_w_backward_abs_90th_prestim',
                    'curvature_midbody_w_backward_abs_IQR_prestim',
                    'curvature_midbody_w_forward_abs_10th_prestim',
                    'curvature_midbody_w_forward_abs_50th_prestim',
                    'curvature_midbody_w_forward_abs_90th_prestim',
                    'curvature_midbody_w_forward_abs_IQR_prestim',
                    'curvature_midbody_w_paused_abs_10th_prestim',
                    'curvature_midbody_w_paused_abs_50th_prestim',
                    'curvature_midbody_w_paused_abs_90th_prestim',
                    'curvature_midbody_w_paused_abs_IQR_prestim']


plotting_df = pd.concat([coiled_meta, coiled_feat], axis=1)

unique_drugs = list(plotting_df['drug_type'].unique())
    
for c, f in enumerate(sig_feats):
    plt.figure()
    sns.swarmplot(x = 'drug_type',
                  y = feat,
                  data = plotting_df,
                  hue = 'imaging_plate_drug_concentration',
                  palette = 'BuGn_r')
    plt.xticks(rotation=45)
    plt.savefig(saveto / '{}_swarmplot.png'.format(feat), dpi=600)
    plt.close(c)
         
    
for c, f in enumerate(speed_feats):
    plt.figure(figsize=[6,4])
    sns.swarmplot(x = 'drug_type',
                  y = f,
                  data = plotting_df,
                  hue = 'imaging_plate_drug_concentration',
                  palette = 'BuGn')
    plt.xticks(rotation=45)
    plt.legend(loc='center left',
               bbox_to_anchor=(1.0,0.5),
               ncol = 1,
               frameon= True)
    plt.tight_layout()
    plt.savefig(saveto / '{}_spiroindolines_swarmplot.png'.format(f), dpi=300)
    plt.close(c+1)

for c, f in enumerate(curvature_feats):
    plt.figure(figsize=[6,4])
    sns.swarmplot(x = 'drug_type',
                  y = f,
                  data = plotting_df,
                  hue = 'imaging_plate_drug_concentration',
                  palette = 'BuGn')
    plt.xticks(rotation=45)
    plt.legend(loc='center left',
               bbox_to_anchor=(1.0,0.5),
               ncol = 1,
               frameon= True)
    plt.tight_layout()
    plt.savefig(saveto / '{}_spiroindolines_swarmplot.png'.format(f), dpi=300)
    plt.close(c+1)
    
# %% levamisol checks
paralysed_meta = [meta[meta['MOA_specific']==paralysed]]
for control in control_drug:
    paralysed_meta.append(meta[meta['drug_type'] == control])
paralysed_meta = pd.concat(paralysed_meta)

paralysed_meta = paralysed_meta[paralysed_meta['worm_strain'] == control_strain]

# filter down the feat mat 
paralysed_feat = feat.loc[paralysed_meta.index,
                          prestim_feats]

# now do statistical tests 
y_classes = list(zip(paralysed_meta.drug_type,
                     paralysed_meta.imaging_plate_drug_concentration))
y_classes = np.array(['{}, {}'.format(x[0], x[1]) for x in y_classes])

sig_feats = significant_features.k_significant_feat(paralysed_feat,
                                                    y_classes,
                                                    k=20)


plotting_df = pd.concat([paralysed_meta,
                         paralysed_feat],
                        axis=1)


for c, f in enumerate(speed_feats):
    plt.figure(figsize=[6,4])
    sns.swarmplot(x = 'drug_type',
                  y = f,
                  data = plotting_df,
                  hue = 'imaging_plate_drug_concentration',
                  palette = 'BuGn')
    plt.xticks(rotation=45)
    plt.legend(loc='center left',
               bbox_to_anchor=(1.0,0.5),
               ncol = 1,
               frameon= True)
    plt.tight_layout()
    plt.savefig(saveto / '{}_paralysed_swarmplot.png'.format(f), dpi=300)
    plt.close(int(c+1))
    
for c, f in enumerate(curvature_feats):
    plt.figure(figsize=[6,4])
    sns.swarmplot(x = 'drug_type',
                  y = f,
                  data = plotting_df,
                  hue = 'imaging_plate_drug_concentration',
                  palette = 'BuGn')
    plt.xticks(rotation=45)
    plt.legend(loc='center left',
               bbox_to_anchor=(1.0,0.5),
               ncol = 1,
               frameon= True)
    plt.tight_layout()
    plt.savefig(saveto / '{}_paralysed_swarmplot.png'.format(f), dpi=300)
    plt.close(int(c+1))


# g = sns.FacetGrid(plotting_df, col='drug_type')
# g = g.map(plt.scatter, 'imaging_plate_drug_concentration', sig_feats[0])

#%% 
# TODO: check out the videos for the DMSO and no compound controls that are 
# more than 2 SD from the control mean for some of the selected features

# select only the control data
control_meta = []
for c in control_drug:
    control_meta.append(meta[meta['drug_type'] == c])
control_meta = pd.concat(control_meta)
control_meta = control_meta[control_meta['worm_strain'] == control_strain]

control_feat = feat.loc[control_meta.index, prestim_feats]

# find wells above and below thresholds
upper_thresh = control_feat > (control_feat.mean(axis=0)
                               + bad_feat_thresh * control_feat.std(axis=0))
lower_thresh = control_feat < (control_feat.mean(axis=0)
                               - bad_feat_thresh * control_feat.std(axis=0))

flagged_df = upper_thresh | lower_thresh
flagged_df['total_counts'] = flagged_df.sum(axis=1)

#make some plots
plt.figure()
sns.distplot(flagged_df['total_counts'],
             bins = 100)
plt.xlim([0, 200])
plt.xlabel('number features per well > {} standard devations from mean'.format(bad_feat_thresh))
plt.savefig(saveto / 'flagged_wells_histogram.png')

plt.figure()
sns.distplot(flagged_df[prestim_feats].sum(axis=0),
             color = 'green')
plt.xlabel('number_wells per feature > {} standard deviations from mean'.format(bad_feat_thresh))
plt.savefig(saveto / 'flagged_features_hisogram.png')

files_to_check = flagged_df[flagged_df['total_counts'] >
                            (flagged_df['total_counts'].mean() +
                             flagged_df['total_counts'].std())]

files_to_check = meta.loc[files_to_check.index,:]
files_to_check[['imgstore_name_prestim', 'well_name']].to_csv(
    saveto / 'files_to_check_{}_{}std.csv'.format(control_strain,
                                                  bad_feat_thresh),
    index=True)

meta.loc[np.random.choice(flagged_df.index, 100), ['imgstore_name_prestim',
                                                   'well_name']].to_csv(
    saveto / 'files_to_check_{}_random.csv'.format(control_strain),
    index=True)

                                                       
# %% differences between N2 and hawaiian
control_drug= 'NoCompound'
strains = ['N2', 'CB4856']

strain_meta = []
for strain in strains:
    strain_meta.append(meta[meta['worm_strain'] == strain])
strain_meta = pd.concat(strain_meta)
strain_meta = strain_meta[strain_meta['drug_type'] == 'NoCompound']

strain_feat = feat.loc[strain_meta.index, prestim_feats]

# now do statistical tests 
strain_classes = np.array(strain_meta.worm_strain)
strain_feats = significant_features.k_significant_feat(strain_feat,
                                                       strain_classes,
                                                       k=20)

plotting_df = pd.concat([strain_feat,
                         strain_meta],
                        axis=1)
# some plots
for c, f in enumerate(speed_feats):
    plt.figure(figsize=[6,4])
    sns.swarmplot(x = 'worm_strain',
                  y = f,
                  data = plotting_df,
                  palette = 'BuGn_r')
    plt.xticks(rotation=45)
    # plt.legend(loc='center left',
    #            bbox_to_anchor=(1.0,0.5),
    #            ncol = 1,
    #            frameon= True)
    # plt.tight_layout()
    plt.savefig(saveto / '{}_N2_hawaiian_swarmplot.png'.format(f), dpi=300)
    plt.close(int(c+1))
    
# some plots
for c, f in enumerate(curvature_feats):
    plt.figure(figsize=[6,4])
    sns.swarmplot(x = 'worm_strain',
                  y = f,
                  data = plotting_df,
                  palette = 'BuGn_r')
    plt.xticks(rotation=45)
    # plt.legend(loc='center left',
    #            bbox_to_anchor=(1.0,0.5),
    #            ncol = 1,
    #            frameon= True)
    # plt.tight_layout()
    plt.savefig(saveto / '{}_N2_hawaiian_swarmplot.png'.format(f), dpi=300)
    plt.close(int(c+1))