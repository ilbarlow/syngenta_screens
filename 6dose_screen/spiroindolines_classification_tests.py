#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:36:24 2020

@author: ibarlow

Read feature summaries from the December SyngentaScreen and try to run
a classifier that can reliably classify spiroindoline and levamisol
vs control (DMSO) for N2 and Hawaiian worms

"""

import pandas as pd
from pathlib import Path
from tierpsytools.read_data.hydra_metadata import read_hydra_metadata, align_bluelight_conditions
from tierpsytools.analysis import significant_features
from tierpsytools.analysis.decomposition import plot_pca
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

from sklearn.decomposition import PCA
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import StratifiedKFold
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import LabelEncoder
import random


CONTROL_DRUGS = ['DMSO', 'NoCompound']
CONTROL_STRAIN = 'CB4856' #'N2'
BAD_FEAT_THRESH = 3 # 3 standard deviations away from the mean
BAD_FEAT_FILTER = 0.1 # 10% threshold for removing bad features
BAD_WELL_FILTER = 0.3 # 30% threshold for bad well

# %% Input
SAVETO = Path('.').resolve() / 'Documents/GitHub/pythonScripts/SyngentaScreens/6dose_screen' / 'classifier_tests'
SAVETO.mkdir(parents=True, exist_ok=True)

root = Path('/Users/ibarlow/Imperial College London/Minga, Eleni - SyngentaScreen')
# root = Path('/Users/em812/OneDrive - Imperial College London/share/SyngentaScreen')
feat_file = root / 'SummaryFiles' / 'features_summaries_compiled.csv'
fname_file = root / 'SummaryFiles' / 'filenames_summaries_compiled.csv'

metadata_file = root / 'AuxiliaryFiles' / 'full_metadata.csv'
# moa_file = './AllCompoundsMoA.csv'
moa_file = root / 'analysis' / 'AllCompoundsMoA.csv'

bad_well_cols = ['is_bad_well_from_gui',
                 'is_bad_well_misplaced_plate',
                 'is_bad_well_ledfailure']

#%% Read data
feat = pd.read_csv(feat_file,
                   comment='#')
fname = pd.read_csv(fname_file,
                    comment='#')

meta = pd.read_csv(metadata_file, index_col=None)

# %% Match metaddata and features
feat, meta = read_hydra_metadata(feat,
                                 fname,
                                 meta)
feat, meta = align_bluelight_conditions(feat,
                                        meta)

# %% Add moa info
moa = pd.read_csv(moa_file,
                  index_col=None)
moa = moa.rename(columns={"CSN": "drug_type"})
meta = pd.merge(meta,
                moa[['MOA_general',
                     'MOA_specific',
                     'drug_type',
                     'MOA_group']],
                on='drug_type',
                how='left'
                )

# %% Preprocess (optional)
# Remove bad wells
bad = meta[bad_well_cols].any(axis=1)

feat = feat.loc[~bad,:]
meta = meta.loc[~bad,:]

# fill in nans for DMSO and NoCompound where concentrations are missing
DMSO_locs = meta[meta['drug_type'] == 'DMSO'].index
meta.loc[DMSO_locs, 'imaging_plate_drug_concentration'] = 0.1

meta['drug_type'].fillna('NoCompound', inplace=True)
NC_locs = meta[meta['drug_type'] == 'NoCompound'].index
meta.loc[NC_locs, 'imaging_plate_drug_concentration'] = 0

# Remove wells missing bluelight conditions
imgst_cols = [col for col in meta.columns if 'imgstore_name' in col]
miss = meta[imgst_cols].isna().any(axis=1)

feat = feat.loc[~miss,:]
meta = meta.loc[~miss,:]

# %% exploratory look to see if high number of 'nan's correlate with outliers
from kneed import KneeLocator

# find wells above and below thresholds
upper_thresh = feat > (feat.mean(axis=0)
                       + BAD_FEAT_THRESH * feat.std(axis=0))
lower_thresh = feat < (feat.mean(axis=0)
                       - BAD_FEAT_THRESH * feat.std(axis=0))

flagged_df = upper_thresh | lower_thresh
flagged_df['total_nan_counts'] = feat.isna().sum(axis=1)
flagged_df['total_outlier_counts'] = flagged_df.sum(axis=1)
flagged_df['percentage_outlier_nonnan_counts'] = [r.total_outlier_counts / 
                                             (flagged_df.shape[1] - 
                                              r.total_nan_counts) for i, r
                                             in flagged_df.iterrows()]
sns.regplot(x='total_nan_counts',
            y='percentage_outlier_nonnan_counts',
            data=flagged_df)
plt.savefig(SAVETO.parent / 'outliers_vs_nanvalues.png')

sns.regplot(x='percentage_outlier_nonnan_counts',
            y='total_nan_counts',
            data=flagged_df,
            logx=True)

plt.scatter(flagged_df['total_nan_counts'],
            np.log2(flagged_df['percentage_outlier_nonnan_counts']))

knee_calculator_df = pd.DataFrame()
knee_calculator_df['total_nan_counts'] = flagged_df['total_nan_counts'].unique()
knee_calculator_df['av_outlier_nonnan'] = [
    flagged_df[flagged_df['total_nan_counts'] ==
               r.total_nan_counts]['percentage_outlier_nonnan_counts'].mean()
    for i,r in knee_calculator_df.iterrows()] 
knee_calculator_df.sort_values(by='total_nan_counts', inplace=True)

# need to make a plot of the          
kneedle = KneeLocator(x=knee_calculator_df['total_nan_counts'].values,
                      y=knee_calculator_df['av_outlier_nonnan'].values,
                      S=1.0,
                      curve='convex',
                      direction='increasing')
kneedle.plot_knee()
plt.xlabel('total_nan_counts')
plt.ylabel('average_outliers_nonnan')
plt.savefig(SAVETO.parent / 'kneeplot.png')

NAN_CUTOFF = round(kneedle.find_knee()[0], -3)

# %% Filter well with too many nans
NAN_CUTOFF = 7000
print('{} features for nan cutoff'.format(NAN_CUTOFF))

# remove wells with too many nans as likely to have lots of outliers
feat = feat[feat.isna().sum(axis=1) < NAN_CUTOFF]
meta = meta.loc[feat.index,:]

# Remove features with too many nans
feat = feat.loc[:, feat.isna().sum(axis=0)/feat.shape[0] < BAD_FEAT_FILTER]

# Impute remaining nans
feat = feat.fillna(feat.mean(axis=0))
 
prestim_feats = [x for x in feat.columns if 'prestim' in x]
bluelight_feats = [x for x in feat.columns if 'bluelight' in x]
pathcurvature_feats = [x for x in feat.columns if 'path_curvature' in x]

prestim_nopath_feats = list(set(prestim_feats) - set(pathcurvature_feats))

# %% remove outlier DMSO data
# select only the control data for N2s
control_meta = []
for c in CONTROL_DRUGS:
    control_meta.append(meta[meta['drug_type'] == c])
control_meta = pd.concat(control_meta)
control_meta = control_meta[control_meta['worm_strain'] == CONTROL_STRAIN]

control_feat = feat.loc[control_meta.index, prestim_nopath_feats]

# find wells above and below thresholds
upper_thresh = control_feat > (control_feat.mean(axis=0)
                               + BAD_FEAT_THRESH * control_feat.std(axis=0))
lower_thresh = control_feat < (control_feat.mean(axis=0)
                               - BAD_FEAT_THRESH * control_feat.std(axis=0))

flagged_df = upper_thresh | lower_thresh
flagged_df['total_counts'] = flagged_df.sum(axis=1)

control_to_drop = flagged_df[flagged_df['total_counts'] > 120].index

feat = feat.drop(index=control_to_drop)
meta = meta.drop(index=control_to_drop)

# %% normalise z score
featZ = pd.DataFrame(data=stats.zscore(feat,
                                       axis=0),
                     columns=feat.columns,
                     index=feat.index)

# %% pca plot
plot_pca(featZ.values,
         labels=meta['drug_type'].values)


# %% look at coiling
coiled = 'spiroindoline-type'
paralysed = 'lev-1 receptor'

# select only coiled worms
selected_meta = [meta[meta['MOA_specific']==coiled]]
selected_meta.append(meta[meta['MOA_specific'] == paralysed])
selected_meta = pd.concat(selected_meta)

# only keep the two highest concentrations of each drug
drug_list = selected_meta['drug_type'].unique()
drug_dict = {d: selected_meta[selected_meta['drug_type'] == d][
    'imaging_plate_drug_concentration'].unique()
    for d in drug_list
    }
for k, v in drug_dict.items():
    drug_dict[k].sort()

# select only the highest concentrations
selected_max=[]
for c, drug in enumerate(drug_list):
    # print (drug, drug_dict[drug][-1])
    selected_max.append(selected_meta[(
        selected_meta['drug_type'] == drug) & 
        (selected_meta['imaging_plate_drug_concentration'] == drug_dict[drug][-1])
        ])
    if len(drug_dict[drug])>1:
        # print (drug, drug_dict[drug][-2])
        selected_max[c] = selected_max[c].append(
            selected_meta[(selected_meta['drug_type'] == drug) & 
            (selected_meta['imaging_plate_drug_concentration'] == drug_dict[drug][-2])])
selected_max = pd.concat(selected_max)

# only N2s
selected_max = selected_max[selected_max['worm_strain'] == CONTROL_STRAIN]

# data very skewed towards control, so rebalance by slimming down the number of controls
n_controls_to_sample = selected_max.shape[0]

# select 40 controls from the dataset and use this dataset of 80 obs for classifier
control_meta = []
for c in CONTROL_DRUGS:
    control_meta.append(meta[(meta['drug_type'] == c) &
                             (meta['worm_strain'] == CONTROL_STRAIN)
                             ].sample(int(n_controls_to_sample/2)))
control_meta = pd.concat(control_meta)
    
subsample_meta = pd.concat([selected_max,
                            control_meta],
                           axis=0)  
subsample_meta['MOA_specific'].fillna('Negative_control',
                                      inplace=True) 

# filter down the feat mat 
subsample_feat = feat.loc[subsample_meta.index,
                          prestim_nopath_feats]

subsample_featZ = pd.DataFrame(data=stats.zscore(subsample_feat),
                               columns=subsample_feat.columns,
                               index=subsample_feat.index)

# %%

# for every feature determine the mean and std of the Z-score for each class
# (DMSO, No compound and spiroindolines and lev-1)

# find mean, std and cov for each of these classses
feat_stats_meta = list(subsample_meta['MOA_specific'].unique())
feat_stats_mean = []
feat_stats_std = []
for f in feat_stats_meta:
    feat_stats_mean.append(subsample_featZ.loc[
        subsample_meta[subsample_meta['MOA_specific'] == f
                       ].index, prestim_nopath_feats
        ].mean().to_frame().transpose())
    feat_stats_std.append(subsample_featZ.loc[
        subsample_meta[subsample_meta['MOA_specific'] == f
                       ].index, prestim_nopath_feats
        ].std().to_frame().transpose())
    
feat_stats_mean = pd.concat(feat_stats_mean)
feat_stats_mean.reset_index(drop=True, inplace=True)
feat_stats_mean['MOA_specific'] = feat_stats_meta

feat_stats_std = pd.concat(feat_stats_std)
feat_stats_std.reset_index(drop=True, inplace=True)
feat_stats_std['MOA_specific'] = feat_stats_meta

plt.figure()
for f in feat_stats_meta:
    sns.distplot(feat_stats_mean[
        feat_stats_mean['MOA_specific'] == f
        ][prestim_nopath_feats],
        label = f)
plt.legend()
plt.xlabel('Z-score mean')
plt.savefig(SAVETO / '{}_Z-score_mean_per_MOA.png'.format(CONTROL_STRAIN))
    
plt.figure()
for f in feat_stats_meta:
    sns.distplot(feat_stats_std[
        feat_stats_std['MOA_specific'] == f
        ][prestim_nopath_feats],
        label = f)
plt.legend()
plt.xlabel('Z-score stdev')
plt.savefig(SAVETO / '{}_Z-score_std_per_MOA.png'.format(CONTROL_STRAIN))

# remove bad wells that have lots of features that are more than 4 standard
# deviations from the mean

# first just index count up how many features are more than 4 standard deviations
# from the mean for that group
bad_well_indexing = []
for f in feat_stats_meta:
    thresholds = {}
    thresholds['max']= feat_stats_mean[
        feat_stats_mean['MOA_specific']==f][prestim_nopath_feats].values +\
    (4 * feat_stats_std[
        feat_stats_std['MOA_specific']==f][prestim_nopath_feats].values)
    thresholds['min'] = feat_stats_mean[
        feat_stats_mean['MOA_specific']==f][prestim_nopath_feats].values -\
    (4 * feat_stats_std[
        feat_stats_std['MOA_specific']==f][prestim_nopath_feats].values)
    
    upper = subsample_feat.loc[
            subsample_meta[subsample_meta['MOA_specific'] == f].index,
            prestim_nopath_feats] > thresholds['max']
    lower = subsample_feat.loc[
            subsample_meta[subsample_meta['MOA_specific'] == f].index,
            prestim_nopath_feats] < thresholds['min']
    
    both = upper | lower
    bad_well_indexing.append(both)
    
bad_well_indexing = pd.concat(bad_well_indexing)
bad_wells = bad_well_indexing[bad_well_indexing.sum(axis=1) >
                              bad_well_indexing.sum(axis=1).mean() +\
                                  bad_well_indexing.sum(axis=1).std()].index

# drop these wells
subsample_feat.drop(index=bad_wells,
                    inplace=True)
subsample_meta.drop(index=bad_wells,
                    inplace=True)
subsample_featZ = pd.DataFrame(data=stats.zscore(subsample_feat),
                               columns=subsample_feat.columns,
                               index=subsample_feat.index)

# now find the bad features
# cov for each feature
feat_stats_cov = []
for f in feat_stats_meta:
    feat_stats_cov.append((subsample_feat.loc[
            subsample_meta[subsample_meta['MOA_specific'] == f
                           ].index, prestim_nopath_feats
            ].std() / subsample_feat.loc[
            subsample_meta[subsample_meta['MOA_specific'] == f
                           ].index, prestim_nopath_feats
            ].mean()).to_frame().transpose())
                
feat_stats_cov = pd.concat(feat_stats_cov)
feat_stats_cov.reset_index(drop=True, inplace=True)
feat_stats_cov['MOA_specific'] = feat_stats_meta
feat_stats_cov = feat_stats_cov.replace([np.inf, -np.inf], np.nan)

bad_feats = list(feat_stats_cov[prestim_nopath_feats].columns[
    np.abs(feat_stats_cov.mean())>10])
# bad_feats.extend(list(feat_stats_cov[prestim_nopath_feats].columns[
#     feat_stats_cov.mean()<-10]))
prestim_nopath_minuscov_feats = list(set(prestim_nopath_feats) - 
                                     set(bad_feats))

plt.figure()
for f in feat_stats_meta:
    sns.distplot(feat_stats_cov[
        feat_stats_cov['MOA_specific'] == f
        ][prestim_nopath_minuscov_feats],
        label = f,
        bins=100)
plt.legend()


# now do statistical tests 
y_classes = list(zip(subsample_meta.drug_type,
                     subsample_meta.imaging_plate_drug_concentration))
y_classes = np.array(['{}, {}'.format(x[0], x[1]) for x in y_classes])
sig_feats = significant_features.k_significant_feat(subsample_feat[prestim_nopath_minuscov_feats],
                                                    y_classes,
                                                    k=50)


# initialise PCA
postprocessingPCA = PCA()
X2= postprocessingPCA.fit_transform(subsample_featZ.values)
cumvar = np.cumsum(postprocessingPCA.explained_variance_ratio_)
thresh = cumvar <= 0.95 #set 95% variance threshold
cut_off = int(np.argwhere(thresh)[-1])

#make a plot
sns.set_style('whitegrid')
plt.figure()
plt.plot(range(0, len(cumvar)), cumvar*100)
plt.plot([cut_off, cut_off], [0, 100], 'k')
plt.xlabel('Number of Principal Components', fontsize =16)
plt.ylabel('variance explained', fontsize =16)
#plt.savefig(os.path.join(directoryA[:-7], 'Figures', 'agarPCvar.png'), dpi =150)
#plt.savefig(os.path.join(directoryA[:-7], 'Figures', 'agarPCvar.svg'),dpi = 150)

#now put the 1:cut_off PCs into a dataframe
PCname = ['PC_%d' %(p+1) for p in range (0,cut_off+1)]
PC_df = pd.DataFrame(data= X2[:,:cut_off+1],
                     columns = PCname,
                     index=subsample_featZ.index)

PC_plotting = pd.concat([PC_df,
                         subsample_meta],
                         axis=1)

plt.figure()
sns.scatterplot(x='PC_1',
                y='PC_2',
                data=PC_plotting,
                hue='drug_type',
                palette='Paired')
plt.axis('equal')
plt.legend(loc='center left', 
                bbox_to_anchor=(1.0,0.5),
                ncol = 1,
                frameon= True)
plt.tight_layout()
plt.xlabel('PC_1 ({}%)'.format(np.round(cumvar[0]*100,2)))
plt.ylabel('PC_2 ({})%'.format(np.round((cumvar[1]-cumvar[0])*100, 2)))
plt.savefig(SAVETO / '{}_PC1_PC2_spiroindoline_lev_only.png'.format(CONTROL_STRAIN))

# %% try a classifier

# select best features to classify between levamisol and controls

# X = [subsample_featZ.loc[
#     (subsample_meta[subsample_meta['MOA_specific'] == coiled]).index,
#                         prestim_nopath_minuscov_feats
#                         ]]

# X.append(subsample_featZ.loc[
#     (subsample_meta[subsample_meta['drug_type'] == 'DMSO']).index,
#                         prestim_nopath_minuscov_feats
#                         ])

# X = pd.concat(X)
# Xmeta = subsample_meta.loc[X.index, :]
# Xy_classes = Xmeta['MOA_specific']

# split the data for training a classifier
Xtrain = subsample_featZ[prestim_nopath_minuscov_feats].sample(120)
Xmeta_train = subsample_meta.loc[Xtrain.index]
y_train_array = Xmeta_train['MOA_specific'].values
# Xtrain_array = Xtrain.values

Xtest = subsample_featZ[prestim_nopath_minuscov_feats].loc[
    set(subsample_featZ.index) - set(Xtrain.index)]
Xtest_meta = subsample_meta.loc[Xtest.index]

# first do feature selection with LDA
from sklearn.feature_selection import SelectFromModel

# fit LDA
clf = LinearDiscriminantAnalysis(solver='svd').fit(Xtrain,
                                                   Xmeta_train['MOA_specific'].values)

importance = np.abs(clf.coef_)[0]
# print(importance)

# sort importance values for each feature and set threshold according to that
idx_msize = importance.argsort()[-X.shape[0]]
threshold = importance[idx_msize] + 0.01

idx_features = (-importance).argsort()[:X.shape[0]-1]
name_features = np.array(prestim_nopath_minuscov_feats)[idx_features]
print('Selected features: {}'.format(name_features))

sfm = SelectFromModel(estimator=clf,
                      threshold=threshold).fit(Xtrain,
                                               Xmeta_train['MOA_specific'])

X_transform = pd.DataFrame(data=sfm.transform(Xtest),
                           index=Xtest.index)

LDA_plot = pd.concat([X_transform,
                      Xtest_meta], 
                     axis=1)
#%%

# try to classify the data with the selected features from LDA
# use a random forest classifier

from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
clf = RandomForestClassifier(max_depth=2, random_state=0, criterion = 'entropy')
clf.fit(Xtrain[name_features], Xmeta_train['MOA_specific'].values)
print(clf.predict(Xtest[name_features]))
print (Xtest_meta['MOA_specific'])


#%% Train random forest
# =============================================================================
# # now the actual test
# 1. Train LDA and set how many features and required threshold
# 2. Get list of features
# 3. Use this feature list to train random forest classifier
# 4. Use cross validation set to find best classfier 
# 5. Test on final TEST set
# =============================================================================

feat_size_tests = range(51, 2401, 50)
sss = StratifiedShuffleSplit(n_splits=20,
                             train_size=0.6)


scoring = {}
featsets = {}
for f in range(10, 1001, 10):
    scoring[f] = []
    featsets[f] = []
    for train_index, cv_index in sss.split(Xtrain.values,
                                          y_train_array):
        X_tr, X_cv = Xtrain.values[train_index], Xtrain.values[cv_index]
        y_tr, y_cv = y_train_array[train_index], y_train_array[cv_index]
        
        
        clf = LinearDiscriminantAnalysis(solver='svd').fit(X_tr,
                                                           y_tr)

        importance = np.abs(clf.coef_)[0]
        # print(importance)
        
        # sort importance values for each feature and set threshold according to that
        idx_msize = importance.argsort()[f]
        threshold = importance[idx_msize] + 0.01
        
        idx_features = (-importance).argsort()[:(f-1)]
        name_features = np.array(prestim_nopath_minuscov_feats)[idx_features]
        # print(len(idx_features))
        featsets[f].append(name_features)
        
        rfm = RandomForestClassifier(max_depth=None,
                                     random_state=0,
                                     criterion = 'gini',
                                     n_estimators=20)
        rfm.fit(X_tr[:,idx_features],
                y_tr)
        scoring[f].append(rfm.score(X_cv[:, idx_features],
                                 y_cv))


scor = []
for s in scoring.keys():
    scor.append(pd.DataFrame({'score': scoring[s],
                              'feat_number': [s] * len(scoring[s])}))
scor = pd.concat(scor)
scor['selection'] = 'LDA'

ax = sns.lineplot(x='feat_number',
                  y='score',
                  data=scor)

#%%
# compare to random feature set
scoring_rand = {}
featsets_rand = {}

for f in range(10, 1001, 10):
    scoring_rand[f] = []
    featsets_rand[f] = []
    seed = f
    for train_index, cv_index in sss.split(Xtrain.values,
                                          y_train_array):
        X_tr, X_cv = Xtrain.values[train_index], Xtrain.values[cv_index]
        y_tr, y_cv = y_train_array[train_index], y_train_array[cv_index]
        
        seed += f
        
        idx_features = random.sample(range(0,
                                           len(prestim_nopath_minuscov_feats)),f)
        
        name_features = np.array(prestim_nopath_minuscov_feats)[idx_features]
        # print(len(idx_features))
        featsets_rand[f].append(name_features)
        
        rfm = RandomForestClassifier(max_depth=None,
                                     random_state=0,
                                     criterion='gini',
                                     n_estimators=20)
        rfm.fit(X_tr[:,idx_features],
                y_tr)
        scoring_rand[f].append(rfm.score(X_cv[:, idx_features],
                                    y_cv))

scor_rand = []
for s in scoring_rand.keys():
    scor_rand.append(pd.DataFrame({'score': scoring_rand[s],
                                   'feat_number': [s] * len(scoring_rand[s])}))
scor_rand = pd.concat(scor_rand)  
scor_rand['selection'] = 'random'

scor = pd.concat([scor, scor_rand], axis=0)
ax = sns.lineplot(x='feat_number',
                  y='score',
                  data=scor,
                  hue = 'selection')


#%%
# final test in on the held out test set using 60 features
# feat50 = featsets_rand[50][0]
feat50 = sig_feats.copy()
rfm = RandomForestClassifier(max_depth=None,
                                random_state=0,
                                criterion = 'gini',
                                n_estimators=20)
rfm.fit(Xtrain[feat50].values,
        Xmeta_train['MOA_specific'])
rfm.score(Xtest[feat50].values,
          Xtest_meta['MOA_specific'])
