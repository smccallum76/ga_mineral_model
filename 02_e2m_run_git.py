# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:03:07 2019

@author: Scott McCallum
"""
import pandas as pd
import numpy as np
import e2m_functions_clean as e2m
from scipy.stats.stats import pearsonr

''' -------------------------------------------------------------------- '''
''' IMPORT PATHS '''
''' -------------------------------------------------------------------- '''
data_path = '..\\..\\GenAlg\\input'
# tuned matrix
tune_path = '..\\..\\GenAlg\\Tuned_fractions'
''' -------------------------------------------------------------------- '''
''' MINERALS TO BE MODELED '''
''' -------------------------------------------------------------------- '''
# list of minerals that are being modeled
minerals = ['Cal', 'Chl', 'Dol', 'Ill', 'MxIS', 'Plg', 'Pyr', 'Qtz']
''' -------------------------------------------------------------------- '''
''' DATA IMPORT '''
''' -------------------------------------------------------------------- '''
litho = pd.read_excel(data_path + 'user elemental and mineral data') 
# copy a few columns for later use
sample = litho['Sample'].copy() # may be deleted depending on user data
depth = litho['Depth'].copy() # may be deleted depending on user data
well = litho['well'].copy() # may be deleted depending on user data
# tuned F matrix (elemental fractions)
F = pd.read_excel(tune_path + '\\ga_tune_score_v1.xlsx')
''' -------------------------------------------------------------------- '''        
''' MATRIX CONFIGURATION '''
''' -------------------------------------------------------------------- '''
F, M_true, E = e2m.matrix_config(litho, F)
''' -------------------------------------------------------------------- '''
''' MATRIX OPERATIONS '''
M = e2m.M_matrix_sqr(E,F)
''' -------------------------------------------------------------------- '''
''' ADD DEPTH COLUMN '''
''' -------------------------------------------------------------------- '''
M['Depth'] = depth
M_true['Depth'] = depth
''' -------------------------------------------------------------------- '''    
''' Bulk Minerals '''
''' -------------------------------------------------------------------- '''
# lumping of minerals:
#   tectosilicates = quartz + plag + (kspar if modeled)
#   carbonates = calcite + dolomite
#   phyllosilicates + sum of all modeled clays
M['Tecto'] = M['Qtz'] + M['Plg']
M['Carb'] = M['Cal'] + M['Dol']
M['Clay'] = M['Chl'] + M['Ill'] + M['MxIS']

M_true['Tecto'] = M_true['Qtz'] + M_true['Plg']
M_true['Carb'] = M_true['Cal'] + M_true['Dol']
M_true['Clay'] = M_true['Chl'] + M_true['Ill'] + M_true['MxIS']
    
''' -------------------------------------------------------------------- '''
''' ERROR STATISTICS Full Mineralogy'''
''' -------------------------------------------------------------------- '''
# RMSE
y_true = M_true[minerals].copy()
y_pred = M[minerals].copy()
rsqr_minerals=[]
for i in minerals:
    rsqr = pearsonr(y_true.loc[:,i], y_pred.loc[:,i])[0]**2
    rsqr_minerals.append(round(rsqr, 3))
    print(i, '->', round(rsqr, 3), 'r squared')

# mean of each mineral for true and predicted
mean_true = y_true.mean()
mean_pred = y_pred.mean()
# mean abs error
mae_mineral = np.abs(y_true - y_pred)
mae_mineral = np.mean(mae_mineral, axis=0)
# mean square error
mse_mineral = (y_true - y_pred)**2
mse_mineral = np.mean(mse_mineral)
# residual sum of squares
rss_mineral = np.sum((y_true - y_pred)**2, axis=0)
# root relative square error
# if score = 1, then perfect model
# if score = 0, then model is no better than using the mean
# if score < 0, then model performs worse than the average
rrse_mineral = np.square(y_true - y_pred).sum(axis=0)
var_mineral = np.mean(y_true, axis=0)
var_mineral = np.square(y_true - var_mineral).sum(axis=0)
# the score for each mineral (if 8 minerals, then 8 scores)
score_mineral = 1-(rrse_mineral/var_mineral)**0.5
# putting all the metrics into a df for ease of reading
mineral_mets = pd.DataFrame()
mineral_mets['mean abs error'] = round(mae_mineral, 1)
mineral_mets['mean square error'] = round(mse_mineral, 1)
mineral_mets['res sum squares'] = round(rss_mineral, 1)
mineral_mets['score'] = round(score_mineral,3)
# r-square
mineral_mets['r-square'] = rsqr_minerals
# mean of true and predicted
mineral_mets['mean_true'] = mean_true
mineral_mets['mean_predicted'] = mean_pred

print('\n', '-----------------------------------------')
''' -------------------------------------------------------------------- '''
''' ERROR STATISTICS Bulk Mineralogy'''
''' -------------------------------------------------------------------- '''
# RMSE
bulk = ['Tecto', 'Carb', 'Clay']
y_true = M_true[bulk].copy()  
y_pred = M[bulk].copy()
rsqr_bulk=[]
for i in bulk:
    rsqr = pearsonr(y_true.loc[:,i], y_pred.loc[:,i])[0]**2
    rsqr_bulk.append(round(rsqr, 3))
    print(i, '->', round(rsqr, 3), 'r squared')
    
# mean of each mineral for true and predicted
mean_true = y_true.mean()
mean_pred = y_pred.mean()
# mean absolute error
mae_bulk = np.abs(y_true - y_pred)
mae_bulk = np.mean(mae_bulk, axis=0)
# mean square error
mse_bulk = (y_true - y_pred)**2
mse_bulk = np.mean(mse_bulk)
# residual sum of squares
rss_bulk = np.sum((y_true - y_pred)**2, axis=0)
# root relative square error
# if score = 1, then perfect model
# if score = 0, then model is no better than using the mean
# if score < 0, then model performs worse than the average
rrse_bulk = np.square(y_true - y_pred).sum(axis=0)
var_bulk = np.mean(y_true, axis=0)
var_bulk = np.square(y_true - var_bulk).sum(axis=0)
# the score for each mineral (if 8 minerals, then 8 scores)
score_bulk = 1-(rrse_bulk/var_bulk)**0.5
# putting all the metrics into a df for ease of reading
bulk_mets = pd.DataFrame()
bulk_mets['mean abs error'] = round(mae_bulk, 1)
bulk_mets['mean square error'] = round(mse_bulk, 1)
bulk_mets['res sum squares'] = round(rss_bulk, 1)
bulk_mets['score'] = round(score_bulk,3)
# r-square
bulk_mets['r-square'] = rsqr_bulk
# mean of true and predicted
bulk_mets['mean_true'] = mean_true
bulk_mets['mean_predicted'] = mean_pred
''' -------------------------------------------------------------------- '''
''' SAVE FILES '''
''' -------------------------------------------------------------------- '''
M['Well'] = well
M['Sample'] = sample
M['Depth'] = depth

#mineral_mets.to_excel('output\\mineral_metrics_rss.xlsx')
#bulk_mets.to_excel('output\\bulk_metrics_rss.xlsx')
#M.to_excel('output\\mineral_predictions_rss.xlsx', index=False)

#mineral_mets.to_excel('output\\mineral_metrics_score.xlsx')
#bulk_mets.to_excel('output\\bulk_metrics_score.xlsx')
#M.to_excel('output\\mineral_predictions_score.xlsx', index=False)