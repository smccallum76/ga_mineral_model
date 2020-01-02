# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:04:21 2019

@author: Scott McCallum
"""
import numpy as np

def matrix_config(litho,F):

    ''' ELEMENTAL FRACTIONS (F) '''
    # import the tuned elemental fractions.  These are the fractions of the key
    # elements (Al, Si, Ca, Mg, Fe, K, Na) in each of the key minerals (chl, 
    # ill, qtz, pyr, cal, dol). 
    
    F = F[['Al', 'Ca', 'Fe', 'K', 'Mg', 'Na', 'Si', 'Unity']].copy()
    ''' MINERALS '''
    # matrix of mineral weight percents (data specific)
    M = litho[['Quartz', 'Plag.', 'Calcite', 'Dolomite', 'Chlorite', 
              'Illite_Mica', 'Mx I_S', 'Pyrite']].copy()
    # rename the columns (data specific)
    M.rename(columns={'Chlorite':'Chl', 'Illite_Mica':'Ill', 'Quartz':'Qtz',
                      'Pyrite':'Pyr', 'Calcite':'Cal', 'Dolomite':'Dol',
                      'Plag.':'Plg', 'Mx I_S':'MxIS'}, inplace=True)
    # reorder cols
    M = M[['Cal', 'Chl', 'Dol', 'Ill', 'MxIS', 'Plg', 'Pyr', 'Qtz']].copy()
    m_sum = M.sum(axis=1)
    # Normalize the data to 100% by looping through each column and divide by 
    # the sum * 100
    for i in M.columns:    
        M[i] = M[i].div(m_sum)*100
    ''' ELEMENTS '''
    E = litho[['Al', 'Ca', 'Fe', 'K', 'Mg', 'Na', 'Si']].copy()
    # add unity to E
    E['Unity'] = 100

    return F, M, E

def M_matrix_sqr(E,F):
    # E: elements as coming from lithoscanner or XRF in parts per hundred (percent)
    # M: minerals as coming from elan or XRD in percent
    # F: tuned elemetal fractions as coming from XRF to XRD or Lithoscanner to ELAN
    # Finv: Inverted F matrix
    
    # E = M.F
    # E.Finv = M or M = E.Finv
    
    # add the unity column
    F['Unity'] = 1
    F = np.array(F) # convert to np array
    # inverse of T
    Finv = np.array(np.linalg.inv(F))

    # solve for M part 1
    M = E.dot(Finv)
    # Add some column names
    M.columns = ['Cal', 'Chl', 'Dol', 'Ill', 'MxIS', 'Plg', 'Pyr', 'Qtz']
    
    ''' CLEAN UP THE NEGS AND NORMALIZE '''
    # kill the negatives
    M[M<0] = 0
    # sum up the mineral weight percents
    m_sum = M.sum(axis=1)
    # Normalize the data to 100% by looping through each column and divide by 
    # the sum * 100
    for i in M.columns:    
        M[i] = M[i].div(m_sum)*100
    return M
  