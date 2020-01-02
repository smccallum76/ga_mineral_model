# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 07:46:18 2019

@author: Scott McCallum
"""

import pandas as pd
import numpy as np
import ga_functions_clean as ga

'''
Tuning workflow for optimizing elemental fractions that will be used to 
predict mineralogy.  

This script requires a square matrix, meaning it needs 1 less element than the
number of minerals that are being predicted.  There can be one less element b/c
the unity makes up for the loss (unity adds one more equation to the system).

If a non-square matrix is used the script must be changed as follows:
    - M_matrix function must be modified
    - E_matrix function must be modified
    - M_matrix_pred must be modified to accomadate non-square inverse mult
'''

''' ------------------------------------------------------------------ '''
''' GA SETTINGS (user defined) '''
''' ------------------------------------------------------------------ '''
pop_size = 100 # number of individuals in population
num_parents = int(pop_size/2)
generations = 1000
mutation_rate = 100 # of rows (chromosomes) that will have one cell (gene) mutated

''' ------------------------------------------------------------------ '''
''' DATA IMPORT '''
''' ------------------------------------------------------------------ '''
litho = pd.read_excel('user elemental and mineral data')
litho.reset_index(inplace=True, drop=True)
depth = litho['Sample'].copy() # USER CAN REMOVE, OR MODIFY FOR DEPTH DATA

''' ------------------------------------------------------------------ '''
''' Configure M, E, and F matrices '''
''' ------------------------------------------------------------------ '''
# The E and M matrix functions should be reviewed prior to running to ensure
# that the proper elements and minerals are being included.  
# Expected minerals = 'Cal', 'Chl', 'Dol', 'Ill', 'MxIS', 'Plg', 'Pyr', 'Qtz'
# Expected elements = 'Al', 'Ca', 'Fe', 'K', 'Mg', 'Na', 'Si'

# E: elements as coming from lithoscanner or XRF in parts per hundred (percent)
# M: minerals as coming from lithoscanner or XRD in percent
# F: tuned elemetal fractions as coming from XRF to XRD or Lithoscanner to ELAN
# E = M.F (linear alg equation)
M_actual = ga.M_matrix(litho)
E = ga.E_matrix(litho)
F = pd.DataFrame(0,index=M_actual.columns, columns=E.columns)
F['Unity'] = 1

''' ------------------------------------------------------------------ '''
''' SETTING MIN/MAX VALUES '''
''' ------------------------------------------------------------------ '''
# A table of ranges for the F matrix need to be built.  The ranges are based
# on a dictionary of values that provides a min and max for each element
# F_max = matrix of max allowable F values
# F_min = matrix of min allowable F values
F_min, F_max = ga.F_variance_table(F)
# define the min max elemental fractions for each mineral of interest
F_min_max_dict = ga.F_min_max_dictionary()
# use the F_min_max_dict to populate the empty F_var matrix
F_min, F_max = ga.F_var_populate(F_min, F_max, F_min_max_dict)

''' ------------------------------------------------------------------ '''
''' INITIAL PARENT POPULATION '''
''' ------------------------------------------------------------------ '''
# make a parent matrix  and convert each to a vector
F_vecs = ga.parents(pop_size, F, F_max, F_min)

''' ------------------------------------------------------------------ '''
''' LOOP THROUGH GENERATIONS '''
''' ------------------------------------------------------------------ '''

best_fit = []
F_vecs = ga.parents(pop_size, F, F_max, F_min) # kill this when plots are finished
for i in range(generations):
    ''' FITNESS '''
    F_wFitness, F_parents = ga.fitness(pop_size, num_parents, 
                                         F_vecs, F_max, M_actual, E)
    ''' FITNESS TRACKING (and reporting) '''
    # calc the min fitness...this needs to be preserved each time it beats 
    # a previous generation (haven't done this yet)
    x = F_wFitness['Fit'].min()  
    best_fit.append(F_wFitness.iloc[0,:])
    print('generation:', i, '  error:', x, ':', F_wFitness.loc[0,'Score'])       
    
    ''' CROSSOVER '''
    next_gen = ga.crossover(F_parents, num_parents)
    ''' MUTATION '''
    F_vecs = ga.mutation(mutation_rate, F, F_max, F_min, num_parents, 
                         pop_size, next_gen)
best_fit = pd.DataFrame(best_fit)


''' ------------------------------------------------------------------ '''
''' OPTIMIZED MODEL '''
''' ------------------------------------------------------------------ '''
#best_fit.sort_values(by=['Fit'], inplace=True)  
best_fit.sort_values(by=['Score'], inplace=True, ascending=False)
best_fit.reset_index(inplace=True, drop=True)
best_model = best_fit.iloc[0,0:-2]

shape = np.shape(F.iloc[:,0:-1])
best_model = best_model.values.reshape(shape)
best_model = pd.DataFrame(best_model)
best_model['Unity'] = 1
best_model.columns = F.columns
best_model.index = F.index