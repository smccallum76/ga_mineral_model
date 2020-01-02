# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 07:35:06 2019

@author: Scott McCallum
"""
import pandas as pd
import numpy as np

def M_matrix(data):
    '''
    Matrix of minerals from XRD
    '''
    # this function organizes the minerals of interest and then normalizes 
    # the minerals to 100%.  The final mineral values are reported in percent
    # matrix of mineral weight percents
    M = data[['Quartz', 'Plag.', 'Calcite', 'Dolomite', 'Chlorite', 
              'Illite_Mica', 'Mx I_S', 'Pyrite']].copy()
    # rename the columns
    M.rename(columns={'Chlorite':'Chl', 'Illite_Mica':'Ill', 'Quartz':'Qtz',
                      'Pyrite':'Pyr', 'Calcite':'Cal', 'Dolomite':'Dol',
                      'Plag.':'Plg', 'Mx I_S':'MxIS'}, inplace=True)
    # reorder cols (alphabetical order)
    M = M[['Cal', 'Chl', 'Dol', 'Ill', 'MxIS', 'Plg', 'Pyr', 'Qtz']].copy()
    
    m_sum = M.sum(axis=1)
    # Normalize the data to 100% by looping through each column and divide by 
    # the sum * 100
    for i in M.columns:    
        M[i] = M[i].div(m_sum)*100
    return M
    
def E_matrix(data):
    '''
    Matrix of elements from XRF
    '''
    # this function organizes the elements of interest and reports the final
    # output in percent
    # put elements in alphabetical order
    E = data[['Al', 'Ca', 'Fe', 'K', 'Mg', 'Na', 'Si']].copy()
    # add unity to E
    E['Unity'] = 100
    return E

def F_variance_table(F):
    # this function builds the empty min max tables called F_min and F_max.  
    # These tables start out with indices as the minerals of interest and the 
    # columns as the elements of interest.  The min/max values are set to zero
    # initially. 
    F_min = F.iloc[:, 0:-1].copy() # just dropping the unity column
    F_max = F.iloc[:, 0:-1].copy() # just dropping the unity column
    for i in F_min.columns:
        # add 'min' and 'max' to the column names
        F_min.rename(columns={i:i+'_min'}, inplace=True)
        F_max.rename(columns={i:i+'_max'}, inplace=True)

    # order the column names alphabetically
    F_min = F_min.reindex(sorted(F_min.columns), axis=1)
    F_max = F_max.reindex(sorted(F_max.columns), axis=1)
    # at this point an all zero matrix has been generated with the correct 
    # column and row names for both the F_min and F_max matrices. 
    return F_min, F_max

def F_min_max_dictionary():
    '''
    Each element has a min and max value that can be obtained.  In order to 
    determine the max value a given element is set to the max number of atoms 
    possible and all other elements are set to the min number of atoms possible.  
    The min value is determined in the exact opposite manner; a given element 
    is set to the min number of atoms possible for each mineral and all other 
    elements are set to the max number of atoms.  It is necessary to work with 
    the number of atoms b/c this will impact the molecular weight, which then 
    impacts the fraction of each element.  Note, b/c dolomite has some Fe 
    dependencies, the min/max values had to be swapped.  This isn't a perfect 
    method, but allowing for some flex in each element allows for accomadation 
    of the fact that not all elements are available for use and neither are all 
    minerals (i.e. there is some unavoidable slop).
    
    If new elements or minerals are added then they can simply be added to 
    the dictionary below. 
    '''
    F_min_max_dict = {
                # min max values
                # CALCITE
                'cal_ca_max' : 0.400479616306954,
                'cal_ca_min' : 0.294236093632016,
                'cal_mg_max' : 0.0765164862824062,
                'cal_mg_min' : 0.0, 
                # CHLORITE
                'chl_al_max' : 0.249809032198329,
                'chl_al_min' : 0.0312286590659182,
                'chl_fe_max' : 0.580279502904112,
                'chl_fe_min' : 0.0,
                'chl_mg_max' : 0.257261937503967,
                'chl_mg_min' : 0.0,
                'chl_si_max' : 0.206924493554328,
                'chl_si_min' : 0.0316076110317201,
                # DOLOMITE
                'dol_ca_max' : 0.217353579175705,
                'dol_ca_min' : 0.185615708794517,
                'dol_mg_max' : 0.131887201735358,
                'dol_mg_min' : 0.0,
                'dol_fe_max' : 0.25864863613208,
                'dol_fe_min' : 0.0,
                # ILLITE
                'ill_al_max' : 0.238759604762864,
                'ill_al_min' : 0.0,
                'ill_ca_max' : 0.00638564313686753,
                'ill_ca_min' : 0.0,
                'ill_fe_max' : 0.476984933423808,
                'ill_fe_min' : 0.0,
                'ill_k_max' :  0.136777052238806,
                'ill_k_min' : 0.0456374088689976,
                'ill_mg_max' : 0.0806639444550909,
                'ill_mg_min' : 0.0,
                'ill_na_max' : 0.00913214259186659,
                'ill_na_min' : 0.0,
                'ill_si_max' : 0.33667935157163,
                'ill_si_min' : 0.0536809068537848,
                # MIXED ILLITE SMECTITE
                'mxis_al_max' : 0.217749588807779,
                'mxis_al_min' : 0.0,
                'mxis_ca_max' : 0.0410231215647742,
                'mxis_ca_min' : 0.0,
                'mxis_fe_max' : 0.436084186933133,
                'mxis_fe_min' : 0.0,
                'mxis_k_max' :  0.111266042514442,
                'mxis_k_min' : 0.0,
                'mxis_mg_max' : 0.0722454920832962,
                'mxis_mg_min' : 0.0,
                'mxis_na_max' : 0.0467814005729913,
                'mxis_na_min' : 0.0,
                'mxis_si_max' : 0.366758062410236,
                'mxis_si_min' : 0.144126957644924,
                # PLAGIOCLASE
                'plg_al_max' : 0.193947236000288,
                'plg_al_min' : 0.102121152477706,
                'plg_ca_max' : 0.159528737462188,
                'plg_ca_min' : 0.0,
                'plg_k_max' : 0.00772132151108829,
                'plg_k_min' : 0.0,
                'plg_na_max' : 0.0981887756191517,
                'plg_na_min' : 0.0,
                'plg_si_max' : 0.301686177639351,
                'plg_si_min' : 0.213542339767224,                                
                # PYRITE
                'pyr_fe_max' : 0.465486489640113,
                'pyr_fe_min' : 0.465486489640113,
                # QUARTZ
                'qtz_si_min' : 0.467465468463971,
                'qtz_si_max' : 0.467465468463971
                } 

    return F_min_max_dict

def F_var_populate(F_min, F_max, F_min_max_dict):
    '''
    The purpose of this function is to populate the F_min and F_max matrices
    with correct values that fall within the prescribed min and max values
    from the dictionary (F_min_max_dict).
    
    F_min = minimum elemental fractions allowable
    F_max = maximum elemental fraction allowable          
    '''
    # F_min builder  
    for r in F_min.index:
        r_lwr = r.lower() # put everything in lower case
        for c in F_min.columns: # loop through the columns
            c_lwr = c.lower() # lower case
            lookup = r_lwr + '_' + c_lwr # lookup value (row and column)
            if lookup in F_min_max_dict:
                # add the lookup value to the F_min matrix.
                F_min.loc[r,c] = F_min_max_dict[lookup]
    # F_max builder            
    for r in F_max.index:
        r_lwr = r.lower() # lower case
        for c in F_max.columns: # loop through the columns
            c_lwr = c.lower() # lower case
            lookup = r_lwr + '_' + c_lwr # lookup value (row and column)
            if lookup in F_min_max_dict:
                # add the lookup value to the F_max matrix. 
                F_max.loc[r,c] = F_min_max_dict[lookup]
    
    return F_min, F_max

def M_matrix_pred(E,F):
    ''' Square matrix '''
    # E: elements as coming from lithoscanner or XRF in parts per hundred (percent)
    # M: minerals as coming from elan or XRD in percent
    # F: tuned elemetal fractions as coming from XRF to XRD or Lithoscanner to ELAN
    # Finv: Inverted F matrix
    
    # E = M.F
    # E.Finv = M or M = E.Finv
    
    # add the unit column
    F['Unity'] = 1
    F = np.array(F) # convert to np array
    # inverse of F
    Finv = np.array(np.linalg.inv(F))
    # solve for M part
    M = E.dot(Finv)
    # Add some column names
    M.columns = ['Cal', 'Chl', 'Dol', 'Ill', 'MxIS', 'Plg', 'Pyr', 'Qtz']
    return M


def parents(pop_size, F, F_max, F_min):
    '''
    The parent population will be user defined.  The output is a list 
    of vectors (F_vecs) where each line of the list is a parent vector (a parent)
    These parent vector are just a line of the fractions (F) values that 
    will be fed to the matrix prediction function in order to calculate 
    fitness.  For example, an 7x7 array would be flattened into a 1x49 vector.
    '''
    F_vecs = [] # empty list of parent vectors
    for i in range(pop_size):
        # loop through each member of the population
        # F_parent is initially just a copy of all but the unity column from 
        # the F matrix, which is a zero matrix.  
        F_parent = F.iloc[:,0:-1].copy() # copy F and drop the unity column
        for r in F_parent.index:
            # loop through each row of the parent
            for c in F_parent.columns:
                # loop through each column of the parent and assign a new value that
                # is randomly selected between the min and max F tables
                high = F_max.loc[r, c+'_max']
                low = F_min.loc[r, c+'_min']
                rando = np.random.uniform(low=low, high=high)       
                F_parent.loc[r,c] = rando
        # After a single parent has been defined flatten the matrix to make
        # it a vector
        vector = list(np.array(F_parent).flatten())
        # append this vector to F_vecs list.
        F_vecs.append(vector)
    return F_vecs


def fitness(pop_size, num_parents, F_vecs, F_max, M_actual, E,metric='rss'):   
    fit_list=[] # empty list to store fitness values
    score_list=[]
    for i in range(pop_size):
        # note, F_max was used to defne shape, but F_min could have been used, or
        # F.iloc[:, 0:-1] could have been used. 
        F_ind = np.reshape(F_vecs[i], np.shape(F_max))
        # convert the F_ind to a df
        F_ind = pd.DataFrame(F_ind)
        # calculate the M predictions based on some matrix of F values
        M_pred = M_matrix_pred(E, F_ind)
        # calculate the fitness for the predicted mineral abundances
        fitness = np.sum(np.square(M_actual - M_pred).sum(axis=1))
        fit_list.append(fitness)
        # score calculation uses root relative square error (rrse), which is the rss
        # divided by the simple model (mean - actual)^2.  In a perfect model
        # rrse would be zero, but it is possible that your model can be worse
        # than using the simple mean, and therefore, greater than 1.  Since
        # the perfect model has rrse=0, the score = 1 - rrse^0.5
        rrse = np.square(M_actual - M_pred).sum(axis=0)
        var = np.mean(M_actual, axis=0)
        var = np.square(M_actual - var).sum(axis=0)
        # the score for each mineral (if 8 minerals, then 8 scores)
        score = 1-(rrse/var)**0.5
        # for grid search, a single score is needed, therefore the avg of all
        # scores is used. The higher this value is, the better.  
        score_rollup = np.mean(score)
        score_list.append(score_rollup)
        
    F_wFitness = pd.DataFrame(F_vecs)
    F_wFitness['Fit'] = fit_list
    F_wFitness['Score'] = score_list
    # sort in the correct way depending on rss or rrse
    if metric == 'rss': # root square error
        F_wFitness.sort_values(by=['Fit'], inplace=True)
        F_wFitness.reset_index(inplace=True, drop=True)
    elif metric == 'rrse': # root relative square error
        F_wFitness.sort_values(by=['Score'], inplace=True, ascending=False)
        F_wFitness.reset_index(inplace=True, drop=True)
    # copy of the upper half of most fit parents that will be used to define 
    # the other half of offspring
    F_parents = F_wFitness.iloc[0:num_parents,0:-2]
    return F_wFitness, F_parents

def crossover(F_parents, num_parents):
    '''
    CROSSOVER
    '''
    next_gen = np.array(F_parents) # offspring will be appended to this
    # initial array of zeros for the offspring
    # offspring will be a matrix of F values where each row is a single 
    # individual
    offspring = np.zeros(shape=np.shape(F_parents))
    for k in range(num_parents):
        crossover = int(np.random.uniform(low=1,high=np.shape(F_parents)[1]))
        # Index of the first parent to mate.
        # code below sourced from, Ahmed Gad:
        # https://towardsdatascience.com/genetic-algorithm-implementation-in-python-5ab67bb124a6
        parent1_idx = k%F_parents.shape[0]
        # Index of the second parent to mate.
        parent2_idx = (k+1)%F_parents.shape[0]
        # The new offspring will have its first half of its genes taken from the 
        # first parent.
        offspring[k, 0:crossover] = F_parents.iloc[parent1_idx, 0:crossover]
        # The new offspring will have its second half of its genes taken from 
        # the second parent.
        offspring[k, crossover:] = F_parents.iloc[parent2_idx, crossover:]
    # next generation of offpsring and parents    
    next_gen = np.append(next_gen, offspring, axis=0)
    return next_gen

def mutation(mutation_rate, F, F_max, F_min, num_parents, pop_size, next_gen):
    '''
    MUTATION
    '''
    # mutation for some randomly selected fraction of genes
    if mutation_rate > 0: # this just allows the user to test 0 mutation
        for j in range(mutation_rate):
            # row index from the F matrix (F_max and F_min)
            mut_row = int(np.random.uniform(low=0, high=np.shape(F.iloc[0,0:-1])[0]-1))
            # col index from the F matrix (F_max and F_min)
            mut_col = int(np.random.uniform(low=0, high=np.shape(F)[1]-1))
            # the highest possible value based on the F_max data
            high = F_max.iloc[mut_row, mut_col]
            # the lowest possible value based on the F_min data
            low = F_min.iloc[mut_row, mut_col]
            # the actual mutation value base on some random number between the 
            # low and high (note, that if the value is between two numbers that are
            # the same then the mut_value will be that value)
            mut_value = np.random.uniform(low=low, high=high)
            # converted index in the flattened array (next gen)
            # this step was confusing as I had to translate a value from a 8x7 array
            # to a vector that was 56 long.  You can't just multiply the row by
            # the col to get the location b/c this blows up on 0x0.  Also, the 
            # discrepancy between zero indexing and the shape of F required the
            # subtraction of 1 from the mut_row and mut_col (above)
            next_gen_mut_col = int(mut_row * np.shape(F.iloc[0,0:-1])[0] + mut_col)
            # num_parents+1 was used below b/c I did not want the mutation to 
            # occur in the parent population, just the offspring
            next_gen_mut_row = int(np.random.uniform(low=num_parents+1, 
                                                     high=pop_size))
            next_gen[next_gen_mut_row, next_gen_mut_col] = mut_value
            F_vecs = next_gen
    else:
        F_vecs = next_gen
        
    return F_vecs
