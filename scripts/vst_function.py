import numpy as np
import pandas as pd
import stattools.resampling as st 
import itertools


def vst(
    x=None, 
    y=None
):
    """
    Vst Statistic for two groups
    x(object): Values from group A
    y(object): Values from group B
    """
    
    ####### Groupby regions #######
    ######## V = within-population variance ######
    vx = np.var(x, axis=0)
    vy = np.var(y, axis=0)
    ########### N =  numbers of individuals sampled from population each cnv ##############    
    #nx = x.drop(columns = ['7-Gender', '10-Region', '11-Country'])
    nx = len(x)
    #ny = y.drop(columns = ['7-Gender', '10-Region', '11-Country'])
    ny = len(y)


    ######## Vt = total variance across all individuals of the pair of populations ########
    if isinstance(x, np.ndarray):
        vt = np.concatenate((x,y)).var()
        
    else:
        vt = pd.concat([x,y]).var()
    ########## Vs ################## 
    ### Vs = (V1*n1+V2*n2)/(n1+n2) 
    ## where V1 is the within-population variance of population 1, 
    ## V2 is the within population variance of population 2, 
    #n1 and n2 are the numbers of individuals sampled from population 1 and 2, respectively.
    v1 = vx*nx
    v2 = vy*ny
    ns = nx+ny
    vs = (v1+v2)/ ns
    ########## Vst #################
    #####(VT−VS)/VT
    
    vst = (vt-vs)/vt
    
        
    return(vst)


def dmean(
    x=None, 
    y=None

):
    """
    Dmean for two groups
    x(object): Values from group A
    y(object): Values from group B
    """
    ####### Groupby regions #######
    ######## M = median within-populatio ######
    mx = np.mean(x, axis=0)
    my = np.mean(y, axis=0)
    
    m = mx - my
    m = abs(m)
    
    return(m)


def permut(lenght_cnvs):
    from itertools import combinations
    from mlxtend.evaluate import permutation_test
    p_value_permutation = []
    combination_regions = list(combinations([0,1,2,3,4,5,6],2))
    for region in combination_regions:
        p_value= []
        p_value_permutation.append(p_value)
        for i in lenght_cnvs:    
            permutation_analysis = permutation_test(dt_groupped[region[0]][i], dt_groupped[region[1]][i], method='approximate',
                           num_rounds=10000,
                           seed=0, 
                          func= vst)#mean_gt
            p_value.append(permutation_analysis.p_value())
    permutation_vst_chm13_deletions_gene_regions = pd.DataFrame(p_value_permutation).set_axis(combination_names)

    return permutation_vst_chm13_deletions_gene_regions
    

    
    
