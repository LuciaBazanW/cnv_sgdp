import numpy as np
import pandas as pd
import stattools.resampling as st 
import itertools


## Function for Vst 
def vst_function(
    x=None, 
    y=None):

    vx = np.var(x, axis=0)
    vy = np.var(y, axis=0)

    nx = len(x)
    ny = len(y)

## Vt 
    a = np.concatenate((x,y))
    vt = np.ndarray.var(a,axis=0)

### Vs = (V1*n1+V2*n2)/(n1+n2) 

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
    

    
    
