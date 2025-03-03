import pandas as pd
from itertools import combinations
import math
import numpy as np
import statistics
from itertools import permutations
from bioinfokit import analys, visuz
import matplotlib.pyplot as plt
from stattools.resampling import PermutationTest
import copy
import random
from scipy import stats

data = pd.read_csv("../data/cnvator_data_sudmant_overlapped.csv")


anotation = pd.read_csv('../data/SGDP_anotation.csv', sep=',', encoding='latin-1',  index_col=1)


cnv = data.pivot_table(index=["Chr", "Start", "End"], 
                    columns='Sample_ID', 
                    values='RD').reset_index()


coordinates = cnv.iloc[:,0:2]
cnvs = cnv.iloc[:,3:284]
cnvs = cnvs.T
cnvs = cnvs.fillna(2)

features = anotation.loc[:,["7-Gender","10-Region", "11-Country",]]
dt = features.merge(cnvs, left_index=True, right_index=True)

####### Groupby regions #######
dt_group = dt.groupby('10-Region')

dt_groupped = []
regions = anotation['10-Region'].unique()
for i in regions:
    dt_groupped.append(dt_group.get_group(str(i)))


def variance_gt(x, y):
    return np.var(x, axis=0) - np.var(y, axis=0)

def mean_gt(x, y):
    return np.mean(x, axis=0) - np.mean(y, axis=0)


def vst_gt(x, y):
    ####### Groupby regions #######
    ######## V = within-population variance ######
    vx = np.var(x, axis=0)
    vy = np.var(y, axis=0)
    ########### N =  numbers of individuals sampled from population each cnv ##############    
    nx = len(x)
    
    ny = len(y)
    
    ######## Vt = total variance across all individuals of the pair of populations ########
    #vt = pd.concat([x,y]).var()
    vt = np.concatenate((x,y)).var()
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
    
#### Permutation ####
p_value_permutation = []

combination_regions = list(combinations([0,1,2,3,4,5,6],2))
for region in combination_regions:
    
    p_value= []
    p_value_permutation.append(p_value)
    for i in range(8650):
        permutation = PermutationTest(dt_groupped[region[0]][i], dt_groupped[region[1]][i], stat=vst_gt, n_perm=9999)
        p_value.append(permutation.p_value())

names = []
regions = dt['10-Region'].unique()
for i in list(combinations(regions,2)):
    names.append(i[0]+str('-')+i[1])


permutation_df = pd.DataFrame(p_value_permutation)

permutation_df = permutation_df.set_axis(names)
permutation_df.to_csv("../permutation_results.csv")
