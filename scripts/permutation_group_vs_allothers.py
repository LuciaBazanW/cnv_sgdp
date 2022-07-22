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

data = pd.read_csv("cnvator_data_sudmant_overlapped.csv")


anotation = pd.read_csv('/branchinecta/jbazanwilliamson/SGDP_anotation.csv', sep=',', encoding='latin-1',  index_col=1)


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
    
    
dt_groupped_no_group = []
for i in regions:
    dt_groupped_no_group.append(dt.drop(dt_group.get_group(i).index))
    
def variance_gt(x, y):
    return np.var(x, axis=0) - np.var(y, axis=0)


#### Permutation ####
p_value_permutation = []

for region in range(7):
    
    p_value= []
    p_value_permutation.append(p_value)
    for i in range(8650):
        permutation = PermutationTest(dt_groupped[region][i], dt_groupped_no_group[region][i], stat=variance_gt, n_perm=9999)
        p_value.append(permutation.p_value())


permutation_df = pd.DataFrame(p_value_permutation)

permutation_df = permutation_df.set_axis(regions)
permutation_df.to_csv("permutation_results_group_vs_allothers.csv")
