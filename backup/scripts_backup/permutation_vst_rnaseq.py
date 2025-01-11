import pandas as pd
from itertools import combinations
import math
import numpy as np
import statistics
from itertools import permutations
from bioinfokit import analys, visuz
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import uniform, randint
from matplotlib_venn import venn2
import copy
import random
from vst_function import *
from stattools.resampling import PermutationTest

#Africa
africa = pd.read_csv('../data/rawcounts_africa.csv').drop(columns=['Name', 'Description'])
africa = africa.T
africa['region'] = 'Africa'

#America
america = pd.read_csv('../data/rawcounts_america.csv').drop(columns=['Name', 'Description'])
america = america.T
america['region'] = 'America'

#Central Asia
centralasia = pd.read_csv('../data/rawcounts_centralasia.csv').drop(columns=['Name', 'Description'])
centralasia = centralasia.T
centralasia['region'] = 'CentralAsia'

#East Asia
eastasia = pd.read_csv('../data/rawcounts_eastasia.csv').drop(columns=['Name', 'Description'])
eastasia = eastasia.T
eastasia['region'] = 'EastAsia'

#SouthAsia
southasia = pd.read_csv('../data/rawcounts_southasia.csv').drop(columns=['Name', 'Description'])
southasia = southasia.T
southasia['region'] = 'SouthAsia'

#Gene anotation
genes = pd.read_csv('../data/rawcounts_africa.csv')
genes = genes.iloc[:,0:2]

## Merging raw counts
regions = pd.concat([africa,america,centralasia, eastasia,southasia])
#regions = regions.set_index('region')
expr_df = regions
anotation = regions['region']
expr_df = expr_df.drop(columns =['region'])
expr_df = expr_df.T
print(expr_df.shape)

### 

## Filter out non-expressed genes

expr_df = expr_df.loc[expr_df.sum(axis=1) > 0, :]
print(expr_df.shape)

## Filter out lowly expressed genes
mask_low_vals = (expr_df > 0.3).sum(axis=1) > 2
expr_df = expr_df.loc[mask_low_vals, :]

print(expr_df.shape)


expr_df = expr_df.T
expr_df.insert(0, 'region', anotation)
expr_df = pd.DataFrame(expr_df)
expr_df

####### Groupby regions #######
dt_group = expr_df.groupby('region')

dt_groupped = []
regions = expr_df['region'].unique()
for i in regions:
    dt_groupped.append(dt_group.get_group(str(i)))
    

#### Names to pair-population
regions = expr_df['region'].unique()
combination_names = []
for i in list(combinations(regions,2)):
    combination_names.append(i[0]+str('-')+i[1])
    


## VST for pair-population
vst_dt = []

combination_regions = list(combinations([0,1,2,3,4],2))
for pair in combination_regions:
        statistic = vst((dt_groupped[pair[0]]), (dt_groupped[pair[1]]))
        vst_dt.append(statistic)

vst_dt = pd.DataFrame(vst_dt).set_axis(combination_names)
vst_dt = vst_dt.T
vst_dt 



gene_names = genes.iloc[vst_dt.index]
vst_dt_annotated = vst_dt.set_index(gene_names['Description'])
vst_dt_annotated

p_value_permutation = []

combination_regions = list(combinations([0,1,2,3,4],2))
for region in combination_regions:
    
    p_value= []
    p_value_permutation.append(p_value)
    for i in vst_dt.index:
        permutation = PermutationTest(dt_groupped[region[0]][i], dt_groupped[region[1]][i], stat=vst, n_perm=9999)#mean_gt
        p_value.append(permutation.p_value())

permutation_df = pd.DataFrame(p_value_permutation)
permutation_df = permutation_df.set_axis(combination_names)
permutation_df = permutation_df.T
permutation_df = permutation_df.set_axis(vst_dt_annotated.index)
permutation_df.to_csv('permutation_rnaseq_vst.csv')

