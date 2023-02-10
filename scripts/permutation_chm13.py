import itertools
from stattools.resampling import PermutationTest
import pandas as pd
import numpy as np
import seaborn
import matplotlib
import statistics
from sklearn.decomposition import PCA
from sklearn.decomposition import KernelPCA
from sklearn.manifold import TSNE
from sklearn.manifold import MDS
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
import math
import glob
from itertools import combinations
from matplotlib import pyplot as plt
from Bio import Phylo
import biotite
from pca_plot import *
from stattools.resampling import PermutationTest
from matplotlib_venn import venn2
from mlxtend.evaluate import permutation_test
import vst_function


def permutation_chm13(lenght_cnvs):
    ## Read and merge files 
    anotation = pd.read_csv('../data/SGDP_anotation.csv', sep=',', encoding='latin-1')

    ids_hg19 = pd.DataFrame()
    ids_hg19['SAMPLE'] = anotation['3-Illumina_ID']
    ids_hg19['REGION'] = anotation['10-Region']
    ids_hg19['COUNTRY'] = anotation['11-Country']
    
    chm13 = pd.read_csv('/Users/luciabazan/Downloads/chm13_gene_regions.csv', index_col=0)
    chm13['LENGHT'] = chm13['END'] - chm13['START']
    chm13 = chm13.drop(columns=['START_GENE', 'END_GENE'])
    chm13 = chm13.drop_duplicates()
    
    # Pivot all cnvs in just one dataframe 
    cnv = chm13.pivot_table(index=["CHR", "START", "END"], 
                    columns='SAMPLE', 
                    values='SCORE').reset_index()
    
    final = cnv.iloc[:,3:284]
    counts = final
    counts = counts.fillna(2)
    counts = counts.T
    cnvs = counts.sort_index()
    
    features = ids_hg19.set_index('SAMPLE').merge(cnvs, left_index=True, right_index=True)
    features = features.loc[:,['REGION', 'COUNTRY']]
    features = features.sort_index()
    
    dt = features.merge(cnvs, left_index=True, right_index=True)
    
    ####### Groupby regions #######
    dt_group = dt.groupby('REGION')

    dt_groupped = []
    regions = features['REGION'].unique()
    for i in regions:
        dt_groupped.append(dt_group.get_group(str(i)))
    
    #### Names to pair-population
    regions = dt['REGION'].unique()
    combination_names = []
    for i in list(combinations(regions,2)):
        combination_names.append(i[0]+str('-')+i[1])
    

   ### Permutation

    p_value_permutation = []
    combination_regions = list(combinations([0,1,2,3,4,5,6],2))
    for region in combination_regions:
        p_value= []
        p_value_permutation.append(p_value)
        
        for i in lenght_cnvs:
            permutation = PermutationTest(dt_groupped[region[0]][i], dt_groupped[region[1]][i], stat=vst_function.vst, n_perm=9999)#mean_gt
            p_value.append(permutation.p_value())
    permutation_vst_chm13_deletions_gene_regions = pd.DataFrame(p_value_permutation).set_axis(combination_names)

    return permutation_vst_chm13_deletions_gene_regions