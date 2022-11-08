import numpy as np
import pandas as pd


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