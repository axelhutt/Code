#!/usr/bin/env python3
import numpy as np
import statistics as stat
import scipy.stats

import classes

def read_in():
    #label_func = 'functional'
    label_func = 'non_functional'
    
    path = './results_AIS/'+label_func+'_electrodes/'
    #path = './results_AIS/functional_electrodes/'
    
    file_list = []
    file_list.append('20110810B') 
    file_list.append('20110817A') 
    file_list.append('20110817B') 
    file_list.append('20110824A')
    file_list.append('20110824B') 
    file_list.append('20110913A') 
    file_list.append('20110913B') 
    
    num_animals = len(file_list)
    num_cond = 3
    num_trials = 20
    
    data = []
    for k in range(num_animals):
        data.append([])
        
        filelabel = path+'AIS_results_'+file_list[k]+'_low_pre.dat'
        data0 = np.loadtxt(filelabel)
        num_chan = (np.shape(data0)[1]-1)//2
        data[k].append([])
        for l in range(num_chan):
            data[k][0].append(data0[:,1+l*2])
        del data0
        
        filelabel = path+'AIS_results_'+file_list[k]+'_low_post.dat'
        data0 = np.loadtxt(filelabel)
        num_chan = (np.shape(data0)[1]-1)//2
        data[k].append([])
        for l in range(num_chan):
            data[k][1].append(data0[:,1+l*2])
        del data0
        
        filelabel = path+'AIS_results_'+file_list[k]+'_high_post.dat'
        data0 = np.loadtxt(filelabel)
        num_chan = (np.shape(data0)[1]-1)//2
        data[k].append([])
        for l in range(num_chan):
            data[k][2].append(data0[:,1+l*2])
        del data0
    
    return data,file_list
D = classes.data()
data,file_list = read_in()
D.compute_two_way_anova_allanimals_AIS(data,file_list)