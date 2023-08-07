#!/usr/bin/env python3
import numpy as np
import statistics as stat
import scipy.stats

import classes

def read_in():
    #path = './results_AIS/non_functional_electrodes/'
    path = './results_AIS/functional_electrodes/'
    
    #file_list.append('20110810B') +
    #file_list.append('20110817A') +
    #file_list.append('20110817B') +
    #file_list.append('20110824A')
    #file_list.append('20110824B') +
    #file_list.append('20110913A') +
    #file_list.append('20110913B') +
    
    
    #label11 = path+'AIS_results_20110817B_low_pre.dat'
    #label12 = path+'AIS_results_20110817B_low_post.dat'
    #label13 = path+'AIS_results_20110817B_high_post.dat'
    #label21 = path+'AIS_results_20110913B_low_pre.dat'
    #label22 = path+'AIS_results_20110913B_low_post.dat'
    #label23 = path+'AIS_results_20110913B_high_post.dat'
    
    #label11 = path+'AIS_results_20110824B_low_pre.dat'
    #label12 = path+'AIS_results_20110824B_low_post.dat'
    #label13 = path+'AIS_results_20110824B_high_post.dat'
    #label21 = path+'AIS_results_20110913A_low_pre.dat'
    #label22 = path+'AIS_results_20110913A_low_post.dat'
    #label23 = path+'AIS_results_20110913A_high_post.dat'
    
    #label11 = path+'AIS_results_20110810B_low_pre.dat'
    #label12 = path+'AIS_results_20110810B_low_post.dat'
    #label13 = path+'AIS_results_20110810B_high_post.dat'
    #label21 = path+'AIS_results_20110817A_low_pre.dat'
    #label22 = path+'AIS_results_20110817A_low_post.dat'
    #label23 = path+'AIS_results_20110817A_high_post.dat'
    
    label11 = path+'AIS_results_20110824A_low_pre.dat'
    label12 = path+'AIS_results_20110824A_low_post.dat'
    label13 = path+'AIS_results_20110824A_high_post.dat'
    label21 = path+'AIS_results_20110824A_low_pre.dat'
    label22 = path+'AIS_results_20110824A_low_post.dat'
    label23 = path+'AIS_results_20110824A_high_post.dat'
    
    
    data_11 = np.loadtxt(label11)
    data_12 = np.loadtxt(label12)
    data_13 = np.loadtxt(label13)
    data_21 = np.loadtxt(label21)
    data_22 = np.loadtxt(label22)
    data_23 = np.loadtxt(label23)
    num_AIS1 = int((np.shape(data_11)[1]-1)//2)
    num_AIS2 = int((np.shape(data_21)[1]-1)//2)
    AIS1 = np.zeros((3,num_AIS1))
    AIS2 = np.zeros((3,num_AIS2))
    pvalue1 = np.zeros((3,num_AIS1))
    pvalue2 = np.zeros((3,num_AIS2))
    
    for k in range(num_AIS1):
        AIS1[0,k] = stat.median(data_11[:,1+k*2])
        AIS1[1,k] = stat.median(data_12[:,1+k*2])
        AIS1[2,k] = stat.median(data_13[:,1+k*2])
        pvalue1[0,k] = (scipy.stats.mannwhitneyu(data_11[:,1+k*2],data_12[:,1+k*2]))[1]
        pvalue1[1,k] = (scipy.stats.mannwhitneyu(data_11[:,1+k*2],data_13[:,1+k*2]))[1]
        pvalue1[2,k] = (scipy.stats.mannwhitneyu(data_12[:,1+k*2],data_13[:,1+k*2]))[1]
        print("data 1    #%d   1-2:%.3f(%.4f) 1-3:%.3f(%.4f)  2-3:%.3f(%.4f)"%(k,\
                AIS1[0,k]-AIS1[1,k],pvalue1[0,k],\
                AIS1[0,k]-AIS1[2,k],pvalue1[1,k],\
                AIS1[1,k]-AIS1[2,k],pvalue1[2,k]))
    for k in range(num_AIS2):
        AIS2[0,k] = stat.median(data_21[:,1+k*2])
        AIS2[1,k] = stat.median(data_22[:,1+k*2])
        AIS2[2,k] = stat.median(data_23[:,1+k*2])
        pvalue2[0,k] = (scipy.stats.mannwhitneyu(data_21[:,1+k*2],data_22[:,1+k*2]))[1]
        pvalue2[1,k] = (scipy.stats.mannwhitneyu(data_21[:,1+k*2],data_23[:,1+k*2]))[1]
        pvalue2[2,k] = (scipy.stats.mannwhitneyu(data_22[:,1+k*2],data_23[:,1+k*2]))[1]
        print("data 2    #%d   1-2:%.3f(%.4f) 1-3:%.3f(%.4f)  2-3:%.3f(%.4f)"%(k,\
                AIS2[0,k]-AIS2[1,k],pvalue2[0,k],\
                AIS2[0,k]-AIS2[2,k],pvalue2[1,k],\
                AIS2[1,k]-AIS2[2,k],pvalue2[2,k]))
    return AIS1,AIS2,pvalue1,pvalue2    
        
        
D = classes.data()
AIS1,AIS2,p1,p2 = read_in()
D.plot_barplots(AIS1,AIS2,p1,p2)