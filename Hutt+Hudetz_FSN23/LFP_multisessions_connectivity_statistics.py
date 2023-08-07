#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import classes

#ana = 'low'#'low' #'medium'#'high'#'medium'#
#ana = 'medium'#'low' #'medium'#'high'#'medium'#
#ana = 'high'#'low' #'medium'#'high'#'medium'#
ana_list = ['low','medium','high']
#path_list = ['./analysis/'+ana_list[0]+'/','./analysis/'+ana_list[1]+'/','./analysis/'+ana_list[22]+'/']
#stype = 'pre'
stype_list = ['pre','stim','post']
#file_list_flag = []
file_list = []

file0='20110810B'
file_list.append(file0)
#file_list_flag.append(1)
file0='20110817A'
file_list.append(file0)
#file_list_flag.append(1)
file0='20110817B'
file_list.append(file0)
#file_list_flag.append(1)
file0='20110824A'
file_list.append(file0)
#file_list_flag.append(1)
file0='20110824B'
file_list.append(file0)
#file_list_flag.append(1)
file0 = '20110913A'
file_list.append(file0)
#file_list_flag.append(1)
file0='20110913B'
file_list.append(file0)
#file_list_flag.append(1)

symmetry_list = []
for l in range(2):# for different fmin
    symmetry_list.append([])
    for k in range(3):
        symmetry_list[l].append([])
        #for m in range(7):
        #    symmetry_list[l][k].append([])
            
D = classes.data()
flag_plots=0
flag_inference = 0
flag_twofmin = 0 # 1

fs = 500.0
flags = [flag_plots,flag_inference,flag_twofmin]



if flag_twofmin == 0:
    fmin = 0.3
    fmax = 2.5
    #fmin=30.0
    #fmax=70.0
    #fmin=0.1
    #fmax=1.5
    #fmin=0.5
    #fmax=4.5
    for ana in ana_list:
        path='./analysis/'+ana+'/'
        for file0 in file_list:
            D.compute_connectivity_statistics(path,ana,stype_list,file0,flags,fmin,fmax,fs)
if flag_twofmin == 1:
    fmin = 0.3
    fmax = 2.5
    for ana in ana_list:
        path='./analysis/'+ana+'/'
        for file0 in file_list:
            D.compute_connectivity_statistics(path,ana,stype_list,file0,flags,fmin,fmax,fs)
    fmin=30.0
    fmax=70.0
    for ana in ana_list:
        path='./analysis/'+ana+'/'
        for file0 in file_list:
            D.compute_connectivity_statistics(path,ana,stype_list,file0,flags,fmin,fmax,fs)

flags[0]=1
D.compute_connectivity_statistics_all(ana_list,fmin,fmax,fs,flags)