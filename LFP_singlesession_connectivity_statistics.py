#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import classes

ana = 'high'#'medium'#'medium'#'high'#'medium'#
path = './analysis/'+ana+'/'
#stype = 'pre'
stype_list = ['pre','stim','post']
#file_list_flag = []
file_list = []

#file0='20110810B'
#file_list.append(file0)
#file_list_flag.append(1)
#file0='20110817A'
#file_list.append(file0)
#file_list_flag.append(1)
#file0='20110817B' ## first bad example in medium
#file_list.append(file0)
#file_list_flag.append(1)
#file0='20110824A'
#file_list.append(file0)
#file_list_flag.append(1)
#file0='20110824B'
#file_list.append(file0)
#file_list_flag.append(1)
#file0 = '20110913A'
#file_list.append(file0)
#file_list_flag.append(1)
file0='20110913B' # good example in medium
#file_list.append(file0)
#file_list_flag.append(1)

flag_plots = 1
flag_inference = 0
flags = [flag_plots,flag_inference]
#fmin=30.0
#fmax=70.0
fmin=0.5
fmax=4.5
#fmin = 0.1
#fmax = 1.5
#fmin = 0.3
#fmax = 2.5
fs=500.0

D = classes.data()
D.compute_connectivity_statistics(path,ana,stype_list,file0,flags,fmin,fmax,fs)
