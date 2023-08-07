#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
 
import classes

ana_list = ['low','medium','high']

#ana = 'medium'
#path = './analysis/'+ana+'/'
#stype = 'pre'
file_list = []
#file_list_flag = []

file_list.append('20110810B')
#file_list_flag.append(1)
file_list.append('20110817A')
#file_list_flag.append(1)
file_list.append('20110817B')
#file_list_flag.append(1)
file_list.append('20110824A')
#file_list_flag.append(1)
file_list.append('20110824B')
#file_list_flag.append(1)
file_list.append('20110913A')
#file_list_flag.append(1)
file_list.append('20110913B')
#file_list_flag.append(1)



ana = ana_list[0] # low
#ana = ana_list[2] # high

#file = file_list[0]
#file = file_list[1]
#file = file_list[2]
file = file_list[3]
#file = file_list[4]
#file = file_list[5]
#file = file_list[6]

stype = 'pre'
#stype = 'post'
fs = 500.0

duration = 1.0 # in seconds

D = classes.data()
D.analyse_InformationTheory(ana,file,stype,fs,duration)
    
        