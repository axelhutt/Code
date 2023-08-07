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

fs = 500.0

D = classes.data()
D.compute_PE_statistics(ana_list,file_list,fs)
    
        


