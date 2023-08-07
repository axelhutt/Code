
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

import classes

parameters = []
parameters.append(['low','pre'])
parameters.append(['low','post'])
parameters.append(['high','post'])

path = './analysis/'
file_='20110913B'

channum = 7#7#12

D = classes.data()
D.extract_set_timeseries(path,parameters,file_,channum) # 0913B, 0824A
    