#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

import classes

ana_list = ['low','high']
#path = './analysis/'+ana+'/'
path = './analysis/'
#stype = 'pre'
file_list = []
stype_list = ['pre','post']
#file_list_flag = []

file_list.append('20110810B')
file_list.append('20110817A')
file_list.append('20110817B')
file_list.append('20110824A')
file_list.append('20110824B')
file_list.append('20110913A')
file_list.append('20110913B')

num_files = len(file_list)

D = classes.data()
power_all_list = []

D.powerbands_all(path,ana_list,stype_list,file_list) # 0913B, 0824A





#for chan in range(D.num_chan):
#   chan0 = chan+1
#   ax = plt.subplot(4,4,chan0)
#   plt.plot(power_all_list[0][0],power_all_list[0][1][chan,:],label='pre')
#   plt.plot(power_all_list[1][0],power_all_list[1][1][chan,:],label='stim')
#   plt.plot(power_all_list[2][0],power_all_list[2][1][chan,:],label='post')
#   ax.legend()
#plt.show()

#D.convert_dat_to_hdf5(dir0,dir1,label,analevel,labeln)
#D.read_data_artred(dir0,dir1,label,analevel)
#D.compute_write_powerspectra()


#file_list_stim = []
#file_list_stim.append('20110810B_low_stim.h5')
#file_list_stim.append('20110817A_low_stim.h5')
#file_list_stim.append('20110817B_low_stim.h5')
#file_list_stim.append('20110824A_low_stim.h5')
#file_list_stim.append('20110824B_low_stim.h5')
#file_list_stim.append('20110913A_low_stim.h5')
#file_list_stim.append('20110913B_low_stim.h5')
#
#file_list_post = []
#file_list_post.append('20110810B_low_post.h5')
#file_list_post.append('20110817A_low_post.h5')
#file_list_post.append('20110817B_low_post.h5')
#file_list_post.append('20110824A_low_post.h5')
#file_list_post.append('20110824B_low_post.h5')
#file_list_post.append('20110913A_low_post.h5')
#file_list_post.append('20110913B_low_post.h5')
#
#
######### low anaesthesia
#analevel = 'low'
#dir0='20110810'
#dir1='B'
#label = '141236' #pre-3.6%
#labeln = analevel+'_pre'

#label = '142258' #stim-3.6%
#labeln = analevel+'_stim'

#label = '143322' #post-3.6%
#labeln = analevel+'_post'
##------------------------


#dir0='20110817'
#dir1='A'

#label = '100125' #pre-3.4%
#labeln = analevel+'_pre'

#label = '101232'  #stim-3.4%
#labeln = analevel+'_stim'

#label = '102258' #post-3.4%
#labeln = analevel+'_post'
##------------------------

#dir0='20110817'
#dir1='B'
#label = '143429' #pre-3.2%
#labeln = analevel+'_pre'

#label = '145103' #stim-3.2%
#labeln = analevel+'_stim'

#label = '150346' #post-3.2%  beta+gamma peak
#labeln = analevel+'_post'
##------------------------


#dir0='20110824'
#dir1='A'
#label = '095341' #pre-3.4%
#labeln = analevel+'_pre'

#label = '100405' #stim-3.4%
#labeln = analevel+'_stim'

#label = '101435' #post-3.4%
#labeln = analevel+'_post'
##------------------------

#dir0='20110824'
#dir1='B'
#label = '142058' #pre-3.8% gamma+beta peak
#labeln = analevel+'_pre'

#label = '143128' #stim-3.8%
#labeln = analevel+'_stim'

#label = '144201' #post-3.8%
#labeln = analevel+'_post'
##------------------------

#dir0='20110913'
#dir='A'
#label = '091038' #pre-3.4%
#labeln = analevel+'_pre'

#label = '092105' #stim-3.4%
#labeln = analevel+'_stim'

#label = '093134'#post-3.4%
#labeln = analevel+'_post'
##------------------------



#dir0='20110913'
#dir1='B'
#label = '124242' #pre-3.6% 
#labeln = analevel+'_pre'

#label = '125305' #stim-3.6%
#labeln = analevel+'_stim'

#label = '130326' #post-3.6%
#labeln = analevel+'_post'
##------------------------


######## medium anaesthesia
#analevel = 'medium'
#dir0='20110810'
#dir1='B'
#label = '130846' #pre-4.5%
#labeln = analevel+'_pre'

#label = '131915' #stim-4.5%
#labeln = analevel+'_stim'

#label = '133032' #post-4.5%
#labeln = analevel+'_post'
##------------------------

#dir0='20110817'
#dir1='A'
#label = '105605' #pre-4.4%
#labeln = analevel+'_pre'

#label = '110638' #stim-4.4% beta+gamma
#labeln = analevel+'_stim'

#label = '111710' #post-4.4% beta+gamma
#labeln = analevel+'_post'
##------------------------


#dir0='20110817'
#dir1='B'
#label = '152321' #pre-4.4%
#labeln = analevel+'_pre'

#label = '153347' #stim-4.4%
#labeln = analevel+'_stim'

#label = '154410' #post-4.4% gamma
#labeln = analevel+'_post'
##------------------------


#dir0='20110824'
#dir1='A'
#label = '105219' #pre-4.6% gamma
#labeln = analevel+'_pre'

#label = '110311' #stim-4.6% gamma
#labeln = analevel+'_stim'

#label = '111336' #post-4.6% gamma
#labeln = analevel+'_post'
##------------------------

#dir0='20110824'
#dir1='B'
#label = '151704' #pre-5.0%
#labeln = analevel+'_pre'

#label = '155128' #stim-5.0% beta+gamma
#labeln = analevel+'_stim'

#label = '160150' #post-5.0%  beta+gamma
#labeln = analevel+'_post'
##------------------------


#dir0='20110913'
#dir1='A'
#label = '095555' #pre-4.6% beta+gamma
#labeln = analevel+'_pre'

#label = '101810' #stim-4.6% gamma
#labeln = analevel+'_stim'

#label = '102834' #post-4.6% gamma
#labeln = analevel+'_post'
##------------------------

#dir0='20110913'
#dir1='B'
#label = '132936' #pre-4.3% gamma
#labeln = analevel+'_pre'

#label = '134008' #stim-4.3%
#labeln = analevel+'_stim'

#label = '135024' #post-4.3%
#labeln = analevel+'_post'
##------------------------


######## high anaesthesia
#analevel = 'high'
#dir0='20110810'
#dir1='B'
#label = '121713' #pre-6.0%
#labeln = analevel+'_pre'

#label = '123022' #stim-6.0%
#labeln = analevel+'_stim'

#label = '124107' #post-6.0%
#labeln = analevel+'_post'
##------------------------



#dir0='20110817'
#dir1='A'
#label = '114515' #pre-6.0%
#labeln = analevel+'_pre'

#label = '115543' #stim-6.0% gamma
#labeln = analevel+'_stim'

#label = '120605' #post-6.0%
#labeln = analevel+'_post'
##------------------------

#dir0='20110817'
#dir1='B'
#label = '133358' #pre-6.0%
#labeln = analevel+'_pre'

#label = '134428' #stim-6.0% beta+gamma
#labeln = analevel+'_stim'

#label = '135458' #post-6.0% beta+gamma
#labeln = analevel+'_post'
##------------------------


#dir0='20110824'
#dir1='A'
#label = '115341' #pre-6.1% gamma
#labeln = analevel+'_pre'

#label = '120418' #stim-6.1% beta+gamma
#labeln = analevel+'_stim'

#label = '121733' #post-6.1% beta+gamma
#labeln = analevel+'_post'
##------------------------


#dir0='20110824'
#dir1='B'
#label = '132740' #pre-6.1%
#labeln = analevel+'_pre'

#label = '133808' #stim-6.1%
#labeln = analevel+'_stim'

#label = '134831' #post-6.1%
#labeln = analevel+'_post'
##------------------------


#dir0='20110913'
#dir1='A'
#label = '105534' #pre-6.0%
#labeln = analevel+'_pre'

#label = '110634' #stim-6.0%
#labeln = analevel+'_stim'

#label = '111651' #post-6.0% gamma
#labeln = analevel+'_post'
##------------------------


#dir0='20110913'
#dir1='B'
#label = '115353' #pre-6.0% gamma
#labeln = analevel+'_pre'

#label = '120413' #stim-6.0%
#labeln = analevel+'_stim'

#label = '121435' #post-6.0%
#labeln = analevel+'_post'
##------------------------
    