import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import scipy.fftpack
import matplotlib.mlab as mlab
import sys

#from spectrum import *
from cwt_modules_Bergner import *
import scipy.stats as stat

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

time_limit = 3000000#12000


print ('Number of arguments:', len(sys.argv), 'arguments.')
if len(sys.argv)>3:
    quit()
print ('Argument List:', str(sys.argv))

## field
#label_field = './out_field.dat'
label_field = str(sys.argv[1])
f_field = open(label_field,'r')
counter = 0
for line in f_field:
    line=line.strip()
    if counter==0:
        columns = line.split()
        num = size(columns)
    counter += 1
    if counter == time_limit:
        break
    
num_times = counter
f_field.close()
print("num_times:",num_times)
data_field = np.zeros((num_times,num))
counter = 0
f_field = open(label_field)
for line in f_field:
    line=line.strip()
    columns = line.split()
    #print(np.shape(columns),counter)
    for k in range(num):
        data_field[counter,k] = float(columns[k])
    counter += 1
    if counter == time_limit:
        break

f_field.close()    
#print(data_field)
#quit()



## spikes
f_spikes = open(str(sys.argv[2]),'r')
#f_spikes = open('./out_spike.dat')
#f_spikes = open('./spike_test.dat')
counter = 0
counter_singletime = 0
counter_times = -1
oldtime = -1
time = -1
spikes_collected = []
spikes_collected_times = []
for line in f_spikes:
    line=line.strip()
    columns = line.split()
    #print
    oldtime = time
    time = int(columns[0])
    if oldtime != time: ## new time
        ## add neuron to list at current time
        #print("time:%d"%time)
        spikes_collected_neurons=[]
        neuron = int(columns[1])
        spikes_collected_neurons.append(neuron)
        
        spikes_collected.append(spikes_collected_neurons)
        spikes_collected_times.append(time)
        
        #print("spikes_collected_neurons:",spikes_collected_neurons)
        #print("spikes_collected:",spikes_collected)
        
        counter_singletime = 0
        counter_times += 1
    if oldtime == time:
        neuron = int(columns[1])
        list_ = spikes_collected[counter_times] 
        #print("time:%d old list:"%time,list_)
        (spikes_collected[counter_times]).append(neuron)
        counter_singletime += 1
        #print("single time:%d new list:"%counter_singletime,spikes_collected[counter_times])
    counter += 1
    #print("time:%d time_limit:%d"%(time,time_limit))
    if time >= time_limit:
        print("break !!!!!!!!!!!!!!!!!!!!!!")
        break
num_spike_times = len(spikes_collected)
print("num_spike_times:%d"%num_spike_times)
#print("final list:",spikes_collected,"   with %d time steps"%len(spikes_collected))
#print("list of times:",spikes_collected_times)
#print("last spike time read:%d"%max(spikes_collected_times))
f_spikes.close()    


dt = 0.0005

#
#
### time-frequency decomposition
print ('\n perform STA analysis....')
#### STA
df = 2.0
num_perseg  = int(1.0/(df*dt)) # time window length in which FFT is computed 
time_window = int(num_perseg*1.0) #2000 # time_window in which the EEG is considered about the respective spike
time_window2 = int(time_window/2)

## count how many spikes can be used for STA
counter_i = 0
spike_times_used = []
for i in range(num_spike_times):
 if spikes_collected_times[i]+time_window2<num_times and spikes_collected_times[i]>time_window2:
    counter_i += 1
    spike_times_used.append(spikes_collected_times[i])
num_spike_times_used = counter_i
max_spiketime = max(spike_times_used)
#print("spike_times_used:",spike_times_used)
print("length of spike_times_used:",len(spike_times_used))
print("last spike time used:%d"%max(spike_times_used))

counter_neurons=np.zeros(num_spike_times_used)

fmin = 1.0
fmax = 60.0
fmin_num=int(fmin/df)
fmax_num=int(fmax/df)

power=np.zeros((num_spike_times_used,int(num_perseg/2)+1)) # num of frequencies = time_window 
power_all=np.zeros((num_spike_times_used,int(num_perseg/2)+1)) # num of frequencies = time_window 
SFC_spikes=np.zeros((int(num_perseg/2)+1,num_spike_times_used)) # num of frequencies = time_window 
SFC0=np.zeros((int(num_perseg/2)+1,max_spiketime+1)) # num of frequencies = time_window 
counter_times=0
fs=1.0/dt
num_overlap = int(num_perseg*0.9)
#print(data_field)
for j in np.arange(num_spike_times_used): # loop over times when a spike occurs
    if j % 200==0:
        print("j=%d (%d)"%(j,num_spike_times_used))
    spikes = spikes_collected[j]
    num_neurons = len(spikes)
    time_init = spike_times_used[j]-time_window2
    range_=np.arange(time_init,time_init+time_window,1)
    
    sta=np.zeros(time_window)
    if num_neurons>0:
        for k in range(num_neurons):
            #neuron_number1_ = np.random.random_integers(0,num_neurons-1,1)
            neuron_number1_ = (0,num_neurons)
            neuron_number1 = spikes[np.squeeze(neuron_number1_)]
            #print (' neuron no. %d ' %(neuron_number1))
            #data_tf1=np.squeeze(data[:,neuron_number1])
            
            data_help=data_field[range_,neuron_number1]
            sta=sta+data_help/float(num_neurons)
            
            [power0,freqs]=mlab.psd(\
                                    data_help,\
                                    NFFT=num_perseg,\
                                    Fs=fs,\
                                    detrend=mlab.detrend_none,\
                                    window=mlab.window_hanning,\
                                    noverlap=num_overlap,\
                                    pad_to=None,\
                                    sides='onesided',\
                                    scale_by_freq=True)
            power_all[j,:]=power_all[j,:]+power0/float(num_neurons)
            
        [power_sta,freqs]=mlab.psd(\
                                sta,\
                                NFFT=num_perseg,\
                                Fs=fs,\
                                detrend=mlab.detrend_none,\
                                window=mlab.window_hanning,\
                                noverlap=num_overlap,\
                                pad_to=None,\
                                sides='onesided',\
                                scale_by_freq=True)
        
        #print(np.shape(SFC0),j,spike_times_used[j])
        SFC_spikes[:,j]=power_sta/power_all[j,:]  
        SFC0[:,spike_times_used[j]]=power_sta/power_all[j,:] 
    
    del sta
    
    #test1 = power_sta[(SFC_spikes[:,j]<0) | (SFC_spikes[:,j]>1)]
    #test2 = power_sta[(SFC_spikes[:,j]<0) | (SFC_spikes[:,j]>1)]
    #num_test = np.size(test1)
    #for ii in range(num_test):
    #print("ii=%d test1=%f test2=%f"%(ii,test1[ii],test2[ii]))
        
    #print("used spike time: %d"%spike_times_used[j])

print("max_spiketime+1=",max_spiketime+1)

## delta_band
label_delta='SFC_delta.dat'
file_delta=open(label_delta,'w+')
f_delta_low=0.0
f_delta_high=4.0
f_delta_low_num=int(f_delta_low/df)
f_delta_high_num=int(f_delta_high/df)
f_delta_diff_num = f_delta_high_num-f_delta_low_num
SFC_delta = np.zeros(max_spiketime+1)
SFC_delta_total = 0.0
SFC_delta_stdv = 0.0
for i in range(max_spiketime+1):
    for k in range(f_delta_diff_num):
        SFC_delta[i]+=SFC0[k+f_delta_low_num,i]/float(f_delta_diff_num)
    SFC_delta_total+=SFC_delta[i]/float(max_spiketime+1)
    SFC_delta_stdv+=SFC_delta[i]**2/float(max_spiketime+1)
    str = '%f \n'%SFC_delta[i]
    file_delta.write(str)
file_delta.close()
print("written ",label_delta)
SFC_delta_stdv = np.sqrt(SFC_delta_stdv-SFC_delta_total**2)
SFC_spike_delta = np.zeros(num_spike_times_used)
for j in range(num_spike_times_used):
    for k in range(f_delta_diff_num):
        SFC_spike_delta[j]+=SFC_spikes[k+f_delta_low_num,j]/float(f_delta_diff_num)
## theta_band
label_theta='SFC_theta.dat'
file_theta=open(label_theta,'w+')
f_theta_low=4.0
f_theta_high=8.0
f_theta_low_num=int(f_theta_low/df)
f_theta_high_num=int(f_theta_high/df)
f_theta_diff_num = f_theta_high_num-f_theta_low_num
SFC_theta = np.zeros(max_spiketime+1)
SFC_theta_total = 0.0
SFC_theta_stdv = 0.0
for i in range(max_spiketime+1):
    for k in range(f_theta_diff_num):
        SFC_theta[i]+=SFC0[k+f_theta_low_num,i]/float(f_theta_diff_num)
    SFC_theta_total+=SFC_theta[i]/float(max_spiketime+1)
    SFC_theta_stdv+=SFC_theta[i]**2/float(max_spiketime+1)
    str = '%f \n'%SFC_theta[i]
    file_theta.write(str)
file_theta.close()
print("written ",label_theta)
SFC_theta_stdv = np.sqrt(SFC_theta_stdv-SFC_theta_total**2)
SFC_spike_theta = np.zeros(num_spike_times_used)
for j in range(num_spike_times_used):
    for k in range(f_theta_diff_num):
        SFC_spike_theta[j]+=SFC_spikes[k+f_theta_low_num,j]/float(f_theta_diff_num)
## alpha_band
label_='SFC_alpha.dat'
file_=open(label_,'w+')
f_alpha_low=8.0
f_alpha_high=12.0
f_alpha_low_num=int(f_alpha_low/df)
f_alpha_high_num=int(f_alpha_high/df)
f_alpha_diff_num = f_alpha_high_num-f_alpha_low_num
SFC_alpha = np.zeros(max_spiketime+1)
SFC_alpha_total = 0.0
SFC_alpha_stdv = 0.0
for i in range(max_spiketime+1):
    for k in range(f_alpha_diff_num):
        SFC_alpha[i]+=SFC0[k+f_alpha_low_num,i]/float(f_alpha_diff_num)
    SFC_alpha_total+=SFC_alpha[i]/float(max_spiketime+1)
    SFC_alpha_stdv+=SFC_alpha[i]**2/float(max_spiketime+1)
    str = '%f \n'%SFC_alpha[i]
    file_.write(str)
file_.close()
print("written ",label_)
SFC_alpha_stdv = np.sqrt(SFC_alpha_stdv-SFC_alpha_total**2)
SFC_spike_alpha = np.zeros(num_spike_times_used)
for j in range(num_spike_times_used):
    for k in range(f_alpha_diff_num):
        SFC_spike_alpha[j]+=SFC_spikes[k+f_alpha_low_num,j]/float(f_alpha_diff_num)
## beta_band
label_='SFC_beta.dat'
file_=open(label_,'w+')
f_beta_low=12.0
f_beta_high=25.0
f_beta_low_num=int(f_beta_low/df)
f_beta_high_num=int(f_beta_high/df)
f_beta_diff_num = f_beta_high_num-f_beta_low_num
SFC_beta = np.zeros(max_spiketime+1)
SFC_beta_total = 0.0
SFC_beta_stdv = 0.0
for i in range(max_spiketime+1):
    for k in range(f_beta_diff_num):
        SFC_beta[i]+=SFC0[k+f_beta_low_num,i]/float(f_beta_diff_num)
    SFC_beta_total+=SFC_beta[i]/float(max_spiketime+1)
    SFC_beta_stdv+=SFC_beta[i]**2/float(max_spiketime+1)
    str = '%f \n'%SFC_beta[i]
    file_.write(str)
file_.close()
print("written ",label_)
SFC_beta_stdv = np.sqrt(SFC_beta_stdv-SFC_beta_total**2)
SFC_spike_beta = np.zeros(num_spike_times_used)
for j in range(num_spike_times_used):
    for k in range(f_beta_diff_num):
        SFC_spike_beta[j]+=SFC_spikes[k+f_beta_low_num,j]/float(f_beta_diff_num)
## gamma_band
label_='SFC_gamma.dat'
file_=open(label_,'w+')
f_gamma_low=30.0
f_gamma_high=50.0
f_gamma_low_num=int(f_gamma_low/df)
f_gamma_high_num=int(f_gamma_high/df)
f_gamma_diff_num = f_gamma_high_num-f_gamma_low_num
SFC_gamma = np.zeros(max_spiketime+1)
SFC_spike_gamma = np.zeros(num_spike_times_used)
SFC_gamma_total = 0.0
SFC_gamma_stdv = 0.0
for i in range(max_spiketime+1):
    for k in range(f_gamma_diff_num):
        freq = float(k+f_gamma_low_num)*df
        help = SFC0[k+f_gamma_low_num,i]
        SFC_gamma[i]+=help/float(f_gamma_diff_num)
        #print("i=%d k=%f  SFC=%f SFC_sum=%f"%(i,freq,help,SFC_gamma[i]))        
    SFC_gamma_total+=SFC_gamma[i]/float(max_spiketime+1)
    SFC_gamma_stdv+=SFC_gamma[i]**2/float(max_spiketime+1)
    str = '%f \n'%SFC_gamma[i]
    file_.write(str)
file_.close()
print("written ",label_)
SFC_gamma_stdv = np.sqrt(SFC_gamma_stdv-SFC_gamma_total**2)
SFC_spike_gamma = np.zeros(num_spike_times_used)
for j in range(num_spike_times_used):
    for k in range(f_gamma_diff_num):
        SFC_spike_gamma[j]+=SFC_spikes[k+f_gamma_low_num,j]/float(f_gamma_diff_num)
        
endtime = max_spiketime*dt   

print("SFC: delta=%f (%f)   theta=%f (%f)    alpha=%f (%f)       beta=%f (%f)       gamma=%f  (%f)"\
    %(  SFC_delta_total,SFC_delta_stdv,\
        SFC_theta_total,SFC_theta_stdv,\
        SFC_alpha_total,SFC_alpha_stdv,\
        SFC_beta_total,SFC_beta_stdv,\
        SFC_gamma_total,SFC_gamma_stdv))

#quit()
#print(SFC)
fign=plt.figure()
ax = fign.add_subplot(111)
#plt.imshow(SFC, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime, fmin_, fmax_])
SFC_out = SFC0[fmin_num:fmax_num,:]
plt.imshow(SFC_out, interpolation='nearest', cmap=plt.cm.jet,origin='lower',vmin=0.0,vmax=1.0,extent=[0,max_spiketime, fmin, fmax])
forceAspect(ax,aspect=1)
plt.colorbar()
#ax = fign.add_subplot(212)
#plt.plot(range(max_spiketime+1),SFC_gamma,'.')
#plt.plot(spike_times_used[:],SFC_spike_gamma[:],'.')
#ax = fign.add_subplot(613)
##plt.plot(range(max_spiketime+1),SFC_beta,'.')
#plt.plot(spike_times_used[:],SFC_spike_beta[:],'.')
#ax = fign.add_subplot(614)
##plt.plot(range(max_spiketime+1),SFC_alpha,'.')
#plt.plot(spike_times_used[:],SFC_spike_alpha[:],'.')
#ax = fign.add_subplot(615)
##plt.plot(range(max_spiketime+1),SFC_theta,'.')
#plt.plot(spike_times_used[:],SFC_spike_theta[:],'.')
#ax = fign.add_subplot(616)
##plt.plot(range(max_spiketime+1),SFC_delta,'.')
#plt.plot(spike_times_used[:],SFC_spike_delta[:],'.')
#fign.subplots_adjust(hspace=0.8)
savefig('result.png')
plt.show()            
            
        
        