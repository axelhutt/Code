import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import scipy.fftpack
import matplotlib.mlab as mlab
from cwt_modules_Bergner import *

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)





## ******************* time-frequency plots
## simulated EEG


print ('read file...')

num = -1
time_counter=0
timecounter_min=0000
timecounter_max=32000
f = open('out.dat','r')
for line in f:
    line=line.strip()
    columns = line.split()
    if time_counter>= timecounter_min and time_counter<timecounter_max:
        if time_counter == 0:
            num = np.size(columns)
        time_counter += 1
    else:
        break
f.close()
num_neurons = num
ntimes0 = time_counter
data0=np.zeros((ntimes0,num_neurons))

time_counter = 0
f = open('out.dat','r')
for line in f:
    line=line.strip()
    columns = line.split()
    if time_counter>= timecounter_min and time_counter<timecounter_max:
        for k in range(num_neurons):
            data0[time_counter,k]=float(columns[k])
        time_counter += 1
    else:
        break
f.close()

## define tiem step
dt=0.0005


# now downsampling 
counter = 0
counter_sample = 0
print ('downsample data...')
sample_factor=1
ntimes=int(ntimes0/sample_factor)
data=np.zeros((ntimes,num_neurons))
for i in range(ntimes):
    data[i,:]=data0[i*sample_factor,:]
dt=dt*sample_factor
dataavg = np.mean(data,axis=1)

print ('\n perform time-frequency analysis....')

timemax=ntimes*dt
endtime=timemax
df=1.0/timemax
fs=1.0/dt
fc=8 # central frequency for Morlet wavelets

fmin=5.0
fmax=80.0#fs/2.0 # Nyquist
num_freqs=int((fmax-fmin)/df)+1
fmax=int(num_freqs*df)+fmin
frequencies = np.arange(fmin,fmax,df)
scales = fc/(frequencies*dt)
#print np.shape(frequencies), num_freqs

power=np.zeros((num_freqs,ntimes))

# choose arbitrarily certain number of neurons
randnum_neurons=20
dangle=np.zeros((randnum_neurons,num_freqs,ntimes))
res1=np.zeros((num_freqs,ntimes))
res2=np.zeros((num_freqs,ntimes))

powerthreshold=0.5

for k in range(randnum_neurons):
    neuron_number1 = np.random.random_integers(0,num_neurons-1,1)
    neuron_number2 = np.random.random_integers(0,num_neurons-1,1)
    #print ' neuron no. %d and no.%d' %(neuron_number1,neuron_number2)
    #print np.max(data[:,neuron_number1]), np.min(data[:,neuron_number1])
    #print np.max(data[:,2]), np.min(data[:,2])
    data_tf1=np.squeeze(data[:,neuron_number1])
    data_tf2=np.squeeze(data[:,neuron_number2])

    #print '1'
    #print np.shape(data_tf) , np.shape(scales)

    res1 = cwt( data_tf1 , scales, morlet, fc,dt)[0]
    res2 = cwt( data_tf2 , scales, morlet, fc,dt)[0]
    #maxpower_freq1=np.max(np.abs(res1),axis=0)
    #maxpower_freq2=np.max(np.abs(res2),axis=0)
    #minpower_freq1=np.min(np.abs(res1),axis=0)
    #minpower_freq2=np.min(np.abs(res2),axis=0)
    #print '1.1'
    
    #print np.shape(angle[neuron_number,:,:]), np.shape(np.angle(res))
    dangle[k,:,:] = np.angle(res1)-np.angle(res2)
    power = power + (np.abs(res1)+np.abs(res2))/2
    #print '2'

powersingle=np.abs(cwt( np.squeeze(data[:,99]) , scales, morlet, fc,dt)[0])
poweravg=np.abs(cwt( dataavg , scales, morlet, fc,dt)[0])
power = power / randnum_neurons
maxpower_freq=np.squeeze(np.max(power,axis=0))
minpower_freq=np.squeeze(np.min(power,axis=0))
diffpower_freq=maxpower_freq-minpower_freq

##### plot result
#powerlog=np.log(power) 

print ('compute PLV.....')

PLV_d = 0.0
f_d_min = 1.0
f_d_max = 4.0
PLV_t = 0.0
f_t_min = 4.0
f_t_max = 8.0
PLV_a = 0.0
f_a_min = 8.0
f_a_max = 12.0
PLV_b = 0.0
f_b_min = 12.0
f_b_max = 25.0
PLV_g = 0.0
f_g_min = 25.0
f_g_max = 80.0
dff = 1.0

windowduration=10
#print ntimes-windowduration
PLV = np.zeros((num_freqs,ntimes-windowduration))
helpcos=np.zeros(num_freqs)
helpsin=np.zeros(num_freqs)
for i in range(ntimes-windowduration):
    helpcos[:]=0.0
    helpsin[:]=0.0
    for w in range(windowduration):
        for k in range(randnum_neurons):
            help = np.squeeze(dangle[k,:,i+w])
            helpcos=helpcos+np.cos(help)/(randnum_neurons*windowduration)
            helpsin=helpsin+np.sin(help)/(randnum_neurons*windowduration)
            #if (w==0) & (k==0) & (i<5):
                #print 'i=%d min PLV=%f  max PLV = %f help=%f %f      help1=%f %f' %(i,np.min(PLV[:,i]),np.max(PLV[:,i]),np.max(help), np.min(help),np.max(help1), np.min(help1))
                #print 'i=%d help=%f     helpcos=%f' %(i,help[num_freqs/2],helpcos[num_freqs/2])    #print np.min(PLV[:,i]), np.max(PLV[:,i])
    help=sqrt(helpcos*helpcos+helpsin*helpsin)
    help1=ma.masked_where((np.squeeze(power[:,i])-minpower_freq[i])/diffpower_freq[i]<powerthreshold,help) 
    #if i<3:
    #    print 'i=%d %f %f' %(i,minpower_freq[i],powerthreshold*maxpower_freq[i])
    #    print power[:,i]
    #    print help1
    PLV[:,i]=help1.filled(0)
    
    if i%100 == 0:
        print ('i=%d PLV=%f %f' %(i,np.max(PLV[:,i]),np.min(PLV[:,i])))
    
    # delta
    num_f = int((f_d_max-f_d_min)/dff)
    k0 = int(f_d_min/dff)
    for k in range(num_f):
        PLV_d+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
    # theta
    num_f = int((f_t_max-f_t_min)/dff)
    k0 = int(f_t_min/dff)
    for k in range(num_f):
        PLV_t+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
    # alpha
    num_f = int((f_a_max-f_a_min)/dff)
    k0 = int(f_a_min/dff)
    for k in range(num_f):
        PLV_a+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
    # beta
    num_f = int((f_b_max-f_b_min)/dff)
    k0 = int(f_b_min/dff)
    for k in range(num_f):
        PLV_b+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
    # gamma
    num_f = int((f_g_max-f_g_min)/dff)
    k0 = int(f_g_min/dff)
    for k in range(num_f):
        PLV_g+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
        

print("gPLV: delta=%f theta=%f alpha=%f beta=%f gamma=%f"%(PLV_d,PLV_t,PLV_a,PLV_b,PLV_g))

fign=plt.figure()
ax = fign.add_subplot(211)
plt.imshow(poweravg[:,0:ntimes-windowduration], interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-dt*(ntimes-windowduration), fmin, fmax])
plt.colorbar()
forceAspect(ax,aspect=2)

ax = fign.add_subplot(212)
plt.imshow(PLV, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-windowduration*dt, fmin, fmax])
#plt.imshow(PLV, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-windowduration*dt, fmin, fmax])
forceAspect(ax,aspect=2)
plt.colorbar()

plt.show()


