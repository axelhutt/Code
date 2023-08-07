import matplotlib
#matplotlib.use("TkAgg")
#import matplotlib.backends.backend_gtkagg
import matplotlib.pyplot as plt
import numpy as np

import scipy.special as special
import params_class as own_classes
import functions as func
import compute_tf as tf



def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)    

# define basic class of parameters
P = own_classes.params()



# initialize signal object to treat signals of dynamics
number_of_time_points = 1000
sample_and_integration_time_step = 0.01
signal = own_classes.signal(  
    number_of_time_points, 
    sample_and_integration_time_step
    )
    
# initialize global dynamics object to investigate spatial mean dynamics
#gd = own_classes.global_dynamics(signal)
#gd.simulation()
#quit()

#signal.compute_timefrequencyplots()
#signal.plot_tfresults()

# initialize field object to simulate network dynamics
field = own_classes.field(signal)
field.initialize_kernels()
matrix = field.compute_connectivity_spectrum()
evalue,evector = np.linalg.eig(matrix)
plt.plot(np.real(evalue),np.imag(evalue),'k.')
f = open("evals.dat","w+")
for k in range(field.num):
    str = '%f %f\n'%(np.real(evalue[k]),np.imag(evalue[k]))
    f.write(str)
f.close()
#plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.jet,origin='lower')#,vmin=0.0,vmax=1.0,extent=[0,max_spiketime, fmin, fmax])
#plt.colorbar()
plt.show()
quit()




num_matrix = np.shape(matrix)[0]
mean = 0.0
var = 0.0
for k in range(num_matrix):
    for l in range(num_matrix):
        mean += matrix[k,l]/float(num_matrix**2)
        var += matrix[k,l]**2/float(num_matrix**2)
var = var - mean**2
print("matrix: mean=%f var=%f"%(mean,var))


sigma = 0.5
num = 2000
p=0.95
pn = float(num)*p
print("pn=%f"%pn)
K = np.zeros((num,num))
diag = np.zeros(num)
for k in range(num):
    for ll in range(num-k-1):
        l = k+ll+1
        rand = np.random.uniform(0.0,1.0)
        if rand<=p:
            K[k,l]=1.0
            #K[k,l]=np.random.normal(0,sigma)#/np.sqrt(num)
for k in range(num):
    for l in range(k):
        ##symmetric
        K[k,l]=K[l,k]
        ## nonsymmetric
        #rand = np.random.uniform(0.0,1.0)
        #if rand<=p:
        #    K[k,l]=1.0
            
    K[k,k]=0.0
    diag[k]=K[k,k]
#print("diagonal: max=%f min=%f"%(np.max(diag),np.min(diag)))
evalue,evector = np.linalg.eig(K)

sumev = np.sum(evalue)
maxevalue = np.max(evalue)
print("max eigenvalue=%f sum_ev=%f"%(np.max(evalue),sumev))
max = -100000.0
for k in range(np.size(evalue)):
    if np.real(evalue[k])>max and np.real(evalue[k])<maxevalue:
        max = np.real(evalue[k])
        kmax2evalue = k
        vecmax2=evector[:,k]

        
mean = 0.0
var = 0.0
for k in range(num):
    for l in range(num):
        mean += K[k,l]/float(num**2)
        var += K[k,l]**2/float(num**2)
var = var - mean**2
stddev = np.sqrt(var)
print("matrix: mean=%f std=%f    %f"%(mean,stddev,np.sqrt(p*(1-p))))        
sigma = stddev

### from Vjačeslav L Girko. Theory of Random Determinants. Mathematics and Its Applications. Kluwer Academic Publishers, Dordrecht, The Netherlands, 1990. :: lambda_max=2*sigma*sqrt(num)
evalue2ndtheory = 2*sigma*np.sqrt(num)#*np.sqrt(p)


print("kmax2evalue=%d  evalue2nd:"%kmax2evalue,evalue[kmax2evalue]," theory:",evalue2ndtheory,"  difference=",np.abs(evalue[kmax2evalue]-evalue2ndtheory)/evalue2ndtheory)  



degree = np.zeros(num,dtype=int)
for k in range(num):
    counter = 0
    for l in range(num):
        if K[k,l]>0.0:
        #if matrix[k,l]>0.0:
                degree[k]+=1
kmax = np.max(degree)
kmin = np.min(degree)
kmean = np.mean(degree)
kkvar = np.var(degree)+kmean**2
hetero = kkvar/kmean**2
print("kmax=%d kmin=%d  mean=%f hetero=%f "%(kmax,kmin,kmean,hetero))
num_bin = 40
dk = float(kmax-kmin)/float(num_bin)
khist = np.zeros(num_bin+1)
bins = np.linspace(kmin,kmax,num_bin+1) 
for k in range(num):
    #print("k=%d"%k)
    bin = int((degree[k]-kmin)/dk)
    khist[bin] += 1
#print("shapes: bins=",np.shape(bins)," hist=",np.shape(khist))
#plt.bar(bins,khist)
plt.plot(np.real(evalue),np.imag(evalue),'k.')
#plt.imshow(K, interpolation='nearest', cmap=plt.cm.jet,origin='lower')#,vmin=0.0,vmax=1.0,extent=[0,max_spiketime, fmin, fmax])
#plt.colorbar()
plt.show()
quit()



fig = plt.figure(1)
plt.plot(np.real(evalue),np.imag(evalue),'k.')

fig = plt.figure(2)
ax = fig.add_subplot(211)
plt.imshow(np.real(evector), interpolation='nearest', cmap=plt.cm.jet,origin='lower')#,vmin=0.0,vmax=1.0,extent=[0,max_spiketime, fmin, fmax])
forceAspect(ax,aspect=2)
plt.colorbar()
ax = fig.add_subplot(212)
plt.imshow(np.imag(evector), interpolation='nearest', cmap=plt.cm.jet,origin='lower')#,vmin=0.0,vmax=1.0,extent=[0,max_spiketime, fmin, fmax])
forceAspect(ax,aspect=2)
plt.colorbar()

for k in range(np.size(evalue)):
    if evalue[k]>0.5:
        vecmax=evector[:,k]
        kmaxevalue=k
        #print("k=%d  evalue:"%k,evalue[k])
        #print("evector_max:",vecmax)
        break
maxevalue = evalue[kmaxevalue]
max = -100000.0
for k in range(np.size(evalue)):
    if np.real(evalue[k])>max and np.real(evalue[k])<maxevalue:
        max = np.real(evalue[k])
        kmax2evalue = k
        vecmax2=evector[:,k]
print("kmax2evalue=%d  evalue2nd:"%kmax2evalue,evalue[kmax2evalue])    
fig = plt.figure(3)
print(vecmax2)
plt.plot(np.real(vecmax),np.imag(vecmax),'r.',np.real(vecmax2),np.imag(vecmax2),'k.')    
#plt.plot(np.real(vecmax2),np.imag(vecmax2),'k.')    

plt.show()

#field.simulation()
#signal.plot_spacetimeresults()
#signal.plot_globalresults()
