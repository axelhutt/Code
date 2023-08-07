import matplotlib
matplotlib.use("TkAgg")
#import matplotlib.backends.backend_gtkagg
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special


def H(V):
    if V>=0.0:
        return 1.0
    else:
        return 0.0

def dge(u):
    df = -u/tau
    return df



class params:
    def __init__(self,num,num_K):
        self.num = num # dimension of OU
        self.num_K = num_K # dimension of network
        self.mu = np.zeros(num)
        self.D = 1.0
        self.tau = 0.0
        self.minu = 0.0
        self.maxu = 0.0
        self.bin_width = 1.0
        self.num_bins = 1
    
    def generate_ER_network(self): ## symmetric ER with zero diagonal
        sigma = 0.5
        num = self.num_K
        p=0.4
        pn = float(num)*p
        print("pn=%f"%pn)
        self.K = np.zeros((num,num))
        #diag = np.zeros(num)
        for k in range(num):
            for ll in range(num-k-1):
                l = k+ll+1
                rand = np.random.uniform(0.0,1.0)
                if rand<=p:
                    self.K[k,l]=1.0/pn
                    #K[k,l]=np.random.normal(0,sigma)#/np.sqrt(num)
        for k in range(num):
            for l in range(k):
                ##symmetric
                self.K[k,l]=self.K[l,k]
                ## nonsymmetric
                #rand = np.random.uniform(0.0,1.0)
                #if rand<=p:
                #    K[k,l]=1.0
            self.K[k,k]=0.0
            #diag[k]=K[k,k]
        #print("diagonal: max=%f min=%f"%(np.max(diag),np.min(diag)))
        #evalue,evector = np.linalg.eig(K)
    
    def compute_OUmean(self):
        self.generate_ER_network()
        evalue,evector0 = np.linalg.eig(self.K)
        evector = evector0*np.sqrt(float(self.num_K))#normalize to norm=1
        E2=np.var(evector,axis=0)+np.mean(evector,axis=0)**2
        #print(np.mean(evector,axis=0))
        #print(np.var(evector,axis=0))
        #print(E2)
        self.ER_lambda = evalue[1:]
        self.ER_meanevec = np.zeros(self.num)
        for k in range(self.num):
            self.ER_meanevec[k]=np.mean(evector[:,k+1])
        print("mean_vec=%g std_vec=%g"%(np.mean(self.ER_meanevec),np.sqrt(np.var(self.ER_meanevec))))
        
        self.product_mean = np.zeros(self.num)
        for k in range(self.num):
            self.product_mean[k] = self.ER_lambda[k]*self.ER_meanevec[k]
        print("mean_product=%f std_product=%f"%(np.mean(self.product_mean),np.sqrt(np.var(self.product_mean))))

        fig = plt.figure(20)
        ax = fig.add_subplot(311)
        plt.plot(range(self.num),self.ER_lambda,'k.')
        ax = fig.add_subplot(312)
        plt.plot(range(self.num),self.ER_meanevec,'k')
        ax = fig.add_subplot(313)
        plt.plot(range(self.num),self.product_mean,'k')
        plt.show()
    
    def distrib_lambda(self,lmin,lmax):
        self.p_lambda = np.zeros(self.num_bins_red)
        z = np.zeros(self.num_bins_red)
        lambdamax2 = lmax**2
        sum = 0.0
        k_lmin=int((lmin-self.minu_red)/self.bin_width_red)
        k_lmax=int((lmax-self.minu_red)/self.bin_width_red)
        num = k_lmax-k_lmin
        print("lmax-lmin=%f num=%d  k_lmin=%d    num_bins=%d"%(lmax-lmin,num,k_lmin,self.num_bins_red))
        if num<=1:
            print("!!!!!!!!!!! lmax-lmin=%f  <= bin_width_red=%f"%(lmax-lmin,self.bin_width_red))
            exit()
        for kk in range(num):
            lambda_=lmin+(lmax-lmin)*float(kk)/float(num)
            k = k_lmin+kk
            self.p_lambda[k] = np.sqrt(lambdamax2-lambda_**2)
            #print("kk=%d lambda_=%f sqrt=%f"%(kk,lambda_,np.sqrt(lambdamax2-lambda_**2)))
            sum += self.p_lambda[k]
        self.p_lambda /= float(sum)
        sum = 0.0
        for k in range(self.num_bins_red):
            z[k]=self.minu_red+float(k)*self.bin_width_red
            #print("z=%f p_lambda=%f"%(z[k],self.p_lambda[k]))
            sum+= self.p_lambda[k]
        print("sum of p_lambda:",sum)
        #plt.figure(10)
        #plt.plot(z,self.p_lambda,'b')
        #plt.show()
    
    def distrib_evectors(self,xmin,xmax):
        self.p_evectors = np.zeros(self.num_bins_red)
        z = np.zeros(self.num_bins_red)
        var = 1.0/float(self.num)
        sum = 0.0
        art_mean = 0.00
        k_xmin=int((xmin-self.minu_red)/self.bin_width_red)
        k_xmax=int((xmax-self.minu_red)/self.bin_width_red)
        num = k_xmax-k_xmin
        print("reconstructed: xmin=%f xmax=%f"%(float(k_xmin)*self.bin_width_red+self.minu_red,float(k_xmax)*self.bin_width_red+self.minu_red))
        for kk in range(num):
            x = xmin + (xmax-xmin)*float(kk)/float(num)
            k = k_xmin + kk
            self.p_evectors[k] = np.exp(-((x-art_mean)**2)/(2.0*var))/(np.sqrt(np.pi*var))
            sum += self.p_evectors[k]
        self.p_evectors /= sum
        sum = 0.0
        for k in range(self.num_bins_red):
            z[k]=self.minu_red+float(k)*self.bin_width_red
            sum+= self.p_evectors[k]
        print("sum of p_evectors:",sum)
        #plt.figure(11)
        #plt.plot(z,self.p_evectors,'r')
        #plt.show()
        
        
    def distrib_mean(self):
        lmin = -0.5
        lmax = -lmin
        xmin = -1.0
        xmax = -xmin
        
        self.num_bins_red = 1000
        self.minu_red = -1
        self.maxu_red = -self.minu_red
        self.bin_width_red = (self.maxu_red-self.minu_red)/float(self.num_bins_red)
        
        self.distrib_lambda(lmin,lmax)
        self.distrib_evectors(xmin,xmax)
        self.mu_analytic = np.zeros(self.num_bins_red)
        self.z = np.zeros(self.num_bins_red)
        sum = 0.0
        for k in range(self.num_bins_red):
            self.mu_analytic[k]=0.0
            self.z[k] = self.minu_red+float(k)*self.bin_width_red
            for kk in range(self.num_bins_red):
                x = self.minu_red+float(kk)*self.bin_width_red
                q=self.z[k]/x
                #print("k=%d kk=%d x=%f q=%f"%(k,kk,x,q))
                absx = np.abs(x)
                if absx<1.e-3: 
                    kkk = self.num_bins_red
                    absx=1.e-3
                else:
                    kkk = int((q-self.minu_red)/self.bin_width_red)
                if kkk<0: kkk=0
                if kkk>=self.num_bins_red: kkk=self.num_bins_red-1
                help = self.p_lambda[kkk]*self.p_evectors[kk]/absx
                self.mu_analytic[k] += help 
                #print("kk=%d kkk=%d help=%f      %f %f    x=%f  q=%f"%(kk,kkk,help,self.p_lambda[kkk],self.p_evectors[kk],x,q))
                #print("k=%d mu_analytic=%f"%(k,self.mu_analytic[k]))
            sum += self.mu_analytic[k]
        sum_new = 0.0
        for k in range(self.num_bins_red):
            self.mu_analytic[k]/=sum
            sum_new += self.mu_analytic[k]
        print("mu analytic:   sum=%f sum_new=%f"%(sum,sum_new))
        #plt.figure(0)
        #plt.plot(self.z,self.p_lambda,'b',self.z,self.p_evectors,'r',self.z,self.mu_analytic,'k-')
        #plt.show()
            
            
    def pdf_analytic(self,v):
        sum = 0.0
        for k in range(self.num_bins_red):
            help = np.exp(-((v-self.mu_analytic[k])**2)/(2.0*P.D/P.tau))/np.sqrt(6.28*P.D/P.tau)
            sum += help/float(self.num_bins_red)
        #print("v=%f sum=%f"%(v,sum))
        return sum
    
    def sample_mu(self,num):
        from numpy.random import choice
        num_samples = num
        samples = choice(self.z, num_samples, p=self.mu_analytic)
        self.mu = samples
        max_samples=np.max(samples)
        min_samples=np.min(samples)
        num_bin = 50
        bwidth=(max_samples-min_samples)/float(num_bin)
        hist = np.zeros(num_bin+1)
        for k in range(num_samples):
            bin = int((samples[k]-min_samples)/bwidth)
            hist[bin]+= 1.0/float(num_samples)
        #plt.figure(13)
        #plt.plot(np.linspace(min_samples,max_samples,num_bin+1),hist,'k')
        #plt.show()

num_K = 200
num = num_K-1
num_time=20000
P = params(num,num_K)

dt = 0.0005
D_2 = 0.0001
tau = 0.005
P.tau = tau
mu1 = 0.0
mu2 = -mu1

D_1 = 0.50
P.D = D_1
percentage = 1.0

mean_p_mean = 5.0
mean_p_std  = 10.5


P.compute_OUmean()


P.distrib_mean()
#P.sample_mu(num)
P.mu = P.product_mean
#quit()
## choose random locations of two different noise processses
u = np.zeros(num)
mu_loc = np.zeros(num,dtype=int)
finished = 0
counter = 0
while finished == 0:
    val = np.random.randint(0,num,1)
    #print(counter,val)
    if mu_loc[val]==0:
        mu_loc[val]=1
        counter +=1
        if counter == int(percentage*num):
            finished = 1

#for k in range(num):
#P.mu[k]=np.random.normal(mean_p_mean,mean_p_std)

## initialize noise process u(t)
for k in range(num-1):
    #
    #if mu_loc[k]==1:
    #    u[k] = mu1
    #else:
    #    u[k] = mu2
    u[k] = P.mu[k]

## compute noise process u(t)
ut = np.zeros((num_time,num))
for i in range(num_time):
    for k in range(num):
        
        du  = dge(u[k])
        muu = P.mu[k]
        D   = D_1
        #if mu_loc[k]==1:
        #    muu = mu1
        #    D   = D_1
        #else:
        #    muu = mu2
        #    D   = D_2
            
        u[k] = u[k] + dt*(du+muu/tau) + np.sqrt(dt*2.0*D)*np.random.normal(0.0,1.0)/tau
        ut[i,k]=u[k]
    #print "%d u=%f" %(k+1,u[k+1])


## compute probability density function of u
min_u=np.min(np.min(ut))
max_u=np.max(np.max(ut))
P.minu = min_u*1.0
P.maxu = max_u*1.0
min_ut=np.min(ut,axis=1)
max_ut=np.max(ut,axis=1)
num_bin = 100#50
P.num_bins = num_bin
bin_width = (max_u-min_u)/float(num_bin)
P.bin_width = bin_width*1.0
bin_widtht = (max_ut-min_ut)/float(num_bin)
print ("min=%f max=%f bin_width=%f" %(min_u,max_u,bin_width))

hist = np.zeros(num_bin+1)
pdf  = np.zeros(num_bin+1)
histt=np.zeros((num_time,num_bin+1))
for i in range(num_time):
    for k in range(num):  
        bin = int((ut[i,k]-min_ut[i])/bin_widtht[i])
        histt[i,bin]+=1.0/float(num)

hist_sum = 0.0
mean_raw=0.0
stddev_raw=0.0
for i in range(num_time):
    for k in range(num):  
        mean_raw+=ut[i,k]/float(num*num_time)
        stddev_raw+=ut[i,k]*ut[i,k]/float(num*num_time)
        bin = int((ut[i,k]-min_u)/bin_width)
        hist[bin]+=1.0/float(num*num_time)
        #print "%d u[k]=%f  bin=%d hist=%f" %(k,u[k],bin,hist[bin])
variance=stddev_raw-mean_raw*mean_raw
for bin in range(num_bin):
    hist_sum+=hist[bin]
print ("raw: mean=%f variance=%f     hist_sum=%f" %(mean_raw,variance,hist_sum))



vv = np.zeros(num_bin+1)
label='out_pdfs_q%.2f.dat'%percentage
fp=open(label,'w')
pdf_sum=0.0
for bin in range(num_bin+1):
    v = float(bin)*bin_width+min_u
    vv[bin]=v
    #pdf[bin]= np.exp(-((v-mu1)**2)/(2.0*variance))*bin_width/np.sqrt(6.28*variance)
    pdf[bin]= P.pdf_analytic(v)*P.bin_width
    pdf_sum = pdf_sum + pdf[bin]
    str = '%f %f %f\n' %(v,hist[bin],pdf[bin])
    fp.write(str)
fp.close()
print ("written %s"%label)
print ('pdf_sum = %f' %pdf_sum)
plt.figure(1)
plt.plot(vv,hist,'b',vv,pdf,'r')
plt.show()

quit()

## compute theoretical distribution and sample distribution of sigmoid 
num1 = 1000
Vmin=-20.0#min_u
Vmax=20.0#max_u
dV=(Vmax-Vmin)/float(num1)
V = np.arange(Vmin,Vmax,dV)
sigmoid = np.zeros(num1)
theory = np.zeros(num1)
for l in range(num1):
    theory[l]=0.5*(percentage*(1.0-special.erf(-(V[l]-mu1)/(np.sqrt(2*D_1/tau))))+ \
                    (1.0-percentage)*(1.0-special.erf(-(V[l]-mu2)/(np.sqrt(2*D_2/tau)))))

sigmoid_mean     = np.zeros(num1)
mean_deviation   = np.zeros(num1)
stddev_deviation = np.zeros(num1)
for i in range(num_time):
    if i%100 == 0:
        print ("i(%d)=%d" %(num_time,i))
    sigmoid[:]=0.0
    for l in range(num1):
        for bin in range(num_bin):
            v = float(bin)*bin_widtht[i]+min_ut[i]
            sigmoid[l]+=H(V[l]+v)*histt[i,bin]
        sigmoid_mean[l]     += sigmoid[l]/float(num_time)
        mean_deviation[l]   += (theory[l]-sigmoid[l])/float(num_time)
        stddev_deviation[l] += pow(theory[l]-sigmoid[l],2)/float(num_time)

label = 'out_sigmoid_q%.2f.dat'%percentage
fp=open(label,'w')
for l in range(num1):
    if stddev_deviation[l]-mean_deviation[l]**2>0:
        hhelp=np.sqrt(stddev_deviation[l]-mean_deviation[l]**2)
    else:
        print ("sqrt < 0: %f" %(stddev_deviation[l]-mean_deviation[l]**2))
        hhelp=0.0
    str = '%f %f  %f     %f %f\n' %(V[l],theory[l],sigmoid_mean[l],mean_deviation[l],hhelp)
    fp.write(str)
fp.close()  
print ("written %s"%label)

## plot results
plt.figure(2)
plt.plot(V,sigmoid,'b',V,theory,'r')
plt.show()



# mean_x = 0
# stddev_x=0.0
# for bin in range(num_bin):
#     x = float(bin)*bin_width+min_u
#     mean_x+=x*hist[bin]#*bin_width
#     stddev_x=stddev_x+x*x*hist[bin]#*bin_width
#     #print "x=%f mean=%f stddev_x=%f   hist=%f %f" %(x,mean_x,stddev_x,hist[bin],hist_sum)
#     hist_sum+=hist[bin]
# stddev_x = np.sqrt(stddev_x-mean_x*mean_x)
# print hist_sum,mean_x,stddev_x

# stddevv = 0.0
# for bin in range(num_bin):
#     x = float(bin)*bin_width+min_u
#     stddevv+=np.power(x-mean_x,2)*hist[bin]*bin_width
#     #print "x=%f stddevv=%f  " %(x,stddevv)
# print "stddev=%f" %np.sqrt(stddevv)


#plt.plot(hist)
#plt.plot(np.arange(min_u,max_u,bin_width),hist[0:num_bin])
#plt.show()