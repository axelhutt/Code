import matplotlib
#matplotlib.use("TkAgg")
#import matplotlib.backends.backend_gtkagg
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
import params_class as own_classes


class field():

    def initialize_kernels(self):
        
        def distance(x):
            L = self.num
            return int(L/2.0-np.abs(L/2.0-np.abs(x)))
        
        
        ## radii=num/2 mandatory for choice of alpha_n,beta_n
        self.radius1 = self.num / 1#4 
        self.radius2 = self.num / 1 
        self.radius3 = self.num / 1 
        self.radius4 = self.num / 1

        self.K1 = np.zeros((self.num,self.num))
        self.K2 = np.zeros((self.num,self.num))
        self.M1 = np.zeros((self.num,self.num))
        self.M2 = np.zeros((self.num,self.num))
        d1 = np.zeros(self.num)
        d2 = np.zeros(self.num)
        d3 = np.zeros(self.num)
        d4 = np.zeros(self.num)
    
        
        #maxK=-1.0 
        #var_sumi  = 0.0
        self.mean_sumK1i = 0.0
        self.mean_sumM1i = 0.0
        self.mean_sumK2i = 0.0
        self.mean_sumM2i = 0.0
        var_sumK1i = 0.0
        var_sumM1i = 0.0
        var_sumK2i = 0.0
        var_sumM2i = 0.0
        
        if self.topo == 0: ## non-symmetrix matrix with random diagonal
            for j in range(self.num):
                if np.abs(j-self.num/2)<=self.radius1:
                    d1[j]=self.alpha1/float(self.num)
                else:
                    d1[j]=self.beta1/float(self.num)                    
                if np.abs(j-self.num/2)<=self.radius2:
                    d2[j]=self.alpha2/float(self.num)
                else:
                    d2[j]=self.beta2/float(self.num)
                if np.abs(j-self.num/2)<=self.radius3:
                    d3[j]=self.alpha3/float(self.num)
                else:
                    d3[j]=self.beta3/float(self.num)
                if np.abs(j-self.num/2)<=self.radius4:
                    d4[j]=self.alpha4/float(self.num)
                else:
                    d4[j]=self.beta4/float(self.num)
            for i in range(self.num):
                for jj in range(self.num-i):
                    j = i+jj
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K1[i,j]=d1[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M1[i,j]=d2[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K2[i,j]=d3[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M2[i,j]=d4[distance(i-j)]
            for i in range(self.num):
                for j in range(i):
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K1[i,j]=d1[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M1[i,j]=d2[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K2[i,j]=d3[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M2[i,j]=d4[distance(i-j)]
                                
        if self.topo == 1: ## non-symmetrix matrix with zero diagonal
            for j in range(self.num):
                if np.abs(j-self.num/2)<=self.radius1:
                    d1[j]=self.alpha1/float(self.num-1)
                else:
                    d1[j]=self.beta1/float(self.num)                    
                if np.abs(j-self.num/2)<=self.radius2:
                    d2[j]=self.alpha2/float(self.num-1)
                else:
                    d2[j]=self.beta2/float(self.num)
                if np.abs(j-self.num/2)<=self.radius3:
                    d3[j]=self.alpha3/float(self.num-1)
                else:
                    d3[j]=self.beta3/float(self.num)
                if np.abs(j-self.num/2)<=self.radius4:
                    d4[j]=self.alpha4/float(self.num-1)
                else:
                    d4[j]=self.beta4/float(self.num-1)
            for i in range(self.num):
                for jj in range(self.num-i-1):
                    j = i+jj+1
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K1[i,j]=d1[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M1[i,j]=d2[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K2[i,j]=d3[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M2[i,j]=d4[distance(i-j)]
            for i in range(self.num):
                for j in range(i):
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K1[i,j]=d1[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M1[i,j]=d2[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K2[i,j]=d3[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M2[i,j]=d4[distance(i-j)]

        if self.topo == 2: ## ER with symmetrix matrix and zero diagonal
            for j in range(self.num):
                if np.abs(j-self.num/2)<=self.radius1:
                    d1[j]=self.alpha1/float(self.num-1)
                else:
                    d1[j]=self.beta1/float(self.num)                    
                if np.abs(j-self.num/2)<=self.radius2:
                    d2[j]=self.alpha2/float(self.num-1)
                else:
                    d2[j]=self.beta2/float(self.num)
                if np.abs(j-self.num/2)<=self.radius3:
                    d3[j]=self.alpha3/float(self.num-1)
                else:
                    d3[j]=self.beta3/float(self.num)
                if np.abs(j-self.num/2)<=self.radius4:
                    d4[j]=self.alpha4/float(self.num-1)
                else:
                    d4[j]=self.beta4/float(self.num-1)
            for i in range(self.num):
                for jj in range(self.num-i-1):
                    j = i+jj+1
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K1[i,j]=d1[distance(i-j)]                        
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M1[i,j]=d2[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.K2[i,j]=d3[distance(i-j)]
                    rand=np.random.uniform(0.0,1.0)
                    if rand<=self.connect_thresh:
                        self.M2[i,j]=d4[distance(i-j)]
            for i in range(self.num):
                for j in range(i):
                    self.K1[i,j]=self.K1[j,i]
                    self.M1[i,j]=self.M1[j,i]
                    self.K2[i,j]=self.K2[j,i]
                    self.M2[i,j]=self.M2[j,i]



        for i in range(self.num):
            sumK1i=0.0
            sumM1i=0.0
            sumK2i=0.0
            sumM2i=0.0
            for j in range(self.num):
                sumK1i=sumK1i+self.K1[i,j]
                sumM1i=sumM1i+self.M1[i,j]
                sumK2i=sumK2i+self.K2[i,j]
                sumM2i=sumM2i+self.M2[i,j]
            self.mean_sumK1i=self.mean_sumK1i+sumK1i/float(self.num)    
            self.mean_sumM1i=self.mean_sumM1i+sumM1i/float(self.num)
            self.mean_sumK2i=self.mean_sumK2i+sumK2i/float(self.num)
            self.mean_sumM2i=self.mean_sumM2i+sumM2i/float(self.num)
            var_sumK1i += sumK1i**2/float(self.num)
            var_sumM1i += sumM1i**2/float(self.num)
            var_sumK2i += sumK2i**2/float(self.num)
            var_sumM2i += sumM2i**2/float(self.num)

            #print "mean_sumK1i=%f mean_sumK2i=%f mean_sumM1i=%f mean_sumM2i=%f" %(mean_sumK1i,mean_sumM1i,mean_sumK2i,mean_sumM2i)
        #stddev=np.sqrt(var_sumi-mean_sumi*mean_sumi)
        var_sumK1i = var_sumK1i - self.mean_sumK1i**2
        var_sumM1i = var_sumM1i - self.mean_sumM1i**2
        var_sumK2i = var_sumK2i - self.mean_sumK2i**2
        var_sumM2i = var_sumM2i - self.mean_sumM2i**2
        
        print ("mean of K1 = %f   theory = %f      %f" %(self.mean_sumK1i,self.Fstat,var_sumK1i))
        print ("mean of M1 = %f   theory = %f      %f" %(self.mean_sumM1i,-self.Mstat,var_sumM1i))
        print ("mean of K2 = %f   theory = %f      %f" %(self.mean_sumK2i,-self.Fstat,var_sumK2i))
        print ("mean of M2 = %f   theory = %f      %f" %(self.mean_sumM2i,self.Mstat,var_sumM2i))
        
        if np.abs(self.mean_sumK1i-self.Fstat)>0.05:
            print ("|mean_sumK1i-self.Fstat|>0.05 !!! %f>0.05 . Quit." %np.abs(self.mean_sumK1i-self.Fstat))
        if np.abs(self.mean_sumM1i+self.Mstat)>0.05:
            print ("|mean_sumM1i+self.Mstat|>0.05 !!! %f>0.05 . Quit." %np.abs(self.mean_sumM1i+self.Fstat))
        if np.abs(self.mean_sumK2i+self.Fstat)>0.05:
            print ("|mean_sumK1i+self.Fstat|>0.05 !!! %f>0.05 . Quit." %np.abs(self.mean_sumK2i+self.Fstat))
        if np.abs(self.mean_sumM2i-self.Mstat)>0.05:
            print ("|mean_sumM2i-self.Mstat|>0.05 !!! %f>0.05 . Quit." %np.abs(self.mean_sumM2i-self.Mstat) )
        

    
    def compute_connectivity_spectrum(self):        
        print ("compute spectrum of connectivity matrices:")
        matrix = -self.M1/self.Mstat
        #evalue,evectors = np.linalg.eig(-self.M1/self.Mstat)
        #print("shape of ev_F",np.shape(ev_F))
        #print ("evalue=",evalue)
        return matrix
     
     
        
    def simulation(self):
        
        def meanfield_integration(i):
            a = self.Vbar[i]
            b = self.Wbar[i]
            self.Vbar[i+1] = a+self.dt*(-a+self.Fstat*self.S1(a-self.meannoise_u)-self.Mstat*self.S2(b)+self.I10)/self.tau_exc\
                                +np.sqrt(self.dt*2.0)*np.random.normal(0.0,np.sqrt(self.noisevariance_u_1/(self.tau_exc*float(self.num))))
            self.Wbar[i+1] = b+self.dt*(-b+self.Mstat*self.S1(a-self.meannoise_u)-self.Fstat*self.S2(b)+self.I20)/self.tau_inh\
                                +np.sqrt(self.dt*2.0)*np.random.normal(0.0,np.sqrt(self.noisevariance_v/(self.tau_inh*float(self.num))))
            
        def dglsyst(uu,vv,duu,dvv):
            for k in range(self.num):
                sumu=0.0
                sumv=0.0
                mean_u=0.0
                stddev_u=0.0
#                 mean_v=0.0
#                 stddev_v=0.0
                mean_sigu = 0.0
#                mean_sigv = 0.0
#                 mean_M1i=0.0
                for l in range(self.num):
                    spike11 = self.H1(uu[l])
                    spike21 = self.H2(vv[l])
                    spike22 = self.H2(vv[l])
                    spike12 = self.H1(uu[l])
                    #print("k=%d l=%d uu=%f spike11=%f"%(k,l,uu[l],spike11))
                    sumu=sumu + self.K1[k,l]*spike11+self.M1[k,l]*spike21
                    sumv=sumv + self.K2[k,l]*spike22+self.M2[k,l]*spike12
                    mean_u+=uu[l]/float(self.num)
                    stddev_u+=uu[l]*uu[l]/float(self.num)
                    mean_sigu+=spike11/float(self.num)
#                     mean_v+=vv[l]/float(self.num)
#                     stddev_v+=vv[l]*vv[l]/float(self.num)
#                     mean_sigv+=self.H2(vv[l])/float(self.num)
                    #mean_K1i+=self.K1[k,l]
                stddev_u=np.sqrt(stddev_u-mean_u*mean_u)
#                 stddev_v=np.sqrt(stddev_v-mean_v*mean_v)
                duu[k] = (-uu[k] + sumu + self.I10)
                dvv[k] = (-vv[k] + sumv + self.I20)
                #if k==0:
                    #print ("sumu=%f -uu=%f duu=%f mean_u=%f mean_sigu=%g theory=%f  var_u=%f     K1=%f" %(sumu,-uu[k],duu[k],mean_u,mean_sigu,self.S1(mean_u),stddev_u**2,self.K1[k,0]))
#                     print "sumv=%f I20=%f -vv=%f dvv=%f mean_v=%f mean_sigv=%f theory=%f  stddev_v=%f" %(sumv,self.I20,-vv[k],dvv[k],mean_v,mean_sigv,self.S2(mean_v),stddev_v)
            return duu,dvv
        
        print("routine simulation.....")
        
        #self.I10=1.45#1.9#1.6
        #self.noisevariance_u=0.4
        
        if self.p1==0.5:
            ## upper branch
            V0=0.69
            W0=4.6
            varu_min = 0.35
            varu_max = 0.35
            ## lower branch
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.55
            #varu_max = 0.55
        
        if self.p1==0.6:
            ## upper branch
            V0=0.69
            W0=4.6
            varu_min = 0.25
            varu_max = 0.35#0.25
            ## lower branch
            #V0=-0.55#
            #W0=0.0#
            #varu_min = 0.33
            #varu_max = 0.33
            
            #V0=0.7
            #W0=4.6
            #varu_min = 0.20
            #varu_max = 0.20
            #V0=0.7
            #W0=4.6
            #varu_min = 0.25
            #varu_max = 0.25
            #V0=0.7
            #W0=4.6
            #varu_min = 0.30
            #varu_max = 0.30
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.35
            #varu_max = 0.35
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.40
            #varu_max = 0.40
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.45
            #varu_max = 0.45
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.50
            #varu_max = 0.50
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.80
            #varu_max = 0.80
            #V0=-0.55
            #W0=0.0
            #varu_min = 1.20
            #varu_max = 1.20
            #V0=-0.55
            #W0=0.0
            #varu_min = 1.80
            #varu_max = 1.80
            #V0=-0.55
            #W0=0.0
            #varu_min = 3.00
            #varu_max = 3.00
            #V0=-0.55
            #W0=0.0
            #varu_min = 5.00
            #varu_max = 5.00
            #V0=-0.55
            #W0=0.0
            #varu_min = 15.00
            #varu_max = 15.00
            
            
        if self.p1==0.8:
            ## upper branch
            V0=1.0
            W0=4.6
            varu_min = 0.20
            varu_max = 0.20
            ## lower branch
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.25
            #varu_max = 0.25
        
        if self.p1==1.0:
            ## upper branch
            V0=1.0
            W0=4.9
            varu_min = 0.10#0.17#0.15
            varu_max = 0.50#0.22#0.15
            ## lower branch
            #V0=-0.55
            #W0=0.0
            #varu_min = 0.20
            #varu_max = 0.20
            
        ## Poisson noise
        #varu_min = 0.0
        #varu_max = 0.0
        lambda_in_min = 0.0#0.08
        lambda_in_max = 0.0#0.11
        
        #V0=-0.55#-0.6#1.3#1.03#-0.18#-0.28
        #W0=0.0#0.32#4.8#0.4#0.34
        #varu_min = 0.33
        #varu_max = 0.33
        

        D1_1 = self.noisevariance_u_1_scale#1.000
        D1_2 = self.noisevariance_u_2_scale#1.000
        D2 = self.noisevariance_v_scale#1.000
        
        ###### mean-field equation variables
        self.Vbar = np.zeros(self.num_time)
        self.Wbar = np.zeros(self.num_time)
        
        ## u is stored as object variable to be able to map it to signal object and visualize it
        self.u = np.zeros((self.num_time,self.num))
        self.signal.fielddynamics = np.transpose(self.u)
        self.signal.dt=self.dt
        self.signal.num_time = self.num_time
        
        du = np.zeros(self.num)

        # V includes the global dynamics and is mapped to the signal object
        self.V = np.zeros(self.num_time)
        self.signal.globaldynamics = self.V
        
        noise = np.zeros(self.num)
        mu_loc = np.zeros(self.num,dtype=int)
        finished = 0
        counter = 0
        label_classes = "./classes.dat"
        f_classes = open(label_classes,"w+")
        while finished == 0:
            val = np.random.randint(0,self.num-1,1)
            if mu_loc[val]==0:
                mu_loc[val]=1 ## self.p1*num elements belong to class 1
                str = "%d \n"%val
                f_classes.write(str)
                counter +=1
                #print ("random location #%d set in location %d" %(counter,val))
                if counter == int((self.num-1)*self.p1):
                    finished = 1
        f_classes.close()
        print("written %s"%label_classes)
        
        sumu_=0.0
        varu_=0.0
        sum_sigu_ = 0.0
        for k in range(self.num):
            if mu_loc[k]==1:
                #self.u[0,k] = V0+np.random.normal(self.mu1,np.sqrt(D1))
                noise[k]=np.random.normal(self.mu1,np.sqrt(D1_1))
            else:
                #self.u[0,k] = V0+np.random.normal(self.mu2,np.sqrt(D1))
                noise[k]=np.random.normal(self.mu2,np.sqrt(D1_2))
            self.u[0,k] = V0+noise[k]
            sum_sigu_ +=self.H1(self.u[0,k])/float(self.num)
            sumu_ = sumu_+noise[k]/float(self.num)
            varu_ = varu_+noise[k]**2/float(self.num)
        varu_=varu_-sumu_**2
        print ("test: mean of u=%f V0=%f diff=%f    var of u=%f mu1=%f mu2=%f   mean sigu=%f theory=%f" %(sumu_+V0,V0,sumu_,varu_,self.mu1,self.mu2,sum_sigu_,self.S1(sumu_+V0)))
        self.Vbar[0]=np.mean(self.u[0,:])
        
        
        v = np.zeros(self.num)
        dv = np.zeros(self.num)
        #W = W0
        #Wana = W0
        for k in range(self.num):
            v[k] = W0+np.sqrt(D2)*np.random.normal(0.0,0.1)
        self.Wbar[0]=np.mean(v)

        string_label       = './out_global_p%.3f_sigma%.2f.dat'%(self.p1,varu_min)
        string_label_spike = './out_spike_p%.3f_sigma%.2f.dat'%(self.p1,varu_min)
        string_label_spike_1 = './out_spike_C1_p%.3f_sigma%.2f.dat' %(self.p1,varu_min)
        string_label_spike_2 = './out_spike_C2_p%.3f_sigma%.2f.dat' %(self.p1,varu_min)
        string_field       = './out_field_p%.3f_sigma%.2f.dat' %(self.p1,varu_min)
        string_field_1     = './out_field_C1_p%.3f_sigma%.2f.dat' %(self.p1,varu_min)
        string_field_2     = './out_field_C2_p%.3f_sigma%.2f.dat' %(self.p1,varu_min)
        f = open(string_label,'w')
        f_spike = open(string_label_spike,'w')
        f_spike_1 = open(string_label_spike_1,'w')
        f_spike_2 = open(string_label_spike_2,'w')
        f_field = open(string_field,'w')
        f_field_1 = open(string_field_1,'w')
        f_field_2 = open(string_field_2,'w')
        
        for i in range(self.num_time-1):
            
            #self.mu1 = 0.4+(0.8-0.4)*float(i)/float(self.num_time)
            #self.mu2 = -self.mu1
            
            self.noisevariance_u_1=(varu_min+(varu_max-varu_min)*float(i)/float(self.num_time))*self.tau_exc
            
            lambda_in        = lambda_in_min+(lambda_in_max-lambda_in_min)*float(i)/float(self.num_time)
            poisson_mean     = lambda_in*self.weight_in*self.tau_exc
            poisson_variance = poisson_mean*self.weight_in/2.0
            self.mu1         = poisson_mean
            self.noisevariance_u_1 += poisson_variance
            
            #self.noisevariance_u_1_scale = self.noisevariance_u_1*self.tau_exc
            D1_1 = self.noisevariance_u_1#1.000
            

            #if i>500:#self.num_time/2:#self.num_time/10:
                #self.mu1 = 0.25
                #self.mu2 = -self.mu1
                #self.sigma1=0.8#+0.4*float(i-self.num_time/2)/float(self.num_time/2)
                #D1 = self.sigma1**2
            #if i%100==0:
            #    print 'i(%d)=%d ' %(self.num_time,i)
            
            du,dv=dglsyst(self.u[i,:],v,du,dv)
            
            sum1 = 0.0
            sum2 = 0.0
            stddev_u=0.0
            stddev_v=0.0
            loc = int(self.num/2)
            #if i%10 == 0:
            #    print ("i=%d(%d)  u=%f du*dt=%f  w=%f dw*dt=%f" %(i,self.num_time-1,self.u[i,loc],du[loc]*self.dt,v[loc],dv[loc]*self.dt))
            
            string = ''
            string_C1 = ''
            string_C2 = ''
            for k in range(self.num):
                if mu_loc[k]==1:#k%2==0: # if k is even
                    #print "k=%d" %k
                    noise[k]=np.random.normal(0,np.sqrt(D1_1))#+np.random.normal(0,np.sqrt(poisson_variance))
                    noise_mean = self.p1*self.mu1
                else:
                    noise[k]=np.random.normal(0,np.sqrt(D1_2))
                    noise_mean = self.mu2
                    #noise_u = np.sqrt(self.dt*2.0*D1)*np.random.normal(0.0,1.0)+self.mu2
                noise_v = np.random.normal(0.0,1.0)
                
                self.u[i+1,k] = self.u[i,k] + (self.dt*du[k] + np.sqrt(self.dt*2.0)*noise[k]+noise_mean*self.dt)/self.tau_exc
                v[k]          = v[k]        + (self.dt*dv[k] + np.sqrt(self.dt*2.0*D2)*noise_v)/self.tau_inh
                
                sum1 = sum1 + self.u[i+1,k]
                sum2 = sum2 + v[k]
                stddev_u += self.u[i+1,k]**2/float(self.num)
                stddev_v += v[k]**2/float(self.num)
                if self.H1(self.u[i+1,k])>0.0 and k%1==0:
                    str = '%d %d\n' %(i+1,k)
                    f_spike.write(str)
                    if mu_loc[k]==1:
                        f_spike_1.write(str)
                    else:
                        f_spike_2.write(str)
                string0 = ' %f ' % self.u[i+1,k]
                string = string + string0
                
                if mu_loc[k]==1:
                    string_C10 = ' %f ' % self.u[i+1,k]
                    string_C1 = string_C1 + string_C10
                else:
                    string_C20 = ' %f ' % self.u[i+1,k]
                    string_C2 = string_C2 + string_C20
                
            self.meannoise_u = np.mean(noise)
            self.meannoise_v = 0.0
            string = string + '\n'
            string_C1 = string_C1 + '\n'
            string_C2 = string_C2 + '\n'
            
            f_field.write(string)
            f_field_1.write(string_C1)
            f_field_2.write(string_C2)
            
            sum1=sum1/float(self.num)
            sum2=sum2/float(self.num)
            stddev_u=np.sqrt(stddev_u-sum1**2)
            stddev_v=np.sqrt(stddev_v-sum2**2)
            mean_u=sum1
            mean_v=sum2
            self.V[i]=mean_u
            
            meanfield_integration(i)
            
#             for k in range(self.num):
#                 noise = np.random.normal(0.0,1.0)
#                 u[k]=u[k]+(-u[k]*self.dt)+np.sqrt(D*self.dt)*noise
#                 mean_u=mean_u+u[k]/float(self.num)
#                 stddev_u=stddev_u+u[k]*u[k]/float(self.num)
#                 noise = np.random.normal(0.0,1.0)
#                 v[k]=v[k]+(-v[k]*self.dt)+np.sqrt(D*self.dt)*noise
#                 mean_v=mean_v+v[k]/float(self.num)
#                 stddev_v=stddev_v+v[k]*v[k]/float(self.num)
#             stddev_u=np.sqrt(stddev_u-mean_u*mean_u)
#             stddev_v=np.sqrt(stddev_v-mean_v*mean_v)
            if i%20==0:
                print ('%d noisevariance=%f  Lin=%f       V=%f (%f)  W=%f  varu=%f varv=%f ' %(i,self.noisevariance_u_1,lambda_in,mean_u,self.Vbar[i],mean_v,stddev_u**2,stddev_v**2))
            #print '%d mean_u=%f  stddev_u=%f ' %(i,mean_u,stddev_u)
            str = '%f %f %f    %f %f     %f %f\n' %(i*self.dt,self.noisevariance_u_1,lambda_in,mean_u,mean_v,self.Vbar[i+1],self.Wbar[i+1])
            #str = '%f  %f  %f\n' %(i*dt,mean_u,stddev_u)
            f.write(str)        
        
        self.V[self.num_time-1]=mean_u
        f.close()
        f_spike.close()    
        f_spike_1.close()    
        f_spike_2.close()    
        f_field.close()
        f_field_1.close()
        f_field_2.close()
        