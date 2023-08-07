import numpy as np
import compute_tf 
import compute_gd
import compute_field
import scipy.special as special

class params:
    def __init__(self):
        self.c = 10.0
        self.c2 = 10.0
        self.theta = 0.0
        self.theta_micro = 0.0
        self.noisevariance_v  = 0.2#1.0#0.5
        self.noisevariance_u_1 =  0.0 # Poisson-noise     #0.436
        self.noisevariance_u_2 =  0.0001#0.32
        self.S10    = 1.7
        self.Fstat  = 2.176321 # but * S10= 3.689
        self.Mstat  = 3.876005 
        self.p1    = 1.0
        self.p2    = 1.0-self.p1
        self.F     = 0.4
        self.M     = 0.4
        
        ### topo = 0: non-symmetric ER with non-zero diagonal 
        ### topo = 1: non-symmetric ER with zero diagonal 
        ### topo = 2: symmetric ER with zero diagonal 
        self.topo  = 0  
        
        ## synaptic time scales
        self.tau_exc = 0.005
        self.tau_inh = 0.020
        
        ## Poisson noise parameters
        self.weight_in = 2.1
        
        self.noisevariance_u_1_scale = self.noisevariance_u_1*self.tau_exc
        self.noisevariance_u_2_scale = self.noisevariance_u_2*self.tau_exc
        self.noisevariance_v_scale = self.noisevariance_v*self.tau_inh
        
        self.I10   = 1.1
        self.I20   = 0.4
        self.V0 = 0.0
        self.w0 = 0.0
        self.valthr = -100.0
        self.deltamu = 0.0 # distance between Gaussian noise means
        self.mu1     = self.deltamu
        self.mu2     = -self.mu1
        #self.dt = 0.005
#  def initialize(self):

    ## firing rate of a single neuron
    def f(self,V):
        return 1.0/((1.0+np.exp(-200.0*(V-self.theta_micro)))*self.dt)

    ## spike generating function
    def H1(self,V):
        help = 0.0
        ## Heaviside function
        # help = self.S10/(1.0+np.exp(-500.0*(V)))
        ## Poisson-distributed spike
        ran = np.random.uniform(0,1)
        rate = self.f(V)
        #print("ran=%f rate*dt=%f"%(ran,rate*self.dt))
        if ran <= rate*self.dt:
            help = self.S10
        return help 
    
    ## spike generating function
    def H2(self,V):
        help = 0.0
        ## Heaviside function
        # help = self.S10/(1.0+np.exp(-500.0*(V)))
        ## Poisson-distributed spike
        ran = np.random.uniform(0,1)
        rate = self.f(V)
        if ran <= rate*self.dt:
            help = 1.0
        return help 


    def S1(self,V):
        ## self.noisevariance_u is correct when it takes the unscaled values
        return ( \
            self.S10*0.5*(1.0-special.erf((self.theta-(V-self.mu1))/(np.sqrt(2*self.noisevariance_u_1/self.tau_exc)))) \
            +self.S10*0.5*(1.0-special.erf((self.theta-(V-self.mu2))/(np.sqrt(2*self.noisevariance_u_1/self.tau_exc)))) \
            )/2.0
    
    def S1_deriv(self,V):
        return ( \
            self.S10*np.exp(-np.power(V-self.mu1-self.theta,2)/(2.0*self.noisevariance_u_1))/(np.sqrt(2*np.pi*self.noisevariance_u_1)) \
            + self.S10*np.exp(-np.power(V-self.mu2-self.theta,2)/(2.0*self.noisevariance_u_1))/(np.sqrt(2*np.pi*self.noisevariance_u_1)) \
            )/2.0

    def S2(self,V):
        return 0.5*(1.0-special.erf((self.theta-V)/(np.sqrt(2*self.noisevariance_v))))

    def S2_deriv(self,V):
        return np.exp(-np.power(V-self.theta,2)/(2.0*self.noisevariance_v))/(np.sqrt(2*np.pi*noisevariance_v))




class signal(params,compute_tf.tf_methods):
    def __init__(self,num_time,dt):
        params.__init__(self)
        self.num_time = num_time
        self.dt       = dt
        self.df       = 1.0/(self.num_time*self.dt)
        self.fmin     = self.df
        self.fmax     = 1.0

    def compute_timefrequencyplots(self):
        print ("compute_timefrequencyplots")
        #data = self.VV
        #self.data = data
        self.power=self.cwt_tf()

    def plot_tfresults(self):
        self.time_max=self.num_time*self.dt
        self.time=np.arange(0.0,self.time_max,self.dt)
        self.power=self.cwt_tf()
        #self.power=self.BartlettWelch()
        #print np.shape(self.time)
        #print np.shape(self.power)
        self.plot_results()

    def plot_spacetimeresults(self):
        self.time_max=self.num_time*self.dt
        self.time=np.arange(0.0,self.time_max,self.dt)
        print (np.shape(self.time))
        print (np.shape(self.fielddynamics))
        self.plot_field_results()


class global_dynamics(params,compute_gd.global_dynamics):

    def __init__(self,signal):
        params.__init__(self) ## inherit class and include member variables in own class
        self.signal   = signal
        self.dt       = self.signal.dt
        self.num_time = self.signal.num_time

    def plot_globalresults(self):
        self.time_max=self.num_time*self.dt
        #del self.time
        self.time=np.arange(0.0,self.time_max,self.dt)
        #del self.data
        self.data=self.signal.globaldynamics
        
        fc=5 # central frequency for Morlet wavelets
        df = 1
        
        ## for this wavelet function, it is f=0.8125*fs/scale
        fc = 0.8125
        fmin=2.0
        fmax=80.0#fs/2.0 # Nyquist
        num_freqs=int((fmax-fmin)/df)
        frequencies = np.arange(fmin,fmax,df)
        scales = fc/(frequencies*dt)
        #print(scales)
        
        wavelet_type = 'morl'
        f = pywt.scale2frequency(wavelet_type, scales)/dt # 
        print("f:",f)
        fmin_ = np.min(f)
        fmax_ = np.max(f)
        coefs, freqs = pywt.cwt(data=self.data, scales=scales, wavelet=wavelet_type, sampling_period=self.dt,method='fft')
        #plt.matshow(coefs) 
        print(freqs)
    
        self.power = np.abs(coefs)
        self.powerlog=np.log(abs(coefs))
        endtime = duration*dt
        
        print("data-shape:",np.shape(power),"  num_segment:%d"%num_segment)
        
        fign=plt.figure()
        ax = fign.add_subplot(111)
        #plt.imshow(power_mean, interpolation='nearest', cmap=plt.cm.jet,origin='lower', vmin=-2, vmax=-0.6, extent=[0,endtime, fmin_, fmax_])
        plt.imshow(self.powerlog, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime, fmin_, fmax_])
        forceAspect(ax,aspect=1)
        plt.colorbar()        
        plt.show()
        
        
        

class field (params,compute_field.field):

    def __init__(self,signal):
        params.__init__(self) ## inherit class and include member variables in own class
        self.num = 200#1000#600#200
        self.dt  = signal.dt
        self.num_time = signal.num_time
        
        ## connectivity parameters
        self.connect_thresh = 0.95#0.95
        self.alpha1  = 1.0*self.Fstat/self.connect_thresh# valid if radii are num/2
        self.beta1   = -0.0
        self.alpha2  = -1.0*self.Mstat/self.connect_thresh
        self.beta2   = 0.0
        self.alpha3  = -1.0*self.Fstat/self.connect_thresh
        self.beta3   = 0.0
        self.alpha4  = 1.0*self.Mstat/self.connect_thresh
        self.beta4   = 0.0
                
        self.signal  = signal
        

            
            