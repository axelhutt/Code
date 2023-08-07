import matplotlib
#matplotlib.use("TkAgg")
#import matplotlib.backends.backend_gtkagg
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
#import params_class as params
import functions as func

class params:
    def __init__(self):
        self.c = 10.0
        self.c2 = 10.0
        self.theta = 0.0
        self.alpha = 200.0 # excitatory synaptic rate
        self.beta  = 50.0  # inhibitory synaptic rate
        self.noisevariance_v  = 0.2/self.beta#1.0#0.5
        self.noisevariance_u_1 = 0.4/self.alpha#0.4#0.32
        self.noisevariance_u_2 = 0.00001/self.alpha#0.4#0.32
        self.p1                = 0.6 ## weight for noise class 1
        self.p2                = 1.0-self.p1 ## weight for noise class 2
        self.S10    = 1.7
        self.Fstat  = 2.176321 # but * S10= 3.689
        self.Mstat  = 3.876005 
        self.F     = 0.4
        self.M     = 0.4
        self.I10   = 1.1
        self.I20   = 0.4
        self.p     = 1.0
        self.R1 = 1.0#0.4
        self.R2 = 1.0#0.6
        self.V0 = 0.0
        self.W0 = 0.0
        self.valthr = -100.0
        self.deltamu = 0.0 # distance between Gaussian noise means
        self.mu1     = self.deltamu
        self.mu2     = -self.mu1
        self.inputrate    = 10.0 # firing rate of synaptic input
        self.inputstrength = 2.1 # synaptic coupling of input
        self.inputtau      = 1.0/self.alpha # time constant of input synapse
        self.input_mean    = self.inputstrength*self.inputrate*self.inputtau
        self.input_variance = self.inputstrength*self.input_mean/2.0
        self.z_e = 0
        #self.dt = 0.005


P = params()    

#### simulation of activity with one parameter set
num_time = 10000
dt = 0.01


P.F=P.Fstat#2.176321 
P.M=P.Mstat#3.876005 
    
### lower state
#Vmin=-0.6
#Vmax=0.5
#Wmin=-0.3
#Wmax=4.3
### center state
#Vmin=-0.2
#Vmax=0.7
#Wmin=0.0
#Wmax=5.0
### upper state
#Vmin=0.6
#Vmax=1.4
#Wmin=4.2
#Wmax=5.0
### center + upper state
Vmin=0.5
Vmax=1.2
Wmin=2.0
Wmax=4.9
### lower+upper state
    #Vmin=-0.7
    #Vmax=1.5
    #Wmin=-0.0
    #Wmax=5.0

V=Vmin 
W=Wmin

inputrate_min = 0.0
inputrate_max = 0.0

Imin=1.10#1.3#P.I10#1.25
Imax=1.10#1.3#P.I10#2.25

I2min = 0.4#0.3
I2max = 0.4#0.3

p0 = 1.0
p1 = 1.0

R20 = 1.0#0.3
R21 = 1.0#0.3*p1
deltamu0 = 0.0#0.05
deltamu1 = 0.0#0.05
noisevar0=0.2/P.alpha
noisevar1=0.4/P.alpha
num=10#100
V,W=func.find_multiple_stationary_solutions(V,W,Vmin,Vmax,Wmin,Wmax,Imin,Imax,I2min,I2max,deltamu0,deltamu1,noisevar0,noisevar1,p0,p1,R20,R21,inputrate_min,inputrate_max,num,P)

