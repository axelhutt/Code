import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
from scipy import signal
#import functions as func
import matplotlib.mlab as mlab
from scipy.signal import butter, lfilter, lfilter_zi, filtfilt
from math import atan

import plots
import routines

class params_CT:
	def __init__(self):
		#self.c = 10.0
		#self.c2 = 10.0
		self.theta_CT = 0.1
		
		self.Dglobal_CT          = 1.e-20#1.e-7

		self.Wee    = 2.21
		self.Wie    = -3.46 
		self.Wei    = 4.58
		self.Wii    = -1.69
		self.Ie     = 0.10#0.1
		self.Ii     = 0.00#0.0
		self.b		= 0.5#0.1
		self.delay  = 0.014#0.010#0.005
		
		self.Ue0 = -0.55
		self.Ui0 = 0.45
		self.Ve0 = 0.66
		self.Vi0 = -0.58
		
		self.alpha_e  = 50.0
		self.alpha_i  = 100.0
		self.a_e      = 5.0#5.0
		self.a_i      = self.a_e
		
		self.mu_CT    = 0.0

		self.DD0 = 1.2#0.1
		self.DD1 = 1.4#1.2
		self.compute_noisevariance_CT_uv(self.DD1)
		self.noisevariance_CT_u1 = self.noisevariance_CT_u
		self.compute_noisevariance_CT_uv(self.DD0)
		self.noisevariance_CT_u0 = self.noisevariance_CT_u
		print("D=%f sigma_e=%f   sigma_i=%f"%(self.DD0,self.noisevariance_CT_u0,self.noisevariance_CT_v))
		


	def compute_noisevariance_CT_uv(self,DD):
		## assumes a_e=a_i
		lambda_plus  = 0.5*(-self.alpha_e-self.a_e+np.sqrt((self.alpha_e-self.a_e)**2+4.0*self.b*self.alpha_e*self.a_e))
		lambda_minus = 0.5*(-self.alpha_e-self.a_e-np.sqrt((self.alpha_e-self.a_e)**2+4.0*self.b*self.alpha_e*self.a_e))
		A 			 =  np.sqrt(2.0*np.pi)*(self.a_e+lambda_plus)/(lambda_plus-lambda_minus)
		B 			 = -np.sqrt(2.0*np.pi)*(self.a_e+lambda_minus)/(lambda_plus-lambda_minus)
		help = -1.4142*DD*(A**2/(2.0*lambda_plus)+2.0*A*B/(lambda_plus+lambda_minus)+B**2/(2.0*lambda_minus))
		self.noisevariance_CT_u = help

		lambda_plus  = 0.5*(-self.alpha_i-self.a_i+np.sqrt((self.alpha_i-self.a_i)**2+4.0*self.b*self.alpha_i*self.a_i))
		lambda_minus = 0.5*(-self.alpha_i-self.a_i-np.sqrt((self.alpha_i-self.a_i)**2+4.0*self.b*self.alpha_i*self.a_i))
		A 			 =  np.sqrt(2.0*np.pi)*(self.a_i+lambda_plus)/(lambda_plus-lambda_minus)
		B 			 = -np.sqrt(2.0*np.pi)*(self.a_i+lambda_minus)/(lambda_plus-lambda_minus)
		help = -1.4142*DD*(A**2/(2.0*lambda_plus)+2.0*A*B/(lambda_plus+lambda_minus)+B**2/(2.0*lambda_minus))
		self.noisevariance_CT_v = help
		
	def Fe(self,V):
		return 0.5*(1.0-special.erf((self.theta_CT-(V-self.mu_CT))/(np.sqrt(2*self.noisevariance_CT_u))))
		
	def Fe_deriv(self,V):
		return np.exp(-np.power(V-self.mu_CT-self.theta_CT,2)/(2.0*self.noisevariance_CT_u))/(np.sqrt(2*np.pi*self.noisevariance_CT_u))
		
	def Fi(self,V):
		return 0.5*(1.0-special.erf((self.theta_CT-V)/(np.sqrt(2*self.noisevariance_CT_v))))

	def Fi_deriv(self,V):
		return np.exp(-np.power(V-self.theta_CT,2)/(2.0*self.noisevariance_CT_v))/(np.sqrt(2*np.pi*noisevariance_CT_v))


class params:
	def __init__(self):
		#self.c = 10.0
		#self.c2 = 10.0
		self.theta = 0.0
		self.noisevariance_v  = 0.5#0.5
		self.noisevariance_u0 = 0.1#0.32
		self.noisevariance_u1 = 0.8#0.2#0.32
		self.noisevariance_u  = self.noisevariance_u0#0.32
		self.Dglobal          = 1.e-4#1.e-5
		self.S10    = 1.7

		self.Fstat  = 2.176321 # but * S10= 3.689
		self.Mstat  = 3.876005 
		self.F     = 0.4
		self.M     = 0.4
		self.I10   = 1.10#1.15#1.45
		self.I20   = 0.4#0.4
		
		## upper state for variance_u=0.4
		#self.V0 = 1.3
		#self.W0 = 4.57
		## lower state for variance_u=0.8
		self.V0 = -0.5
		self.W0 = 0.5
		
		self.alpha  = 200.0
		self.beta 	= 50.0
		
		self.valthr = -100.0
		self.deltamu = 0.0 # distance between Gaussian noise means
		self.mu1     = self.deltamu
		self.mu2     = -self.mu1
		self.mu 	 = 0.0
		#self.dt = 0.005

	def S1(self,V):
		return ( \
			self.S10*0.5*(1.0-special.erf((self.theta-(V-self.mu1))/(np.sqrt(2*self.noisevariance_u)))) \
			+self.S10*0.5*(1.0-special.erf((self.theta-(V-self.mu2))/(np.sqrt(2*self.noisevariance_u)))) \
			)/2.0
	
	def S1_deriv(self,V):
		return ( \
			self.S10*np.exp(-np.power(V-self.mu1-self.theta,2)/(2.0*self.noisevariance_u))/(np.sqrt(2*np.pi*self.noisevariance_u)) \
			+ self.S10*np.exp(-np.power(V-self.mu2-self.theta,2)/(2.0*self.noisevariance_u))/(np.sqrt(2*np.pi*self.noisevariance_u)) \
			)/2.0

	def S2(self,V):
		return 0.5*(1.0-special.erf((self.theta-V)/(np.sqrt(2*self.noisevariance_v))))

	def S2_deriv(self,V):
		return np.exp(-np.power(V-self.theta,2)/(2.0*self.noisevariance_v))/(np.sqrt(2*np.pi*noisevariance_v))


class externalfield (params,params_CT,plots.plots,routines.routines):

	def __init__(self,t):
		params.__init__(self) ## inherit class and include member variables in own class
		params_CT.__init__(self) ## inherit class and include member variables in own class

		self.dt                 = 0.0005
		self.fs  			    = 1.0/self.dt
		self.weight			    = 0.5#0.333 # 1.0: CT only ; 0.0: intracortical only
		self.stimulation_period = 50
		self.sample_factor	    = 4
		self.dt_s               = self.dt*float(self.sample_factor)
		self.fs_s				= 1.0/self.dt_s
		
		self.shortest_interval  = int(2.0*10.0*self.fs) # 2 segments, df=0.1
		self.num_time           = int(self.stimulation_period*1.2*self.shortest_interval)
		self.num_time_s 		= int(self.num_time/self.sample_factor)
		self.shortest_interval_s = int(2.0*10.0*self.fs_s) # 2 segments, df=0.1

		self.DT0			 = 0.2*np.pi# need to be be < pi/2
		self.DT1			 = (np.pi/2-self.DT0)*2.0
		print("DT0=%f DT1=%f"%(self.DT0,self.DT1))
		
		self.initialtime	 = int(self.delay/self.dt)
		self.t               = t
		#print("finished in __init__")

		self.expcond		 = np.zeros(self.num_time_s,dtype=np.int8)
		self.EEG 			 = np.zeros(self.num_time_s)
		self.effective_function = np.zeros(self.num_time_s)
		
		self.Vexcite		 = 0.0

		self.C_Uu			 = 0.01
		
		self.signal_low 	 = []
		self.counter_low = 0

	def initialize(self):
						
		## ---------------  for cortical model 
		D1 = self.noisevariance_u0#1.000
		D2 = self.noisevariance_v#1.000
			
		## u is stored as object variable to be able to map it to signal object and visualize it
		self.u = np.zeros(self.num_time)  
		self.u[0] = self.V0+0.0*np.random.normal(0,1.0)*np.sqrt(D1)     
		self.du = 0.0
		self.v = self.W0+0.0*np.sqrt(D2)*np.random.normal(0.0,1.0)
		self.dv = 0.0

		string_label='./out_global.dat' 
		string_field       = './out_field.dat'
		string_noise = './out_noise.dat'
		self.f = open(string_label,'w')
		self.f_field = open(string_field,'w')
		self.f_noise = open(string_noise,'w')
	
	
		## ---------------  for cortico-thalamic model 
		sigma_e = self.noisevariance_CT_u0#1.000
		sigma_i = self.noisevariance_CT_v#1.000

		## u is stored as object variable to be able to map it to signal object and visualize it
		self.Ue = np.zeros(self.num_time+self.initialtime)  
		self.Ue[0:self.initialtime+1] = self.Ue0+0.0*np.random.normal(0,1.0,self.initialtime+1)     
		self.dUe = 0.0
		self.Ui = np.zeros(self.num_time+self.initialtime)
		self.Ui[0:self.initialtime+1] = self.Ui0+0.0*np.random.normal(0,1.0,self.initialtime+1)     
		self.dUi = 0.0
		self.Ve = self.Ve0+0.0*np.random.normal(0,1.0)     
		self.dVe = 0.0
		self.Vi = self.Vi0+0.0*np.random.normal(0,1.0)     
		self.dVi = 0.0

		self.icounter = 0
		

if  __name__ == '__main__':
	# initialize field object to simulate network dynamics

	
	field = externalfield(0.0)
	field.initialize()
	for i in range(field.num_time-1):
	    field.run(field.Vexcite,i)
	field.close()
	field.writeout()
	field.powerspectrum()
	field.plot_ts_power()
	field.plot_power_segments()
	
	field.bandpassfilter()
	field.plot_ts_power_bandpass()
	