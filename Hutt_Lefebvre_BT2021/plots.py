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

class plots:
		
	def plot_ts_power_bandpass(self):
		time   = np.linspace(0,self.num_time*self.dt,self.num_time) # time in hours
		time_s = np.linspace(0,self.num_time_s*self.dt_s,self.num_time_s) # time in hours
		
		fig=plt.figure(1)
		rows=3
		cols=1
		plt.subplot(rows,cols,1)
		plt.plot(time_s,self.expcond,'k')    
		plt.subplot(rows,cols,2)
		plt.plot(time,self.u,'b')    
		plt.subplot(rows,cols,3)
		plt.plot(time,self.Ue[self.initialtime:self.initialtime+self.num_time],'r')    
		#plt.xticks(np.arange(1,P.Nnat_ass,step),np.arange(0,end_time,tstep))
		
		#fig=plt.figure(2)
		#rows=3
		#cols=1
		#plt.subplot(rows,cols,1)
		#plt.plot(self.freqs_p_n,self.p_n,'b')
		#plt.subplot(rows,cols,2)
		#plt.plot(self.freqs_p_n,self.p_CT_n,'r')
		#plt.subplot(rows,cols,3)
		#plt.plot(self.freqs_p_n,self.p_EEG_n,'k')
		
		fig=plt.figure(3)
		rows=4
		cols=1
		min_ = min(np.min(self.EEG_gammaband),np.min(self.EEG_alphaband))
		max_ = max(np.max(self.EEG_gammaband),np.max(self.EEG_alphaband))
		yticks_pos = [min_,0,max_]
		min_label = '%.3f'%min_
		max_label = '%.3f'%max_
		yticks_label = [min_label,'0',max_label]
		
		plt.subplot(rows,cols,1)
		plt.plot(time_s,self.effective_function,'k')
		
		plt.subplot(rows,cols,2)
		plt.plot(time_s,self.EEG_gammaband,'b')
		plt.yticks(yticks_pos,yticks_label)
		
		plt.subplot(rows,cols,3)
		plt.plot(time_s,self.EEG_alphaband,'r')
		plt.yticks(yticks_pos,yticks_label)
		
		plt.subplot(rows,cols,4)
		plt.plot(time_s,self.EEG_fullband,'k')
		yticks_pos = [np.min(self.EEG_fullband),0,np.max(self.EEG_fullband)]
		yticks_full_label = ['%.3f'%np.min(self.EEG_fullband),'0','%.3f'%np.max(self.EEG_fullband)]
		plt.yticks(yticks_pos,yticks_full_label)
		
		label = 'result_timeseries.dat'
		f = open(label,"w")
		for k in range(self.num_time_s):
			str = '%f    %f %f %f %f\n'%(time_s[k],self.effective_function[k],self.EEG_gammaband[k],self.EEG_alphaband[k],self.EEG_fullband[k])
			f.write(str)
		f.close()
		
		
		fig=plt.figure(4)
		rows=2
		cols=1
		plt.subplot(rows,cols,1)
		plt.plot(self.freqs_p_gammaband,self.p_gammaband,'b')
		plt.subplot(rows,cols,2)
		plt.plot(self.freqs_p_alphaband,self.p_alphaband,'r')

		fig=plt.figure(5)
		rows=2
		cols=1
		plt.subplot(rows,cols,1)
		plt.plot(time_s[:self.num_time_s-self.instpower_window],self.expcond[:self.num_time_s-self.instpower_window],'k')
		plt.subplot(rows,cols,2)
		plt.plot( time_s[:self.num_time_s-self.instpower_window],self.instpower_gamma,'b',\
				  time_s[:self.num_time_s-self.instpower_window],self.instpower_alpha,'r')
		#plt.subplot(rows,cols,2)
		#plt.plot(self.freqs_p_alphaband,self.p_alphaband,'r')
		label = 'result_instpower.dat'
		f = open(label,"w")
		for k in range(self.num_time_s-self.instpower_window):
			str = '%f    %f %f %f \n'%(time_s[k],self.effective_function[k],self.instpower_gamma[k],self.instpower_alpha[k])
			f.write(str)
		f.close()

		plt.show() 

	def plot_ts_power(self):
		time = np.linspace(0,self.num_time*self.dt,self.num_time) # time in hours
		time_s = np.linspace(0,self.num_time_s*self.dt_s,self.num_time_s) # time in hours
		
		fig=plt.figure(1)
		rows=2+1
		cols=1
		plt.subplot(rows,cols,1)
		plt.plot(time_s,self.expcond,'k')    
		plt.subplot(rows,cols,2)
		plt.plot(time,self.u,'b')    
		plt.subplot(rows,cols,3)
		plt.plot(time,self.Ue[self.initialtime:self.initialtime+self.num_time],'r')    
		#plt.xticks(np.arange(1,P.Nnat_ass,step),np.arange(0,end_time,tstep))
		
		label = 'plot_timeseries.png'
		plt.savefig(label)
		
#		fig=plt.figure(2)
#		rows=3
#		cols=1
#		plt.subplot(rows,cols,1)
#		plt.plot(self.freqs_p_n,self.p_n,'k')
#		plt.subplot(rows,cols,2)
#		plt.plot(self.freqs_p_n,self.p_CT_n,'b')
#		plt.subplot(rows,cols,3)
#		plt.plot(self.freqs_p_n,self.p_EEG_n,'r')
		
		plt.show() 

	def plot_power_segments(self):
		fig=plt.figure(11)
		cols = 2
		rows = max(self.list_length_down,self.list_length_up)
		sum_power_up = self.power_up[0]*0.0
		sum_power_down = self.power_down[0]*0.0
		counter_up = 0
		for k in range(self.list_length_up):
			plotnum = 2*k+1
			plt.subplot(rows,cols,plotnum)	
			plt.plot(self.freqs_up[k],self.power_up[k],'b')    
			argmax_gamma = np.argmax(self.power_up[k])
			if len(self.power_up[k])==len(self.power_up[0]):
				sum_power_up = sum_power_up + self.power_up[k]
				counter_up += 1
		sum_power_up /= float(counter_up)
		print("number of valid UP states for power averaging: %d"%counter_up)
		
		
		counter_down = 0
		print("self.list_length_down=%d number of power spectra:%d"%(self.list_length_down,len(self.power_down)))
		for k in range(len(self.power_down)):
			plotnum = 2*k+2
			plt.subplot(rows,cols,plotnum)	
			#plt.plot(self.freqs_down[k],np.log(self.power_down[k]),'k')  
			plt.plot(self.freqs_down[k],self.power_down[k],'r')  
			argmax_alpha = np.argmax(self.power_down[k])
			print("alpha-peak at %f"%(self.freqs_down[k])[argmax_alpha])  
			if len(self.power_down[k])==len(self.power_down[0]):
				sum_power_down = sum_power_down + self.power_down[k]
				counter_down += 1
		sum_power_down /= float(counter_down)
		print("number of valid DOWN states for power averaging: %d"%counter_down)
		#plt.xticks(np.arange(1,P.Nnat_ass,step),np.arange(0,end_time,tstep))
		label = 'plot_power_in_segments.png'
		plt.savefig(label)
		
		fig=plt.figure(12)
		cols = 2
		rows = 1
		plt.subplot(rows,cols,1)	
		#plt.plot(self.freqs_up[0],np.log(sum_power_up),'k') 
		plt.plot(self.freqs_up[0],sum_power_up,'b') 
		plt.subplot(rows,cols,2)	
		#plt.plot(self.freqs_down[0],np.log(sum_power_down),'k') 
		plt.plot(self.freqs_down[0],sum_power_down,'r') 
		label = 'plot_avgpower_segments.png'
		plt.savefig(label)
		
		fig=plt.figure(13)
		cols = 1
		rows = 1
		plt.subplot(rows,cols,1)	
		#plt.plot(self.freqs_up[0],np.log(sum_power_up),'k',self.freqs_down[0],np.log(sum_power_down),'b') 
		plt.plot(self.freqs_up[0],sum_power_up,'b',self.freqs_down[0],sum_power_down,'r') 
		label = 'plot_avgpower_UP+DOWN.png'
		plt.savefig(label)
		
		plt.show()	