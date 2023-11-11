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

class routines:
	
	def dglsyst(self,i):

		####### ----------- cortico-thalamic loop -----------------
		#i+self.initialtime : current time moment
		Ue		 = self.Ue[i+self.initialtime]
		Ui		 = self.Ui[i+self.initialtime]
		Ue_delay = self.Ue[i] # Ue(t-tau), so self.Ue[i+self.initialtime-self.initialtime]
		Ui_delay = self.Ui[i] # Ue(t-tau)

		self.dUe = self.alpha_e*(-Ue + self.b*Ue + self.Wee*self.Fe(Ue_delay) + self.Wie*self.Fi(Ui_delay) + self.C_Uu*self.u[i] + self.Ie)
		self.dUi = self.alpha_i*(-Ui + self.b*Ui + self.Wei*self.Fe(Ue_delay) + self.Wii*self.Fi(Ui_delay) + self.Ii)
		self.dVe = self.a_e*(-self.Ve-Ue)
		self.dVi = self.a_i*(-self.Vi-Ui)

		####### ----------- cortical dynamics ------------------
		self.du = self.alpha*(-self.u[i] + self.Fstat*self.S1(self.u[i]) - self.Mstat*self.S2(self.v) + self.I10)
		self.dv = self.beta *(-self.v - self.Fstat*self.S2(self.v) + self.Mstat*self.S1(self.u[i]) + self.I20)
		#print("u:  %f %f %f %f   du=%f   %f"%(-self.u[i],self.Fstat*self.S1(self.u[i]),-self.Mstat*self.S2(self.v),self.I10,self.du/self.alpha))
		#print("v:  %f %f %f %f   dv=%f"%(-self.v,-self.Fstat*self.S2(self.v),self.Mstat*self.S1(self.u[i]),self.I20,self.dv/self.beta))

		#help1 = -V-(P.R1*P.F*S1(V,P)-P.R1*P.M*S2(W,P)+P.R1*P.I10)
		#help2 = -W-(-P.R2*P.F*S2(W,P)+P.R2*P.M*S1(V,P)+P.R2*P.I20)


	def own_switch0(self,phi0):
		phi = phi0%(2.0*np.pi)
		if phi>=0 and phi<self.DT0:
			help = 1.0
		if phi>=self.DT0 and phi<self.DT0+self.DT1:
			help = 1-2.0*(phi-self.DT0)/self.DT1
		if phi>=self.DT0+self.DT1 and phi<3.0*self.DT0+self.DT1:
			#print("phi0=%f self.DT0+self.DT1=%f 3*self.DT0+self.DT1=%f"%(phi0,self.DT0+self.DT1,3.0*self.DT0+self.DT1))
			help = -1.0
		if phi>=3.0*self.DT0+self.DT1 and phi<2.0*np.pi-self.DT0:	
			help = -1.0+2.0*(phi-(2.0*np.pi-self.DT1-self.DT0))/self.DT1
		if phi>=2.0*np.pi-self.DT0 and phi<2.0*np.pi:	
			help = 1.0
		return help
		
		
	def run(self,Vexcite,i):
					
		arg = self.stimulation_period*2*np.pi*float(i)/float(self.num_time)

		## jump between 0 and 1
		#switch0 = np.cos(arg)
		#switch0 = signal.square(arg)
		
		switch0 = (2*1.0/np.pi)*atan(np.sin(arg)/0.01);
				
		#switch0 = self.own_switch0(arg)
		#print("i=%d arg=%f switch0=%f"%(i,arg,switch0))
		switch = (1+switch0)/2.0
		#switch = 0.0
		
		self.noisevariance_u=self.noisevariance_u0+(self.noisevariance_u1-self.noisevariance_u0)*switch
		
		DD = self.DD0+(self.DD1-self.DD0)*switch
		self.compute_noisevariance_CT_uv(DD)
		#print("i=%d  DD=%f sigma_e=%f sigma_i=%f"%(i,DD,self.noisevariance_CT_u,self.noisevariance_CT_v))
		#self.noisevariance_CT_u=self.noisevariance_CT_u0+(self.noisevariance_CT_u1-self.noisevariance_CT_u0)*(1+signal.square(arg))/2.0
		
		self.dglsyst(i) ## computes self.du and self.dv
		if i%10000==0:
			print ("i=%d (%d) u=%f du=%f  w=%f dw=%f Ue=%f dUe=%f Ui=%f dUi=%f   Ve=%f dVe=%f Vi=%f dVi=%f" %(i,\
			self.num_time,\
			self.u[i],self.du,self.v,self.dv,\
			self.Ue[i+self.initialtime],self.dUe,self.Ui[i+self.initialtime],self.dUi,\
			self.Ve,self.dVe,self.Vi,self.dVi))

		### step forward in time
		gn = np.sqrt(2.0*self.Dglobal*self.dt)*np.random.normal(0.0,1.0)
		gn_CT = np.sqrt(2.0*self.Dglobal_CT*self.dt)*np.random.normal(0.0,1.0)
		#print("gn=%f gn_CT=%f"%(gn,gn_CT))
		
		self.u[i+1] 				  = self.u[i] + self.dt*self.du + self.alpha*gn
		self.v      				  = self.v   + self.dt*self.dv
		self.Ue[i+1+self.initialtime] = self.Ue[i+self.initialtime] + self.dt*self.dUe + self.alpha_e*gn_CT
		self.Ui[i+1+self.initialtime] = self.Ui[i+self.initialtime] + self.dt*self.dUi 
		self.Ve						  = self.Ve + self.dt*self.dVe		
		self.Vi						  = self.Vi + self.dt*self.dVi		

		self.Vexcite = self.u[i+1]
		if i%self.sample_factor == 0:
			
			self.EEG[self.icounter]  = (1.0-self.weight)*self.u[i]+self.weight*self.Ue[self.initialtime+i]
			
			## experimental condition
			if switch>=0.5:
				self.expcond[self.icounter]=1
			else:
				self.expcond[self.icounter]=0
				self.signal_low.append(self.Ue[self.initialtime+i])
				#self.signal_low.append(self.EEG[self.icounter])
				self.counter_low += 1
			
			self.effective_function[self.icounter]=switch
			
			self.icounter += 1
				
		i = i + 1
		
		return i, self.Vexcite       
	
	
	def close(self):
		self.f.close()
		self.f_field.close()   

	def writeout(self):
		label = "./timeseries_D%.2f-%.2f.dat"%(self.noisevariance_u0,self.noisevariance_u1)
		f = open(label,"w+")
		for k in range(self.num_time):
			t = self.dt*float(k)
			str = '%f   %f %f %f\n'%(t,self.u[k],self.Ue[k+self.initialtime],self.u[k]+self.Ue[k+self.initialtime])
			f.write(str)
		f.close()
		print("written %s"%label)

		label = "./timeseries_low.dat"
		duration = len(self.signal_low)
		f = open(label,"w+")
		for k in range(duration):
			t = self.dt*self.sample_factor*float(k)
			str = '%f   %f \n'%(t,self.signal_low[k])
			f.write(str)
		f.close()
		print("written %s"%label)
		
		df = 0.1
		T=1.0/df
		nfft=int(T/self.dt_s)
		if nfft<duration:
			[power0_n,freqs_n]=mlab.psd(self.signal_low, NFFT=nfft,Fs=self.fs_s, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(nfft*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
			num_freqs_n=np.shape(freqs_n)[0]
			df_n=freqs_n[2]-freqs_n[1]
			fmin_n=   2.0
			fmax_n = 30.0
			fmin_num_n=int(fmin_n/df_n)
			fmax_num_n=int(fmax_n/df_n)
			#plt.plot(freqs_n[fmin_num_n:fmax_num_n],power0_n[fmin_num_n:fmax_num_n])
			#plt.show()

		
	def powerspectrum(self):
		## ********************* power spectrum
		deltaf = 0.1
		T=1.0/deltaf
		num_perseg_s = int(T/self.dt_s)
		num_perseg   = int(T/self.dt)

		
		#assume initial state is up
		state = self.expcond[0]
		if state != 1:
			print("ERROR !!!!!!! initial state is not the UP-state")
			quit()
		state_up = []
		state_down = []
		initial = 0
		for i in range(1,self.num_time_s):
			diff = self.expcond[i-1]-self.expcond[i] 
			if diff == 1 : # change of state from UP to DOWN
				final = i
				array = [initial, final]
				state_up.append(array)
				state = 0
				initial = i+1 				
			if diff == -1 : # change of state from UP to DOWN
				final = i
				array = [initial, final]
				state_down.append(array)
				state = 1
				initial = i+1 
		if state == 1:
			final = self.num_time_s
			array = [initial, final]
			state_up.append(array)
		else:
			final = self.num_time_s
			array = [initial, final]
			state_down.append(array)	 
		print("UP state:",state_up) 
		print("DOWN state:",state_down) 

		self.list_length_up = len(state_up)
		self.power_up = []
		self.freqs_up = []
		label = "./result_power_segments_up.dat"
		f = open(label,"w+")
		for k in range(self.list_length_up):
			initial   = (state_up[k])[0]
			final     = (state_up[k])[1]
			deltaf = 0.1
			T=1.0/deltaf
			nfft=int(T/self.dt_s)
			diff = final-initial
			print("UP-state #%d (%d)    nfft=%d   final-initial=%d"%(k,self.list_length_up,nfft,diff))
			if nfft>diff:
				print(" ERROR !!!!!!! in UP-state:   nfft=%d window length=%d ."%(nfft,diff))
				#quit()
				continue
			#nfft      = final-initial
			EEG_array = self.EEG[initial:final] 		
			[power0_n,freqs_n]=mlab.psd(EEG_array, NFFT=nfft,Fs=self.fs_s, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(nfft*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
			num_freqs_n=np.shape(freqs_n)[0]
			print("freqs_n.shape:",np.shape(freqs_n))
			df_n=freqs_n[2]-freqs_n[1]
			fmin_n=2.0
			fmax_n = 50.0
			fmin_num_n=int(fmin_n/df_n)
			fmax_num_n=int(fmax_n/df_n)
			num_n=fmax_num_n-fmin_num_n
			self.freqs_up.append(freqs_n[fmin_num_n:fmax_num_n])
			self.power_up.append(power0_n[fmin_num_n:fmax_num_n])
			str = ""
			for kk in range(num_n):
				str += "%f "%freqs_n[fmin_num_n+kk]
			for kk in range(num_n):
				str += "%f "%power0_n[fmin_num_n+kk]
			str += "\n"
			f.write(str)
		f.close()
		
		self.list_length_down = len(state_down)
		self.power_down = []
		self.freqs_down = []
		label = "./result_power_segments_down.dat"
		f = open(label,"w+")
		for k in range(self.list_length_down):
			initial   = (state_down[k])[0]
			final     = (state_down[k])[1]
			T=1.0/deltaf
			nfft=int(T/self.dt_s)
			print("DOWN-state #%d    nfft=%d"%(k,nfft))
			if nfft>final-initial:
				print(" ERROR !!!!!!! in DOWN_state : nfft=%d window length=%d ."%(nfft,final-initial))
				#quit()
				continue
			EEG_array = self.EEG[initial:final] 		
			[power0_n,freqs_n]=mlab.psd(EEG_array, NFFT=nfft,Fs=self.fs_s, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(nfft*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
			num_freqs_n=np.shape(freqs_n)[0]
			df_n=freqs_n[2]-freqs_n[1]
			fmin_n=2.0
			fmax_n = 50.0
			fmin_num_n=int(fmin_n/df_n)
			fmax_num_n=int(fmax_n/df_n)
			num_n=fmax_num_n-fmin_num_n
			#plt.plot(freqs_n[fmin_num_n:fmax_num_n],power0_n[fmin_num_n:fmax_num_n])
			#plt.show()
			self.freqs_down.append(freqs_n[fmin_num_n:fmax_num_n])
			self.power_down.append(power0_n[fmin_num_n:fmax_num_n])
			str = ""
			for kk in range(num_n):
				str += "%f "%freqs_n[fmin_num_n+kk]
			for kk in range(num_n):
				str += "%f "%power0_n[fmin_num_n+kk]
			str += "\n"
			f.write(str)
		f.close()

		## Welch PSD
		#[power0_n,freqs_n]=mlab.psd(self.u, NFFT=num_perseg,Fs=fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=num_perseg*0.8, pad_to=None, sides='onesided', scale_by_freq=True)
		[power0_n,freqs_n]=mlab.psd(self.u, NFFT=num_perseg,Fs=self.fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(num_perseg*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
		[power0_CT_n,freqs_CT_n]=mlab.psd(self.Ue[self.initialtime:self.initialtime+self.num_time], NFFT=num_perseg,Fs=self.fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(num_perseg*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
		[power0_EEG_n,freqs_CT_n]=mlab.psd(self.EEG, NFFT=num_perseg_s,Fs=self.fs_s, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(num_perseg_s*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
		num_freqs_n=np.shape(freqs_n)[0]
		df_n=freqs_n[2]-freqs_n[1]
		fmin_n=2.0
		#fmin_n=freqs_n[0]
		fmax_n = 50.0
		#fmax_n=freqs_n[num_freqs_n-1]
		fmin_num_n=int(fmin_n/df_n)
		fmax_num_n=int(fmax_n/df_n)
		#print("fmin=%f fmax=%f df=%f num_freq=%d"%(fmin_n,fmax_n,df_n,num_freqs_n))
		#print("fmin_num=%d fmax_num=%d "%(fmin_num_n,fmax_num_n))
		
		self.freqs_p_n=freqs_n[fmin_num_n:fmax_num_n]
		self.p_n=power0_n[fmin_num_n:fmax_num_n]
		self.p_CT_n=power0_CT_n[fmin_num_n:fmax_num_n]
		self.p_EEG_n=power0_EEG_n[fmin_num_n:fmax_num_n]
		#print("frequencies: ",self.freqs_p_n)
		#print("power values cortex: ",self.p_EEG_n)


		## Multitaper PSD
#		power0_n__ = pmtm(data_all0, NW=2.0, show=False)
#		power0_n=power0_n__[0:n/2]
#		freqs_ = scipy.fftpack.fftfreq(n, dt)
#		df_n=1/(dt*n)
#		fmin_n=df_n
#		fmin_num_n=int(fmin_n/df_n)-1
#		fmax_num_n=int(fmax_n/df_n)-1
#		freqs_n = df_n*np.arange(fmin_n,fmax_n,df_n)
#
#		fignn=plt.figure()
#		ax = fignn.add_subplot(411)
#		freqs_p_n=freqs_n[fmin_num_n:fmax_num_n]
#		#print freqs_p_n
#		p_n=power0_n[fmin_num_n:fmax_num_n]
#		plt.plot(freqs_p_n,p_n)

		label = "./ps_D%.3f.dat"%self.noisevariance_u0
		f = open(label,"w+")
		length = np.size(self.freqs_p_n)
		for k in range(length):
			str = '%f %f\n'%(self.freqs_p_n[k],self.p_EEG_n[k])
			f.write(str)
		f.close()
		print("written %s"%label)
		
	def bandpassfilter(self):
		
		def butter_bandpass(lowcut, highcut, fs, order=5):
			nyq = 0.5 * fs
			low = lowcut / nyq
			high = highcut / nyq
			print("low:%f high:%f   lowcut=%f highcut=%f fs=%f"%(low,high,lowcut,highcut,fs))
			b, a = butter(order, [low, high], btype='band')
			return b, a

		self.lowcut_alpha = 6.0#0.05#0.1
		self.highcut_alpha = 12.0#0.7#4.0
		self.lowcut_gamma = 25.0#0.05#0.1
		self.highcut_gamma = 45.0#0.7#4.0
		self.lowcut_full = 1.0#0.05#0.1
		self.highcut_full = 60.0#0.7#4.0
			
		b, a = butter_bandpass(self.lowcut_alpha, self.highcut_alpha, self.fs_s, order=4)
		#zi = lfilter_zi(b, a)
		#self.EEG_alphaband,zf = lfilter(b, a, self.EEG,zi=zi*self.EEG[0])
		self.EEG_alphaband = lfilter(b, a, self.EEG)
		#self.EEG_alphaband = filtfilt(b, a, self.EEG)
		print("self.EEG_alphaband:",self.EEG_alphaband)

		b, a = butter_bandpass(self.lowcut_gamma, self.highcut_gamma, self.fs_s, order=4)
		self.EEG_gammaband = filtfilt(b, a, self.EEG)

		b, a = butter_bandpass(self.lowcut_full, self.highcut_full, self.fs_s, order=4)
		self.EEG_fullband = filtfilt(b, a, self.EEG)
		
		#### power spectra of bandpass-filtered EEG
		deltaf = 0.1
		T=1.0/deltaf
		num_perseg_s=int(T/self.dt_s)
		
		## Welch PSD
		[power_gammaband,freqs_n]=mlab.psd(self.EEG_gammaband, NFFT=num_perseg_s,Fs=self.fs_s, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(num_perseg_s*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
		num_freqs_n=np.shape(freqs_n)[0]
		df_n=freqs_n[2]-freqs_n[1]
		fmin_n=self.lowcut_gamma
		#fmin_n=freqs_n[0]
		fmax_n=self.highcut_gamma
		#fmax_n=freqs_n[num_freqs_n-1]
		fmin_num_n=int(fmin_n/df_n)
		fmax_num_n=int(fmax_n/df_n)
		#print("fmin=%f fmax=%f df=%f num_freq=%d"%(fmin_n,fmax_n,df_n,num_freqs_n))
		#print("fmin_num=%d fmax_num=%d "%(fmin_num_n,fmax_num_n))
		self.freqs_p_gammaband=freqs_n[fmin_num_n:fmax_num_n]
		self.p_gammaband=power_gammaband[fmin_num_n:fmax_num_n]

		[power_alphaband,freqs_n]=mlab.psd(self.EEG_alphaband, NFFT=num_perseg_s,Fs=self.fs_s, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(num_perseg_s*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
		fmin_n=self.lowcut_alpha
		#fmin_n=freqs_n[0]
		fmax_n=self.highcut_alpha
		#fmax_n=freqs_n[num_freqs_n-1]
		fmin_num_n=int(fmin_n/df_n)
		fmax_num_n=int(fmax_n/df_n)
		#print("fmin=%f fmax=%f df=%f num_freq=%d"%(fmin_n,fmax_n,df_n,num_freqs_n))
		#print("fmin_num=%d fmax_num=%d "%(fmin_num_n,fmax_num_n))
		self.freqs_p_alphaband=freqs_n[fmin_num_n:fmax_num_n]
		self.p_alphaband=power_alphaband[fmin_num_n:fmax_num_n]

		[power_fullband,freqs_n]=mlab.psd(self.EEG_fullband, NFFT=num_perseg_s,Fs=self.fs_s, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=int(num_perseg_s*0.8), pad_to=None, sides='onesided', scale_by_freq=True)		
		fmin_n=self.lowcut_full
		fmax_n=self.highcut_full
		fmin_num_n=int(fmin_n/df_n)
		fmax_num_n=int(fmax_n/df_n)
		self.freqs_p_fullband=freqs_n[fmin_num_n:fmax_num_n]
		self.p_fullband=power_fullband[fmin_num_n:fmax_num_n]


		### instantaenous power of bandpass-filtered data
		self.instpower_window = int(self.shortest_interval_s*0.1)
		duration = self.num_time_s-self.instpower_window
		self.instpower_gamma  = np.zeros(duration)
		self.instpower_alpha  = np.zeros(duration)
		#print("length of window:%d  duration=%d"%(self.instpower_window,duration))
		for i in range(duration):
			sum_gamma = 0.0
			sum_alpha = 0.0
			for j in range(self.instpower_window):
				t = i+j
				sum_gamma += self.EEG_gammaband[t]**2
				sum_alpha += self.EEG_alphaband[t]**2
			sum_gamma /= np.sqrt(sum_gamma/float(self.instpower_window))
			sum_alpha /= np.sqrt(sum_alpha/float(self.instpower_window))
			self.instpower_gamma[i]=sum_gamma
			self.instpower_alpha[i]=sum_alpha
			if i%1000==0:
				print("i=%d (%d) sum_gamma=%f sum_alpha=%f"%(i,duration,sum_gamma,sum_alpha))
		#print("instpower_gamma:",self.instpower_gamma)
				