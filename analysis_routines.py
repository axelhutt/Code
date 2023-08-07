#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

class tfmethods():
    
    def forceAspect(self,ax,aspect=1):
        im = ax.get_images()
        extent =  im[0].get_extent()
        ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    
    def Wavelettransform_Morlet(self,y):
        from scipy.fft import fft, ifft
        
        def FT_realmotherwavelet(k):
            nu = k*self.df
            return np.sqrt(np.pi)*np.exp(-0.25*(self.sigma-2.0*np.pi*nu)**2)
        
        def FT_complexmotherwavelet(k):
            nu = k*self.df
            return np.sqrt(np.pi)*(np.exp(-0.25*(self.sigma-2.0*np.pi*nu)**2)-np.exp(-0.25*self.sigma**2)*np.exp(-np.pi**2*nu**2))
        
        def FT_motherwavelet(k):
            #return FT_realmotherwavelet(k)
            return FT_complexmotherwavelet(k)
        
        num_time = np.size(y)
        
        self.sigma = 10.0
        center_frequency = self.sigma/(2.0*np.pi)
        self.df    = (self.fmax-self.fmin)/float(self.num_scales)#1.0/self.T
        
        #self.num_scales = 20 # number of scales
        self.freqs = np.linspace(self.fmin,self.fmax,self.num_scales)
        a = center_frequency/(self.freqs)
        DFT_convolved = np.zeros(num_time,dtype=complex)
        WT = np.zeros((self.num_scales,num_time),dtype=complex)
        
        #print("all prepared......")
        DFT = fft(y)
        #print("size of y:",np.shape(y))
        #print("size of DFT:",np.shape(DFT))
        nums2 = self.num_scales//2
        for l in range(self.num_scales):
            #if l%nums2 == 0:
            #    print("l=%d (%d)"%(l,self.num_scales))
            for k in range(num_time//2):
                DFT_convolved[k] = DFT[k]*FT_motherwavelet(a[l]*k)
                kk = k+num_time//2
                DFT_convolved[kk] = DFT[kk]*FT_motherwavelet(a[l]*(-num_time//2+k))
            WT[l,:] = ifft(DFT_convolved)
            #print("scale %d WT:"%l,WT[l,0:10]," DFT:",DFT[0:10])
        return WT	
    
 
    def Compute_phaseinference_band(self,y1,y2,fmin,fmax,window1,window2,chan_list):
        #import numpy.ma as ma
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        ## define time step
        dt=self.dt
        self.fs = 1.0/self.dt 
        
        ntimes1 = int(window1/self.dt)
        ntimes2 = int(window2/self.dt)
        #ntimes=self.num_time
        
        if ntimes1<ntimes2:
            ntimes2=ntimes1
        else:
            ntimes1=ntimes2
        ntimes = ntimes1
        
        data0_1 = np.zeros((ntimes1,self.num_chan))
        data0_2 = np.zeros((ntimes2,self.num_chan))
        data0_1[:,:]=y1[:ntimes1,:] 
        data0_2[:,:]=y2[:ntimes2,:] 
        del y1
        del y2
        endtime1=ntimes1*dt
        endtime2=ntimes2*dt
        
        #data = np.zeros((ntimes,self.num_chan))
        anasig1 = np.zeros((ntimes1,self.num_chan),dtype=complex)
        anasig2 = np.zeros((ntimes2,self.num_chan),dtype=complex)
        from scipy.signal import hilbert
        for k in range(self.num_chan):
            data=self.Compute_filter_butterworth(data0_1[:,k],fmin,fmax)
            anasig1[:,k] = hilbert(data)
            data=self.Compute_filter_butterworth(data0_2[:,k],fmin,fmax)
            anasig2[:,k] = hilbert(data)
        del data0_1
        del data0_2
        #amplitude = np.abs(analytical_signal)
       
        #chan_list = [[0,1],[0,14],[0,15],\
        #    [1,2],[1,3],[1,4],[1,13],[1,14],\
        #    [2,3],[2,4],[2,5],[2,11],[2,13],\
        #    [3,4],[3,5],[3,8],[3,11],[3,13],
        #    [4,5],[4,8],[4,11],\
        #    [5,6],[5,8],[5,11],\
        #    [6,7],[6,8],\
        #    [7,8],[9,11],\
        #    [10,11],\
        #    [11,12],[11,13],\
        #    [12,13],[13,14],[14,15],\
        #    ]           
        import subprocess
        pvalue_list = []
        #for chan1 in range(self.num_chan-1):
        #    chan01 = chan1+1
        #    for chan02 in range(self.num_chan-chan01):
        #        chan2 = chan01+chan02
        #        chan_list.append([chan1,chan2])
        for k in range(len(chan_list)):
            chan1=chan_list[k][0]
            chan2=chan_list[k][1]
            #print("chan01=%d chan02%d   chan1=%d chan2=%d"%(chan01,chan02,chan1,chan2))
            dangle1 = np.angle(anasig1[:,chan1])-np.angle(anasig1[:,chan2])
            dangle2 = np.angle(anasig2[:,chan1])-np.angle(anasig2[:,chan2])
            f = open("data.dat","w")
            for i in range(ntimes1):
                str = '%f %f\n'%(dangle1[i],dangle2[i])
                f.write(str)
            f.close()
            res = subprocess.call("Rscript ./WWTest.R", shell=True)
            f = open("out.dat","r")
            p=f.read()
            #print("p-value:",p)
            pvalue_list.append(float(p))
        
        #alpha = 0.05/len(chan_list)
        #for k in range(len(chan_list)):
        #    if pvalue_list[k]<alpha:
        #        print("%d %d   - p=%g"%(chan_list[k][0],chan_list[k][1],pvalue_list[k]))
        
        return [chan_list,pvalue_list]


    ####### compute p-values in each subset (condition) separately
    def Compute_phaseinference_band_intra(self,y1,fmin,fmax,window1,chan_list):
        #import numpy.ma as ma
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        ## define time step
        dt=self.dt
        self.fs = 1.0/self.dt 
        
        ntimes1 = int(window1/self.dt)
        ntimes = ntimes1
        
        data0_1 = np.zeros((ntimes1,self.num_chan))
        data0_1[:,:]=y1[:ntimes1,:] 
        del y1
        endtime1=ntimes1*dt
        
        #data = np.zeros((ntimes,self.num_chan))
        anasig1 = np.zeros((ntimes1,self.num_chan),dtype=complex)
        from scipy.signal import hilbert
        for k in range(self.num_chan):
            data=self.Compute_filter_butterworth(data0_1[:,k],fmin,fmax)
            anasig1[:,k] = hilbert(data)
        del data0_1
                   
        import subprocess
        pvalue_list = []
        chan_pval_list = []
        for k in range(len(chan_list)):
            chan1=chan_list[k][0]
            chan2=chan_list[k][1]
            dangle1 = np.angle(anasig1[:,chan1])-np.angle(anasig1[:,chan2])
            for kk0 in range(len(chan_list)-k-1):
                kk = k+kk0+1
                chan11=chan_list[kk][0]
                chan22=chan_list[kk][1]
                chan_pval_list.append([k,kk])
                dangle2 = np.angle(anasig1[:,chan11])-np.angle(anasig1[:,chan22])
                print("k=%d  %d %d     kk=%d %d %d "%(k,chan1,chan2,kk,chan11,chan22))
                f = open("data.dat","w")
                for i in range(ntimes1):
                    str = '%f %f\n'%(dangle1[i],dangle2[i])
                    f.write(str)
                f.close()
                res = subprocess.call("Rscript ./WWTest.R", shell=True)
                f = open("out.dat","r")
                p=f.read()
                #print("p-value:",p)
                pvalue_list.append(float(p))
                #print(chan_pval_list)
        #alpha = 0.05/len(chan_list)
        #for k in range(len(chan_list)):
        #    if pvalue_list[k]<alpha:
        #        print("%d %d   - p=%g"%(chan_list[k][0],chan_list[k][1],pvalue_list[k]))
            
        return [chan_pval_list,pvalue_list]
    

    def compute_permutation_test_PLV(self,data1,data2,PLV):

        data = np.concatenate((data1,data2))
        num1 = np.size(data1)
        num2 = np.size(data2)
        if num1 != num2:
            print("different length of data !!! Abort !!!")
            quit()
        data_new1 = np.zeros(np.size(data1))
        data_new2 = np.zeros(np.size(data2))
        del data1
        del data2
        num_perm = 200
        PLV_shuffled = np.zeros(num_perm)
        counter = 0
        print("PLV=%f"%PLV)
        #data_orig = data.copy()
        data_shuffled = data.copy()
        p=0
        counter = 0
        for perm in range(num_perm):
            np.random.shuffle(data_shuffled)
            data_new1 = data_shuffled[0:num1]
            data_new2 = data_shuffled[num1:num1+num2]
            dangle = data_shuffled[0:num1]-data_shuffled[num1:num1+num2]#data_new1-data_new2
            #print("shuffled data:",data_shuffled)
            helpcos = 0.0
            helpsin = 0.0
            for i in range(num1):
                help0 = dangle[i]
                #print("i=%d help0=%f"%(i,help0))
                helpcos += np.cos(help0)/float(num1)
                helpsin += np.sin(help0)/float(num1)
            PLV_shuffled[perm]=np.sqrt(helpcos*helpcos+helpsin*helpsin)
            if PLV<PLV_shuffled[perm] :
                counter+= 1
                p=float(counter)/float(num_perm)
                #print("mean>pe_shuffled[perm]: p=%f"%p)
            print("perm=%d PLV=%f   p=%f"%(perm,PLV_shuffled[perm],p))
        print("final p-value=%f"%p)
        quit()

    def Compute_phasecoherence_band(self,y,fmin,fmax,window,stype):
        import numpy.ma as ma
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        ## define time step
        dt=self.dt
        self.fs = 1.0/self.dt 
        
        ntimes = int(window/self.dt)
        #ntimes=self.num_time
        data0 = np.zeros((ntimes,self.num_chan))
        data0[:,:]=y[:ntimes,:] 
        
        endtime=ntimes*dt
        
        #data = np.zeros((ntimes,self.num_chan))
        anasig = np.zeros((ntimes,self.num_chan),dtype=complex)
        sigamp = np.zeros((ntimes,self.num_chan))
        from scipy.signal import hilbert
        for k in range(self.num_chan):
            data=self.Compute_filter_butterworth(data0[:,k],fmin,fmax)
            anasig[:,k] = hilbert(data)
            sigamp[:,k] = np.abs(anasig[:,k])
            #print(data[10000:])
            #quit()
        del data0
        
        #amplitude = np.abs(analytical_signal)
        
        
        powerthreshold=0.99
        #print("compute wavelet transform.....")
        
        windowduration_time=self.dt#2.0#10.0
        windowduration = int(windowduration_time/self.dt)
        imax_essential = ntimes-windowduration
        if imax_essential <= 0:
            print("time window (%f) is shorter than window duration (%f)! Abort !"%(ntimes,windowduration_time))
            quit()
            
        ## jumps (time jump windowstep*dt) of time windows of duration windowduration_time
        windowstep = windowduration
        num_timewindows = int(imax_essential/windowduration)
        print("windowduration:%d windowstep:%d ntimes:%d   num_timewindows=%d"%(windowduration,windowstep,ntimes,num_timewindows))
        
        results_d = []
        results   = []
        
        num_couplechans = (self.num_chan*(self.num_chan-1))//2
        connectivity_matrix = np.identity(self.num_chan)#np.zeros((self.num_chan,self.num_chan))
        connectivityphase_matrix = np.identity(self.num_chan)*0.0#np.zeros((self.num_chan,self.num_chan))
        PLV = 0.0#np.zeros(windowduration)
        MeanPhase = 0.0#np.zeros(windowduration)
        for chan1 in range(self.num_chan-1):
            chan01 = chan1+1
            for chan02 in range(self.num_chan-chan01):
                chan2 = chan01+chan02
                print("chan01=%d chan02%d   chan1=%d chan2=%d"%(chan01,chan02,chan1,chan2))
        
                PLV_d 		= 0.0
                PLV_d_var 	= 0.0
                
                PLV=0.0
                MeanPhase=0.0
                
                dangle = np.angle(anasig[:,chan1])-np.angle(anasig[:,chan2])
                
                print ('compute PLV.....')
                helpcos = 0.0
                helpsin = 0.0
                for i in range(num_timewindows):
                    #if i%100==0:
                    #    print("i=%d(%d)"%(i,num_timewindows))
                    ii = i*windowstep # means ii=i
                    #print("i=%d ii=%d w=%d ii+w=%d shape=%d"%(i,ii,w,ii+w,np.shape(dangle)[1]))
                    help0 = dangle[ii]
                    helpcos += np.cos(help0)/float(num_timewindows)
                    helpsin += np.sin(help0)/float(num_timewindows)
                    #print("w=%d dangle:"%i,dangle[ii])
                PLV=np.sqrt(helpcos*helpcos+helpsin*helpsin) ## circular variance
                #self.compute_permutation_test_PLV(np.angle(anasig[:,chan1]),np.angle(anasig[:,chan2]),PLV)
                MeanPhase=np.arctan2(helpsin,helpcos) ## mean phase
                MeanPhase_stddev = np.sqrt(-2*np.log(PLV))    
                PLV_d       = PLV
                
#               for w in range(windowduration):
#                   if w%500 == 0:
#                       print("w=%d(%d)"%(w,windowduration))
#                   helpcos=0.0
#                   helpsin=0.0
#                   for i in range(num_timewindows):
#                       #if i%100==0:
#                       #    print("i=%d(%d)"%(i,num_timewindows))
#                       ii = i*windowstep
#                       #print("i=%d ii=%d w=%d ii+w=%d shape=%d"%(i,ii,w,ii+w,np.shape(dangle)[1]))
#                       help0 = dangle[ii+w]
#                       helpcos += np.cos(help0)/float(num_timewindows)
#                       helpsin += np.sin(help0)/float(num_timewindows)
#                       #print("w=%d dangle:"%w,dangle[:,ii+w])
#                   PLV[w]=np.sqrt(helpcos*helpcos+helpsin*helpsin) ## circular variance
#                   MeanPhase[w]=np.arctan2(helpsin,helpcos) ## mean phase
#                   
#                   
#                   PLV_d       += PLV[w]/float(windowduration)
#                   PLV_d_var	+= PLV[w]**2/float(windowduration)
#                   hcos        += np.cos(MeanPhase[w])/float(windowduration)
#                   hsin        += np.sin(MeanPhase[w])/float(windowduration)
#
#               MeanPhase_avg = np.arctan2(hsin,hcos)
#               R = np.sqrt(hcos**2+hsin**2)
#               MeanPhase_avg_stddev = np.sqrt(-2*np.log(R))
#               PLV_d_var = PLV_d_var - PLV_d**2
#               
                #z = windowduration*(R**2)
                z = windowduration*(PLV**2)
                p = np.exp(-z)
                
                print("cond:%s   chan1:%d chan2:%d   gPLV=%f    MeanPhase=%f(%f)  #trials:%d   p:%f"%(\
                        stype,\
                        chan1,chan2,\
                        PLV,\
                        MeanPhase,MeanPhase_stddev,
                        num_timewindows,p\
                ))
                
                
                connectivity_matrix[chan1,chan2]=PLV
                connectivity_matrix[chan2,chan1]=PLV
                connectivityphase_matrix[chan1,chan2]=MeanPhase
                connectivityphase_matrix[chan2,chan1]=-MeanPhase
                results_d.append([PLV,MeanPhase,MeanPhase_stddev,num_timewindows])
                
        results.append(connectivity_matrix)
        results.append(connectivityphase_matrix)
        results.append(results_d)
        results.append(num_timewindows)
        results.append(sigamp)
        
        return results
    
    def Compute_filter_butterworth(self,signal,lowcut,highcut,order=3):
        from scipy.signal import butter, lfilter
        def butter_bandpass(lowcut, highcut, fs, order=5):
            return butter(order, [lowcut, highcut], fs=fs, btype='band')
        
        def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
            b, a = butter_bandpass(lowcut, highcut, fs, order=order)
            y = lfilter(b, a, data)
            return y
        if lowcut==0.1 or lowcut==0.3:
            signal_filtered = butter_bandpass_filter(signal, lowcut, highcut, self.fs, order=3)
        else:
            signal_filtered = butter_bandpass_filter(signal, lowcut, highcut, self.fs, order=4)
        
        return signal_filtered
    
    def Compute_phasecoherence(self,y1,y2,):
        import numpy.ma as ma
        #import matplotlib.pyplot as plt
        #import scipy.fftpack
        #import matplotlib.mlab as mlab
        #from cwt_modules_Bergner import *
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        ## define tiem step
        dt=self.dt
        
        # now downsampling 
        counter = 0
        counter_sample = 0
        #print ('downsample data...')
        ntimes=self.num_time
        data=np.zeros((ntimes,2))
        data[:,0]=y1.copy()
        data[:,1]=y2.copy()
        dataavg = np.mean(data,axis=1)
        
        #print ('\n perform time-frequency analysis....')
        
        timemax=ntimes*dt
        endtime=timemax

        self.fmin=0.5
        self.fmax=30.0#fs/2.0 # Nyquist
        self.num_scales = 80
        dff=(self.fmax-self.fmin)/float(self.num_scales)
        
        powerthreshold=0.99
        #print("compute wavelet transform.....")
        
        res1 = self.Wavelettransform_Morlet(data[:,0])
        res2 = self.Wavelettransform_Morlet(data[:,1])
        dangle = np.angle(res1)-np.angle(res2)
        power = (np.abs(res1)+np.abs(res2))/2
        
        poweravg=np.abs(self.Wavelettransform_Morlet(dataavg))
        maxpower_freq=np.squeeze(np.max(power,axis=0))
        minpower_freq=np.squeeze(np.min(power,axis=0))
        diffpower_freq=maxpower_freq-minpower_freq
        
        ##### plot result
        #print ('compute PLV.....')
        
        PLV_d 		= 0.0
        PLV_d_var 	= 0.0
        f_d_min 		= 0.5
        f_d_max 		= 4.0
        
        PLV_t 		= 0.0
        PLV_t_var 	= 0.0
        f_t_min 		= 4.0
        f_t_max 		= 8.0
        
        PLV_a 		= 0.0
        PLV_a_var 	= 0.0
        f_a_min 		= 8.0
        f_a_max 		= 12.0
        
        PLV_b 		= 0.0
        PLV_b_var 	= 0.0
        f_b_min 		= 12.0
        f_b_max 		= 25.0
        
        PLV_g 		= 0.0
        PLV_g_var 	= 0.0
        f_g_min 		= 25.0
        f_g_max 		= 29.0
        
        num_freqs = np.shape(dangle)[0]
        
        num_f = int((f_d_max-f_d_min)/dff)
        k0 = int(f_d_min/dff)
        #print("d: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_t_max-f_t_min)/dff)
        k0 = int(f_t_min/dff)
        #print("t: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_a_max-f_a_min)/dff)
        k0 = int(f_a_min/dff)
        #print("a: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_b_max-f_b_min)/dff)
        k0 = int(f_b_min/dff)
        #print("b: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_g_max-f_g_min)/dff)
        k0 = int(f_g_min/dff)
        #print("g: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        
        windowduration=100
        windowstep = 1 ### !=1: to be implemented
        imax_essential = (ntimes-windowduration)//windowstep
        #print ntimes-windowduration
        PLV = np.zeros((num_freqs,imax_essential))
        helpcos=np.zeros(num_freqs)
        helpsin=np.zeros(num_freqs)
        imax2 = imax_essential//5
        for i in range(imax_essential):
            #if i%imax2 == 0:
            #    print("%d(%d)"%(i,imax_essential))
            ii = i*windowstep
            helpcos[:]=0.0
            helpsin[:]=0.0
            for w in range(windowduration):
                    help = np.squeeze(dangle[:,ii+w])
                    helpcos=helpcos+np.cos(help)/(windowduration)
                    helpsin=helpsin+np.sin(help)/(windowduration)
            help=np.sqrt(helpcos*helpcos+helpsin*helpsin)
            ### threshold PLV by average power
            if self.FLAG_THR == 1:
                help1=ma.masked_where((np.squeeze(power[:,i])-minpower_freq[i])/diffpower_freq[i]<powerthreshold,help) 
                PLV[:,i]=help1.filled(0)
            else:
                PLV[:,i]=help
                
            # delta
            num_f = int((f_d_max-f_d_min)/dff)
            k0 = int(f_d_min/dff)
            for k in range(num_f):
                PLV_d		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_d_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # theta
            num_f = int((f_t_max-f_t_min)/dff)
            k0 = int(f_t_min/dff)
            for k in range(num_f):
                PLV_t		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_t_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # alpha
            num_f = int((f_a_max-f_a_min)/dff)
            k0 = int(f_a_min/dff)
            for k in range(num_f):
                PLV_a		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_a_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # beta
            num_f = int((f_b_max-f_b_min)/dff)
            k0 = int(f_b_min/dff)
            for k in range(num_f):
                PLV_b		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_b_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # gamma
            num_f = int((f_g_max-f_g_min)/dff)
            k0 = int(f_g_min/dff)
            for k in range(num_f):
                PLV_g		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_g_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
                
        PLV_d_var = PLV_d_var - PLV_d**2
        PLV_t_var = PLV_t_var - PLV_t**2
        PLV_a_var = PLV_a_var - PLV_a**2
        PLV_b_var = PLV_b_var - PLV_b**2
        PLV_g_var = PLV_g_var - PLV_g**2
        #print("gPLV: delta=%f(%f) theta=%f(%f) alpha=%f(%f) beta=%f(%f) gamma=%f(%f)  number of values=%d"%(\
        #        PLV_d,np.sqrt(PLV_d_var),\
        #        PLV_t,np.sqrt(PLV_t_var),\
        #        PLV_a,np.sqrt(PLV_a_var),\
        #        PLV_b,np.sqrt(PLV_b_var),\
        #        PLV_g,np.sqrt(PLV_g_var),
        #        num_f*(ntimes-windowduration)\
        #    ))
        
        results = []
        results.append([PLV_d,np.sqrt(PLV_d_var),\
            PLV_t,np.sqrt(PLV_t_var),\
            PLV_a,np.sqrt(PLV_a_var),\
            PLV_b,np.sqrt(PLV_b_var),\
            PLV_g,np.sqrt(PLV_g_var),\
            num_f*(ntimes-windowduration)])
        results.append(poweravg[:,0:ntimes-windowduration])
        results.append(PLV)
        
#       fign=plt.figure()
#       ax = fign.add_subplot(211)
#       plt.imshow(poweravg[:,0:ntimes-windowduration], interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
#       plt.colorbar()
#       forceAspect(ax,aspect=2)
#       
#       ax = fign.add_subplot(212)
#       plt.imshow(PLV, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-windowduration*dt, self.fmin, self.fmax])
#       #plt.imshow(PLV, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-windowduration*dt, fmin, fmax])
#       forceAspect(ax,aspect=2)
#       plt.colorbar()
#       
#       plt.show()
#       
        return results
    
    def compute_phases(self,y,fmin,fmax,timewindow0,stype,chanlist):
        #import numpy.ma as ma
        ## define time step
        dt=self.dt
        self.fs = 1.0/self.dt 
        
        ntimes = int(timewindow0/self.dt)
        data0 = np.zeros((ntimes,self.num_chan))
        data0[:,:]=y[:ntimes,:] 
        endtime=ntimes*dt
        del y
        
        chan_index = np.zeros(self.num_chan,dtype=int)
        for c_in_list in chanlist:
            chan_index[c_in_list[0]]=1
            chan_index[c_in_list[1]]=1
        
        angle= np.zeros((ntimes,self.num_chan))
        from scipy.signal import hilbert
        for k in range(self.num_chan):
            if chan_index[k]==1:
                data=self.Compute_filter_butterworth(data0[:,k],fmin,fmax)
                anasig = hilbert(data)
                #angle[:,k] = np.random.uniform(0.0,6.28,len(data))#np.angle(anasig)
                angle[:,k] = np.angle(anasig)
        del data0
        
        phase_list = []
        for chans in chanlist:
            dphase = angle[:,chans[0]]-angle[:,chans[1]]
            #print(chans[0],chans[1],dphase[:3])
            phase_list.append(dphase)
        
        return phase_list
                
    
    def Compute_WindowedDFT(self,y,fmin,fmax,deltaf,num_windows):
        
        num_time = len(y)
        #num_windows =10
        num_onewindow = num_time//num_windows 
        num_perseg = num_onewindow
        fs = 1.0/self.dt
        
        #df_n = fs/float(num_perseg)
        #print("--------- df_n=%f"%df_n)
        #fmin_num=int(fmin/df_n)
        #fmax_num=int(fmax/df_n)
        #num_freqs_n = int((fmax-fmin)/df_n)
        
        #print("shape of data y:",np.shape(y))
        #power_tf=np.zeros((num_windows,num_freqs_n))
        
        power_list = []
        for num in range(num_windows):
            #print("num %d(%d)"%(num,num_windows))
            y_window = y[num*num_onewindow:(num+1)*num_onewindow-1]
            # Welch PSD
            [power,freqs]=mlab.psd(y_window, NFFT=num_perseg,Fs=fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=num_perseg*0.9, pad_to=None, sides='onesided', scale_by_freq=True)
            power_list.append(power)
            #power_tf[num,:] = power[fmin_num:fmax_num]
        #freqs_n = freqs_n[fmin_num:fmax_num]
        
        num_freqs=np.shape(freqs)[0]
        df_n=freqs[2]-freqs[1]
        #print("num_freqs:%d df_n=%f"%(num_freqs,df_n))
        #fmin_n = 0.05
        #fmax_n = 70.0
        if freqs[-1]<=fmax:
            fmax=freqs[-1]
            print("fmax adapted to range of possible frequencies !")
        fmin_num_n=int(fmin/df_n)
        fmax_num_n=int(fmax/df_n)
        num_freqs_n = (fmax_num_n-fmin_num_n)
        #print("fmin:",fmin_num_n," fmax:",fmax_num_n," df:",df_n," num_freqs:",num_freqs_n)
        #num_freqs_n = (fmax_num_n-fmin_num_n)
        power_tf=np.zeros((num_windows,num_freqs_n))
        freqs_n  = freqs[fmin_num_n:fmax_num_n]
        #print(":",fmin_num_n,fmax_num_n)
        #print(freqs)
        #print(freqs_n)
        #print("shape(freqs):",np.shape(freqs_n))
        #p_true     = power_mean[fmin_u_num_n:fmax_u_num_n]
        for num in range(num_windows):
            p = power_list[num]
            #print(np.shape(p))
            power_tf[num,:]=p[fmin_num_n:fmax_num_n]
        #print("shape(freqs):",np.shape(freqs_n))
        
        results = [] 
        results.append(freqs_n)
        results.append(power_tf)
        return results
    

    def Compute_powerspectrum_band_freqs(self,y,fmin,fmax):
        deltaf = 0.1
        T=10.0/deltaf
        fs = 1.0/self.dt
        num_perseg   = 1*int(T/self.dt)
        
        num_time = len(y)
        num_windows =10
        num_onewindow = num_time//num_windows 
        power = 0.0
        df_n = fs/float(num_perseg)
        # Welch PSD
        [power,freqs]=mlab.psd(y, NFFT=num_perseg,Fs=fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=num_perseg*0.9, pad_to=None, sides='onesided', scale_by_freq=True)
        
        res = []
        fmin_num=int(fmin/df_n)
        fmax_num=int(fmax/df_n)
        freqs_ = freqs[fmin_num:fmax_num]
        power_ = power[fmin_num:fmax_num]
        res.append(freqs_)
        res.append(power_)
        return res
    

    def Compute_powerspectrum_band(self,y,fmin,fmax):
        deltaf = 0.1
        T=10.0/deltaf
        fs = 1.0/self.dt
        num_perseg   = 1*int(T/self.dt)
        
        num_time = len(y)
        num_windows =10
        num_onewindow = num_time//num_windows 
        power = 0.0
        df_n = fs/float(num_perseg)
        # Welch PSD
        [power,freqs]=mlab.psd(y, NFFT=num_perseg,Fs=fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=num_perseg*0.9, pad_to=None, sides='onesided', scale_by_freq=True)
        
        fmin_num=int(fmin/df_n)
        fmax_num=int(fmax/df_n)
        power_avg = np.mean(power[fmin_num:fmax_num])
        return power_avg
    
    def Compute_powerspectrum(self,y):
        deltaf = 0.1
        T=10.0/deltaf
        fs = 1.0/self.dt
        num_perseg   = 1*int(T/self.dt)
        
        num_time = len(y)
        num_windows =10
        num_onewindow = num_time//num_windows 
        fmin_d = 0.1
        fmax_d = 4.0
        fmin_t = 4.0
        fmax_t = 8.0
        fmin_a = 8.0
        fmax_a = 12.0
        fmin_b = 12.0
        fmax_b = 25.0
        power_d     = 0
        power_t     = 0
        power_a     = 0
        power_b     = 0
        power_d_var = 0.0
        power_t_var = 0.0
        power_a_var = 0.0
        power_b_var = 0.0
        
        power_mean = 0.0
        power_var = 0.0
        df_n = fs/float(num_perseg)
        #print("--------- df_n=%f"%df_n)
        
        #print("shape of data y:",np.shape(y))
        
        for num in range(num_windows):
            #print("num %d(%d)"%(num,num_windows))
            y_window = y[num*num_onewindow:(num+1)*num_onewindow-1]
            # Welch PSD
            [power,freqs_n]=mlab.psd(y, NFFT=num_perseg,Fs=fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=num_perseg*0.9, pad_to=None, sides='onesided', scale_by_freq=True)
            
            fmin_num=int(fmin_d/df_n)
            fmax_num=int(fmax_d/df_n)
            power_d += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_d_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            fmin_num=int(fmin_t/df_n)
            fmax_num=int(fmax_t/df_n)
            power_t += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_t_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            fmin_num=int(fmin_a/df_n)
            fmax_num=int(fmax_a/df_n)
            power_a += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_a_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            fmin_num=int(fmin_b/df_n)
            fmax_num=int(fmax_b/df_n)
            power_b += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_b_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            power_mean += power/float(num_windows)
            power_var  += np.power(power,2)/float(num_windows)
            
        power_d_var = power_d_var-power_d**2
        power_t_var = power_t_var-power_t**2
        power_a_var = power_a_var-power_a**2
        power_b_var = power_b_var-power_b**2
        power_var = power_var-np.power(power_mean,2)
        
        #print("power_d=%d(%g)  power_t=%d(%g)  power_a=%d(%g)  power_b=%d(%g)  "%(\
        #        power_d,power_d_var,\
        #        power_t,power_t_var,\
        #        power_a,power_a_var,\
        #        power_b,power_b_var))  
        
        num_freqs_n=np.shape(freqs_n)[0]
        df_n=freqs_n[2]-freqs_n[1]
        
        fmin_n = 0.05
        fmax_n = 70.0
        fmin_u_num_n=int(fmin_n/df_n)
        fmax_u_num_n=int(fmax_n/df_n)
        freqs_p_n  = freqs_n[fmin_u_num_n:fmax_u_num_n]
        p_true     = power_mean[fmin_u_num_n:fmax_u_num_n]
        var_true   = power_var[fmin_u_num_n:fmax_u_num_n]
        
        results = [] 
        r1 = [power_d,power_d_var,power_t,power_t_var,power_a,power_a_var,power_b,power_b_var,var_true]
        results.append(r1)
        results.append(freqs_p_n)
        results.append(p_true)
        return results
    
    def clustering_KDE(self,data0):# one-dimensional data
        
                
        #################################################
        ########### kernel density method ###############
        #################################################
        from sklearn.neighbors import KernelDensity
        #print(data0)
        #print(np.shape(data0))
        #print(type(data0))
        #data = data0[:,np.newaxis]
        data = np.zeros((len(data0),1))
        data[:,0]=data0
        num_kernels = 10
        bw_min= 0.001
        bw_max = 0.1
        kde_list = []
        #for nk in range(num_kernels):
        #bw = bw_min+float(nk)*
        #kde = KernelDensity(kernel="gaussian", bandwidth=0.05).fit(data)
        #kde = KernelDensity(kernel="tophat", bandwidth=0.05).fit(data)
        kde = KernelDensity(kernel="epanechnikov", bandwidth=0.1).fit(data)
        #kde = KernelDensity(kernel="tophat", bandwidth=0.005).fit(data)
        #data_uniform0 = np.linspace(np.min(data0),np.max(data0),len(data0))
        data_uniform0 = np.linspace(0.0,1.0,len(data0))
        data_uniform = data_uniform0[:,np.newaxis]
        log_dens = kde.score_samples(data)
        log_dens_uniform = kde.score_samples(data_uniform)
        return data_uniform,np.exp(log_dens),np.exp(log_dens_uniform)


    def clustering_phase_KDE(self,data0,bandwidth):# one-dimensional data
        
        
        #################################################
        ########### kernel density method ###############
        #################################################
        from sklearn.neighbors import KernelDensity
        #print(data0)
        #print(np.shape(data0))
        #print(type(data0))
        #data = data0[:,np.newaxis]
        data = np.zeros((len(data0),1))
        data[:,0]=data0
        kde = KernelDensity(kernel="cosine", bandwidth=bandwidth).fit(data)
        #data_uniform0 = np.linspace(np.min(data0),np.max(data0),len(data0))
        data_uniform0 = np.linspace(-np.pi,np.pi,len(data0))
        data_uniform = data_uniform0[:,np.newaxis]
        log_dens = kde.score_samples(data)
        log_dens_uniform = kde.score_samples(data_uniform)
        return data_uniform,np.exp(log_dens),np.exp(log_dens_uniform)


    def clustering_hierarchical(self,data0,num_cluster):
        
        #################################################
        ########### hierarchical clustering method ###############
        #################################################
        from sklearn.cluster import AgglomerativeClustering
        from scipy.cluster.hierarchy import dendrogram, linkage
        
        ##########################
        ######## artificial data #
        ###########################
        #num_simuldata = 100
        #data_simul = np.zeros(num_simuldata)
        #data_simul[10:num_simuldata//5] = np.random.uniform(0.55,0.65,num_simuldata//5-10)
        #data_simul[num_simuldata//5:] = np.random.uniform(0.0,0.1,4*num_simuldata//5)
        #data_simul[:10] = np.random.uniform(0.95,1.0,10)
        #data0 = data_simul
        
        
        
        data = np.zeros((len(data0),1))
        data[:,0]=data0
        data_random = data.copy()
        
        if num_cluster>0:
            flag_numclusters_auto = 0
        else:
            flag_numclusters_auto = 1
            
        data_copy = data.copy()
        
        if flag_numclusters_auto == 1:
        
            def compute_gapstatistic(labels0,num_cluster,num_pairs):
                W = 0.0
                cluster_channel_list = []
                for c in range(num_cluster):
                    cluster_channel_list.append([])
                D = np.zeros(num_cluster)
                for chan in range(num_pairs):
                    cluster_channel_list[labels0[chan]].append(chan)
                
                for c in range(num_cluster):
                    ll_list = cluster_channel_list[c]
                    #print("c=%d  ll_list:"%c,ll_list)
                    numchans = len(ll_list)
                    #for e1 in range(numchans):
                    #    for e20 in range(numchans-e1-1):
                    for e1 in range(numchans):
                        for e2 in range(numchans):
                            #e2 = e1+e20+1
                            D[c]+=(data0[ll_list[e1]]-data0[ll_list[e2]])**2
                            #print("%d %d D=%f"%(ll_list[e1],ll_list[e2],D[c]))
                            
                    W+=D[c]/float(2*numchans)
                if W>0:
                    logW=np.log(W)
                else:
                    logW=-1000000.0
                return logW
            
            num_pairs = len(data0)
            
            ### compute uniformly distributed data
            max_ = np.max(data0)
            min_ = np.min(data0)
            num_samples = 40
            data0_uniform = np.random.uniform(min_,max_,(num_samples,num_pairs))
            
            num_cluster_max = 40
            if num_cluster_max>num_pairs:
                num_cluster_max=num_pairs
                
            W  = np.zeros(num_cluster_max)
            EW = np.zeros(num_cluster_max)
            s  = np.zeros(num_cluster_max)
            Gap= np.zeros(num_cluster_max)
            Gap_stat= np.zeros(num_cluster_max)
            Gapdiff = np.zeros(num_cluster_max)
            for nc in range(num_cluster_max):
                num_cluster = nc+1
                hierarchical_cluster = AgglomerativeClustering(n_clusters=num_cluster, affinity='euclidean', linkage='ward')
                labels0 = hierarchical_cluster.fit_predict(data)
                W[nc]=compute_gapstatistic(labels0,num_cluster,num_pairs)
#               cluster_channel_list = []
#               for c in range(num_cluster):
#                   cluster_channel_list.append([])
#               D = np.zeros(num_cluster)
#               for chan in range(num_pairs):
#                   cluster_channel_list[labels0[chan]].append(chan)
#               W[nc]=0.0
#               for c in range(num_cluster):
#                   ll_list = cluster_channel_list[c]
#                   numchans = len(ll_list)
#                   print(ll_list)
#                   for e1 in range(numchans):
#                       for e20 in range(numchans-e1-1):
#                           e2 = e1+e20+1
#                           #print("%d %d"%(ll_list[e1],ll_list[e2]))
#                           D[c]+=(data0[ll_list[e1]]-data0[ll_list[e2]])**2
#                   W[nc]+=D[c]/float(2*numchans)
#               W[nc]=np.log(W[nc])
                help = np.zeros(num_samples)
                for k in range(num_samples):
                    data_random[:,0]=data0_uniform[k,:]
                    hierarchical_cluster = AgglomerativeClustering(n_clusters=num_cluster, affinity='euclidean', linkage='ward')
                    labels0 = hierarchical_cluster.fit_predict(data_random)
                    help[k] = compute_gapstatistic(labels0,num_cluster,num_pairs)
                    EW[nc]+=help[k]/float(num_samples)
                EWstddev = 0.0
                for k in range(num_samples):
                    EWstddev+= (help[k]-EW[nc])**2/float(num_samples)
                EWstddev = np.sqrt(EWstddev)
                
                s[nc] = EWstddev*np.sqrt(1.0+1.0/float(num_samples))
                Gap[nc] = EW[nc]-W[nc] 
                Gap_stat[nc]=Gap[nc]-s[nc]
                if nc>2:
                    Gapdiff[nc] = Gap[nc-1]>Gap_stat[nc]
                    #print("nc=%d  EW=%f W=%f  Gap=%f Gap_stat=%f  Gapdiff=%f  %f"%(nc,EW[nc],W[nc],Gap[nc],Gap[nc]-s[nc],Gapdiff[nc],s[nc]))
            for nc in range(num_cluster_max):
                if nc>2 and Gapdiff[nc]>0:
                    print("optimal number of clusters is :",nc-1)
                    break    
            
            #print("data comparison:",data_copy-data)
            
            num_cluster = nc-1
            #hierarchical_cluster_n = AgglomerativeClustering(n_clusters=num_cluster, affinity='euclidean', linkage='ward')
            hierarchical_cluster_n = AgglomerativeClustering(n_clusters=num_cluster, affinity='euclidean', linkage='single')
            labels0_n = hierarchical_cluster_n.fit_predict(data)
            labels = labels0_n+1
            
            #fig2 = plt.figure(2)
            ##plt.plot(range(num_cluster_max),Gapdiff,'ko-',range(num_cluster_max),EW,'r',range(num_cluster_max),W,'b')
            #plt.plot(range(num_cluster_max),Gapdiff,'ko-')
            #plt.show()
#           quit()
            
            #print("cluster labels:",labels)
            
        if flag_numclusters_auto == 0:
            #print("Hierarchical clustering with pre-defined number of cluster : %d"%num_cluster)            
#           p = 100
#           n = len(data0)
#           linkage_data = linkage(data, method='ward', metric='euclidean')
#           print("linkage_data:",linkage_data[n-p-2:,:])
#           dendrogram(linkage_data,p=p,truncate_mode='lastp')
#           plt.show()
#           quit()
            
            #distance_threshold = 0.02#0.1#0.08
            #hierarchical_cluster = AgglomerativeClustering(n_clusters=None, distance_threshold=distance_threshold, affinity='euclidean', linkage='ward',compute_full_tree=True)
            #labels0 = hierarchical_cluster.fit_predict(data)
            #labels = labels0+1
            #num_cluster = 1+np.amax(hierarchical_cluster.labels_)
            #print(f"Number of clusters = {1+np.amax(hierarchical_cluster.labels_)}")
            #num_cluster = 3
            
            hierarchical_cluster = AgglomerativeClustering(n_clusters=num_cluster, affinity='euclidean', linkage='ward')
            labels0 = hierarchical_cluster.fit_predict(data)
            labels = labels0+1
            #print("cluster labels:",labels)
            
        #print("number of clusters:",num_cluster)
        return num_cluster,labels
        