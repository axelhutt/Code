#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pylt

import analysis_routines as ar

class data(ar.tfmethods):
    def __init__(self):
        self.num_times = 0
        self.dt_orig = 0.001
        self.filelabel_connectivity="./connectivities+phases"
            
        self.MWTest_conds_list = []
        for k in range(3):
            self.MWTest_conds_list.append([])
        
        self.projections_list = []
        for k in range(3):
            self.projections_list.append([])
        
        self.symmetry_list = []
        for l in range(2):# for different fmin
            self.symmetry_list.append([])
            for k in range(3):
                self.symmetry_list[l].append([])
                
        self.mutual_connections_data_list = []
        self.mutual_connections_not_data_list = []
        for k in range(3):
            self.mutual_connections_data_list.append([])
            self.mutual_connections_not_data_list.append([])
        
        self.connect_notmutual123_list = []
        for k in range(3):
            self.connect_notmutual123_list.append([])

        self.mutual_connections_data_list
        self.mutual_connections_topphases_allssesions_low = np.zeros((16,16))
        self.mutual_connections_allssesions_list = []
        for k in range(3):
            self.mutual_connections_allssesions_list.append(np.zeros((16,16)))
        self.mutual_connections_allanaconds_list = []
        for k in range(7):
            self.mutual_connections_allanaconds_list.append(np.zeros((16,16)))
            #self.mutual_connections_allssesions_low = 
        self.connection_mutual_list = []
        self.connection_mutual_sig_list = []
        self.connection_mutual_transitions_list = []
        self.topcluster_meanphase1_singlesession_list = []
        self.topcluster_meanphase2_singlesession_list = []
        self.topcluster_meanphase3_singlesession_list = []
        self.topcluster_stdphase1_singlesession_list = []
        self.topcluster_stdphase2_singlesession_list = []
        self.topcluster_stdphase3_singlesession_list = []
        for k in range(3):
            self.connection_mutual_list.append([])
            self.connection_mutual_sig_list.append([])
            self.connection_mutual_transitions_list.append([])
            self.topcluster_meanphase1_singlesession_list.append([])
            self.topcluster_stdphase1_singlesession_list.append([])
            self.topcluster_meanphase2_singlesession_list.append([])
            self.topcluster_stdphase2_singlesession_list.append([])
            self.topcluster_meanphase3_singlesession_list.append([])
            self.topcluster_stdphase3_singlesession_list.append([])
        self.topcluster_meanphase1_sessionpool_list = []
        self.topcluster_meanphase2_sessionpool_list = []
        self.topcluster_meanphase3_sessionpool_list = []
        
    def read_data(self):
        
        ######### low anaesthesia
        path = './'
        #dir0='20110810'
        #dir=dir0+'B/'
        #label = '141236' #pre-3.6%
        #label = '142258' #stim-3.6%
        #label = '143322' #post-3.6%
        
        #dir0='20110817'
        #dir=dir0+'A/'
        #label = '100125' #pre-3.4%
        #label = '101232'  #stim-3.4%
        #label = '102258' #post-3.4%
        
        #dir0='20110817'
        #dir=dir0+'B/'
        #label = '143429' #pre-3.2%
        #label = '145103' #stim-3.2%
        #label = '150346' #post-3.2%  beta+gamma peak
        
        #dir0='20110824'
        #dir=dir0+'A/'
        #label = '095341' #pre-3.4%
        #label = '100405' #stim-3.4%
        #label = '101435' #post-3.4%
        
        #dir0='20110824'
        #dir=dir0+'B/'
        #label = '142058' #pre-3.8% gamma+beta peak
        #label = '143128' #stim-3.8%
        #label = '144201' #post-3.8%
        
        #dir0='20110913'
        #dir=dir0+'A/'
        #label = '091038' #pre-3.4%
        #label = '092105' #stim-3.4%
        #label = '093134'#post-3.4%
        
        #dir0='20110913'
        #dir=dir0+'B/'
        #label = '124242' #pre-3.6% 
        #label = '125305' #stim-3.6%
        #label = '130326' #post-3.6%
        
        
        ######## medium anaesthesia
        #path = './'
        #dir0='20110810'
        #dir=dir0+'B/'
        #label = '130846' #pre-4.5%
        #label = '131915' #stim-4.5%
        #label = '133032' #post-4.5%
        
        #dir0='20110817'
        #dir=dir0+'A/'
        #label = '105605' #pre-4.4%
        #label = '110638' #stim-4.4% beta+gamma
        #label = '111710' #post-4.4% beta+gamma
        
        #dir0='20110817'
        #dir=dir0+'B/'
        #label = '152321' #pre-4.4%
        #label = '153347' #stim-4.4%
        #label = '154410' #post-4.4% gamma
        
        #dir0='20110824'
        #dir=dir0+'A/'
        #label = '105219' #pre-4.6% gamma
        #label = '110311' #stim-4.6% gamma
        #label = '111336' #post-4.6% gamma
        
        #dir0='20110824'
        #dir=dir0+'B/'
        #label = '151704' #pre-5.0%
        #label = '155128' #stim-5.0% beta+gamma
        #label = '160150' #post-5.0%  beta+gamma
        
        #dir0='20110913'
        #dir=dir0+'A/'
        #label = '095555' #pre-4.6% beta+gamma
        #label = '101810' #stim-4.6% gamma
        #label = '102834' #post-4.6% gamma
        
        #dir0='20110913'
        #dir=dir0+'B/'
        #label = '132936' #pre-4.3% gamma
        #label = '134008' #stim-4.3%
        #label = '135024' #post-4.3%
        
        
        ######## high anaesthesia
        #path = './'
        #dir0='20110810'
        #dir=dir0+'B/'
        #label = '121713' #pre-6.0%
        #label = '123022' #stim-6.0%
        #label = '124107' #post-6.0%
        
        #dir0='20110817'
        #dir=dir0+'A/'
        #label = '114515' #pre-6.0%
        #label = '115543' #stim-6.0% gamma
        #label = '120605' #post-6.0%
        
        #dir0='20110817'
        #dir=dir0+'B/'
        #label = '133358' #pre-6.0%
        #label = '134428' #stim-6.0% beta+gamma
        #label = '135458' #post-6.0% beta+gamma
        
        #dir0='20110824'
        #dir=dir0+'A/'
        #label = '115341' #pre-6.1% gamma
        #label = '120418' #stim-6.1% beta+gamma
        #label = '121733' #post-6.1% beta+gamma
        
        #dir0='20110824'
        #dir=dir0+'B/'
        #label = '132740' #pre-6.1%
        #label = '133808' #stim-6.1%
        #label = '134831' #post-6.1%
        
        #dir0='20110913'
        #dir=dir0+'A/'
        #label = '105534' #pre-6.0%
        #label = '110634' #stim-6.0%
        #label = '111651' #post-6.0% gamma
        
        dir0='20110913'
        dir=dir0+'B/'
        #label = '115353' #pre-6.0% gamma
        #label = '120413' #stim-6.0%
        label = '121435' #post-6.0%
        
        infile = path+dir+dir0+'-'+label+'-001.dat'
        self.outfile = path+dir+dir0+'-'+label+'_artred.dat' ## artifact reduced data
        self.data_orig=np.loadtxt(infile)
        
        

    def read_data_artred(self,dir0,dir1,label,analevel):
        
        path = './'
        outpath = './analysis/'+analevel+'/'
                
        dir=dir0+dir1
        label_file_in = path+dir+dir0+'-'+label+'_artred'
        label_file_out = outpath+dir+dir0+'-'+label+'_artred'
        infile = label_file_in+'.dat' ## artifact reduced data
        outfile = label_file_out+'.dat' ## artifact reduced data
        self.outfile = label_file+'power.dat' ## artifact reduced data
        self.outfile_bands = label_file+'powerbands.dat' ## artifact reduced data
        
        data0 = np.loadtxt(infile)
        self.data = np.zeros((np.shape(data0)[0],np.shape(data0)[1]-1))
        self.data=data0[:,1:]
        self.num_time = np.shape(self.data)[0]
        self.num_chan = np.shape(self.data)[1]
        self.dt = data0[1,0]-data0[0,0]
        self.times = data0[:,0]
        del data0
        
    def compute_write_powerspectra(self):
        bands_list = []
        power_list = []
        for kk in range(self.num_chan):
            res = self.Compute_powerspectrum(self.data[kk,:])
            bands_list.append(res[0])            
            if kk==0:
                power_f = res[1]
            power_p = res[2]
            power_list.append(power_p)
            
        ### write out power in bands
        fpowerbands = open(self.outfile_bands,"w+")
        num_bands = 4
        for k in range(num_bands):
            str_pb = ''
            for kk in range(self.num_chan):
                str_pb += '\t%f'%(bands_list[kk])[k*2]
            str_pb += '\n'
            fpowerbands.write(str_pb)
        fpowerbands.close()
        print("written %s"%self.outfile_bands)

        ###### write out power spectra
        fpower = open(self.outfile,"w+")
        lmax = len(power_f)
        for l in range(len(power_f)):
            if power_f[l]>50.0:
                lmax=l-1
                break
        for freq in range(lmax):
            str_p = ''
            for kk in range(self.num_chan):
                str_p += '\t%f'%(power_list[kk])[l]
            str_p += '\n'
            fpower.write(str_p)
        fpower.close()
        print("written %s"%self.outfile)    
            
            
    def analyze_data(self):
        num_times = np.shape(self.data_orig)[0]
        self.channel_indices = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]  #could be channel_indexes = [0]
        self.num_chan = len(self.channel_indices)
        if np.shape(self.data_orig)[1] != self.num_chan+1:
            print("%d != %d !!!! Quit !!!"%(np.shape(self.data_orig)[1],self.num_chan))
            quit()
        # now downsampling 
        sample_factor=5
        self.num_time=int(num_times/sample_factor)
        self.signal = np.zeros(self.num_time)
        self.data = np.zeros((self.num_chan,self.num_time))
        for k in range(self.num_chan):
            for i in range(self.num_time):
                self.data[k,i]=self.data_orig[i*sample_factor,self.channel_indices[k]]
        for i in range(self.num_time):
            self.signal[i]=self.data_orig[i*sample_factor,16]
        self.dt=self.dt_orig*sample_factor
        endtime = self.num_time*self.dt
        self.times = np.linspace(0.0,endtime,self.num_time)
        self.times_stored = self.times.copy()
        print("shape:",np.shape(self.data))
        del self.data_orig
        
        ## remove stimulation artifacts
        flag_stim = 0
        counter_stim = 0
        stimulations = []
        for i in range(self.num_time):
            if self.signal[i]>10000.0 and flag_stim == 0: ## start of stimulation
                flag_stim = 1
                counter_stim += 1
                stim_start = i-200
                print("stimulation start (%d): %f"%(counter_stim,stim_start*self.dt))
            if self.signal[i]<-10000.0 and flag_stim == 1: ## end of stimulation
                flag_stim = 0
                stim_end = i+1000
                stimulations.append([stim_start,stim_end])
                print("stimulation end (%d): %f"%(counter_stim,stim_end*self.dt))
            #if flag_stim == 0:
            #    data_final_list.append(self.data[:,i])
        data_final_list = []
        stimulation_counter = 0
        flag_stim = 0
        
        if len(stimulations)>0:
            print("stimulation found !!!!!!!!!!!")
            for i in range(self.num_time):
                if i==(stimulations[stimulation_counter])[0] and flag_stim==0 :
                    flag_stim = 1
                    print("data point at stim start# %d (t=%f):%f"%(stimulation_counter,i*self.dt,self.data[0,i]))
                if i==(stimulations[stimulation_counter])[1] and flag_stim==1 :
                    flag_stim = 0
                    print("data point at stim end# %d (t=%f):%f"%(stimulation_counter,i*self.dt,self.data[0,i]))
                    stimulation_counter +=1
                if flag_stim == 0:
                    data_final_list.append(self.data[:,i])
                if stimulation_counter == 10:
                    break
        else:
            print("stimulation NOT found !!!!!!!!!!!")
            for i in range(self.num_time):
                data_final_list.append(self.data[:,i])
                
        self.data_stored = self.data.copy()
        del self.data
        del self.times 
        self.num_time = len(data_final_list)
        endtime = self.num_time*self.dt
        self.times = np.linspace(0.0,endtime,self.num_time)
        self.data = np.zeros((self.num_chan,self.num_time))
        for i in range(self.num_time):
            self.data[:,i]=data_final_list[i]
        print("new length of data:",self.num_time)
        
        ## final artifact reduction step which is independent of AINP channel
        print("maximum:",np.max(self.data,axis=1))
        mean = np.mean(self.data,axis=1)
        stddev = np.sqrt(np.var(self.data,axis=1))
        print("mean:",mean,"  stddev:",stddev)
        for i in range(self.num_time):
            for k in range(self.num_chan):
                #if np.abs(self.data[k,i])>20000.0:
                #    print("data[%d,%d]=%f 10*stddev=%f"%(k,i*self.dt,self.data[k,i],10*stddev[k]))
                if np.abs(self.data[k,i])>5*stddev[k]:
                    self.data[k,i]=mean[k] 
                    #print("outlier!!!!!!!!!!!!! at time %f in channel %d"%(i*self.dt,k))
        print("maximum:",np.max(self.data,axis=1))
    
    def remove_powerlinepeak(self,path,ana,stype_list,file0):
        import h5py
        
        num_stype = 3
        counter = 0
        counter_all = 0
        label_list = []
        power_list = []
        label_list = []
        freqs_list = []
        for stype in stype_list:
            print(stype)
            #label = path+file0+'_'+ana+'_'+stype+'.h5'
            label = path+file0+'_'+ana+'_'+stype+'_fs500.0_'+'.h5'
            print("datafile:",label)
            label_list.append(stype)
            data=self.read_data_h5(label) #data=[times,chans]
            data_out = np.zeros((self.num_time,self.num_chan+1))
            data_out[:,1:] = data.copy()
            data_out[:,0] = self.time
            for k in range(self.num_chan):
                data_out[:,1+k]=self.notch_filter(1.0,60.0,data[:,k])
        
            #labelout = path+file0+'_'+ana+'_'+stype+'_notch'+'.h5'
            labelout = path+file0+'_'+ana+'_'+stype+'_fs500.0_'+'notch'+'.h5'
            h5f = h5py.File(labelout, 'w')
            label_data = 'LFP_anaesthesia_level:'+ana+'_'+stype
            ret = h5f.create_dataset(label_data, data=data_out)
            print(ret)
            h5f.close()
            print("converted %s to %s"%(label,labelout))
        
            power_spectra = []
            for kk in range(self.num_chan):
                print("channel %d(%d)"%(kk,self.num_chan))
                res = self.Compute_powerspectrum(data[:,kk])
                res_new = self.Compute_powerspectrum(data_out[:,1+kk])
                if kk == 0:
                    freqs = res[1]
                    freqs_new = res_new[1]
                power_rel = res[2]
                power_rel_new = res_new[2]
                power_spectra.append([power_rel,power_rel_new])
            freqs_list.append([freqs,freqs_new])    
            power_list.append(power_spectra)
            #counter += 1
            #counter_all += 1
        
        print("num of freqs:",len(freqs_list[0]))
        fig1 = plt.figure(file0,figsize=(10,7))
        for chan in range(self.num_chan):
            chan0 = chan+1
            ax = plt.subplot(4,4,chan0)
            plt.plot(freqs_list[0][1],np.log(power_list[0][chan][1][:]),color='k',label='pre')            
            plt.plot(freqs_list[1][1],np.log(power_list[1][chan][1][:]),color='r',label='stim')
            plt.plot(freqs_list[2][1],np.log(power_list[2][chan][1][:]),color='b',label='post')
            if chan==0:
                ax.legend(loc=1,fontsize=8)
        fig1.subplots_adjust(hspace=0.3,wspace=0.3)
        plt.show()
        
        
        
        
    def notch_filter(self,bw,f0,data0):
        ## notch filter at 60Hz
        import scipy.signal as scisig
        #bw = 1.0
        Q = f0/bw
        fs = 1.0/self.dt
        b, a = scisig.iirnotch(f0, Q, fs=fs)
        data = data0.copy() 
        data = scisig.filtfilt(b, a, data0)
        return data
    
    def read_data_h5(self,label):
        import h5py
        h5f = h5py.File(label, 'r')
        datasetkey = list(h5f.keys())
        #print("key:",datasetkey)
        dataset = h5f[datasetkey[0]]
        self.num_time = dataset.shape[0]
        self.num_chan = dataset.shape[1]-1
        #print("time points:%d  channels:%d"%(self.num_time,self.num_chan))
        dataset0 = h5f[datasetkey[0]][()]
        dataset = np.zeros((self.num_time,self.num_chan))
        dataset = dataset0[:,1:]
        self.time = dataset0[:,0]
        self.dt = self.time[1]-self.time[0]
        #print("time :",self.time)
        h5f.close()
        #quit()
        return dataset
    
    def average_powerspectra(self,path,ana,stype,file_list):
        import h5py
        import statistics as stat 
        import scipy.stats as ss
        
        def read_data(label):
            h5f = h5py.File(label, 'r')
            datasetkey = list(h5f.keys())
            print("key:",datasetkey)
            dataset = h5f[datasetkey[0]]
            self.num_time = dataset.shape[0]
            self.num_chan = dataset.shape[1]-1
            print("time points:%d  channels:%d"%(self.num_time,self.num_chan))
            dataset0 = h5f[datasetkey[0]][()]
            dataset = np.zeros((self.num_time,self.num_chan))
            dataset = dataset0[:,1:]
            self.time = dataset0[:,0]
            self.dt = self.time[1]-self.time[0]
            h5f.close()
            return dataset
        
        num_files = len(file_list)
        counter = 0
        counter_all = 0
        power_list = []
        bands_list = []
        label_list = []
        fs = 500.0
        for file in file_list:
            print(file)
            if fs==200.0:
                label = path+file+'_'+ana+'_'+stype+'_notch.h5'
            if fs==500.0:
                label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
            print("datafile:",label)
            label_list.append(file)
            data=read_data(label)
            
            power_bands = []
            power_spectra = []
            for kk in range(self.num_chan):
                res = self.Compute_powerspectrum(data[kk,:])
                if counter == 0 and kk == 0:
                    freqs = res[1]
                poweravg = 1.0#np.mean(res[2])
                power_rel = res[2]/poweravg
                power_bands.append(res[0])
                power_spectra.append(power_rel)
                
            power_list.append(power_spectra)
            bands_list.append(power_bands)
            counter += 1
            counter_all += 1
            
        num_files_new = num_files    
        num_freqs = len(freqs)
        power_all = np.zeros((self.num_chan,num_freqs,num_files_new))
        for k in range(num_files_new):
            p = power_list[k]
            for chan in range(self.num_chan):
                power_all[chan,:,k] = p[chan]
        shape = np.shape(power_all)
        power_median = np.zeros((shape[0],shape[1]))
        power_sem = np.zeros((shape[0],shape[1]))
        for k in range(shape[0]):
            for l in range(shape[1]):
                power_median[k,l] = stat.median(power_all[k,l,:])
                power_sem[k,l] = ss.sem(power_all[k,l,:])
        power_avg = np.mean(power_all,axis=2)
        
#       plt.figure(1)
#       chan = 0
#       ax = plt.subplot(1,1,1)
#       plt.plot(freqs,power_avg[chan,:],'k')
#       for k in range(num_files_new):
#           p = power_list[k]
#           p_ = p[chan]
#           plt.plot(freqs,p_[:],'--',label=label_list[k])
#       ax.legend()
#       plt.show()
        
        return freqs,power_avg[:,:],power_median,power_sem

    def read_data(self,label):
        import h5py
        h5f = h5py.File(label, 'r')
        datasetkey = list(h5f.keys())
        print("key:",datasetkey)
        dataset = h5f[datasetkey[0]]
        self.num_time = dataset.shape[0]
        self.num_chan = dataset.shape[1]-1
        print("time points:%d  channels:%d"%(self.num_time,self.num_chan))
        dataset0 = h5f[datasetkey[0]][()]
        dataset = np.zeros((self.num_time,self.num_chan))
        dataset = dataset0[:,1:]
        self.time = dataset0[:,0]
        self.dt = self.time[1]-self.time[0]
        h5f.close()
        return dataset
    
    def compute_permutation_test_pe(self,data,mean):
        import ordpy
        num_perm = 200
        pe_shuffled = np.zeros(num_perm)
        counter = 0
        print("mean pe=%f"%mean)
        data_orig = data.copy()
        data_shuffled = data.copy()
        p=-1
        for perm in range(num_perm):
            np.random.shuffle(data_shuffled)
            #print("shuffled data:",data_shuffled)
            pe_shuffled[perm] = ordpy.permutation_entropy(data_shuffled,dx=3,taux=1,tauy=1,normalized=True)
            if mean>pe_shuffled[perm] :
                counter+= 1
                p=float(counter)/float(num_perm)
                print("mean>pe_shuffled[perm]: p=%f"%p)
            print("%d %f"%(perm,pe_shuffled[perm]))
        print("final p-value=%f"%p)
        quit()
        
    def compute_PE(self,ana_list,file_list,fs):
        import ordpy
        stype_list = ['pre','stim','post']
        
        num_chan = 16
        #fs = 500.0
        dt = 1.0/fs
        timewindow=520.0
        num_trials = 1
        
        final_results = []
        
        
        for ana_ in range(len(ana_list)):
            ana = ana_list[ana_]
            path = './analysis/'+ana+'/'
            final_results.append([])
            for file0 in file_list:
                #file='20110913B'
                
                PE = np.zeros((3,16))
                PE_var = np.zeros((3,16))
                counter = 0
                for stype in stype_list:
                    #print(stype)
                    #label = path+file0+'_'+ana+'_'+stype+'_notch.h5'
                    if fs==200.0:
                        label = path+file0+'_'+ana+'_'+stype+'_notch.h5'
                    if fs==500.0:
                        label = path+file0+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
                    #print("datafile:",label)
                    #label_list.append(stype)
                    data=self.read_data_h5(label)
                    num_time = np.shape(data)[0]
                    #data_list.append(data)
                    #print("shape of data:",np.shape(data))
                    if num_time<(timewindow/dt):
                        print("num_time<(timewindow/dt): %d < %d    : choose duration %d"%(num_time,int(timewindow/dt),num_time))
                        timewindow0=num_time*dt
                    else:
                        timewindow0 = timewindow
                    #timewindow_list.append(timewindow0)
                    
                    num_trial = int(timewindow0/(dt*num_trials))
                    num_timewindow0 = int(timewindow0/dt)
                    #print("num_trial=%d"%num_trial)
                    pe_ = np.zeros(num_trials)
                    for k in range(num_chan):
                        mean = 0.0
                        var = 0.0
                        loop_counter = 0
                        for trial in range(num_trials):
                            if (trial+1)*num_trial>int(timewindow/self.dt):
                                break
                            datatrial = data[trial*num_trial:(trial+1)*num_trial,k]
                            pe_[trial] = ordpy.permutation_entropy(datatrial,dx=3,taux=1,tauy=1,normalized=True)
                            mean += pe_[trial]
                            loop_counter += 1
                            #print("trial %d     %f %g        %d - %d "%(trial,pe,np.sqrt(var),trial*num_trial,(trial+1)*num_trial))
                        mean /=float(loop_counter)
                        print("mean=%f"%mean)
                        #self.compute_permutation_test_pe(data[:num_timewindow0,k],mean)
                        for trial in range(loop_counter):
                            var += (pe_[trial]-mean)**2/float(loop_counter)
                        #print("mean=%f var=%f  loop_counter=%d"%(mean,var,loop_counter))
                        PE[counter,k] = mean
                        PE_var[counter,k] = var
                        print("%s %s %s  PE[%d,%d]=%f (%f)"
                                %(ana,file0,stype,counter,k,PE[counter,k],np.sqrt(PE_var[counter,k])))
                    counter += 1
                final_results[ana_].append(PE)
                
                label0 = '_%.1f_'%fs
                if fs==200.0:
                    label='PE'+label0+file0+'_'+ana+'.dat'
                if fs==500.0:
                    label='PE'+label0+file0+'_'+ana+'_fs500.0'+'.dat'
                f = open(label,"w+")
                for chan in range(num_chan):
                    str = '%d  '%chan
                    for k in range(3):
                        str += '  %f  '%PE[k,chan]
                    str += '\n'
                    print(str)
                    f.write(str)
                f.close()
                print("written %s"%label)
        
        
    def compute_PE_statistics(self,ana_list,file_list,fs):
        import statistics as stat
        import scipy.stats
        from matplotlib import ticker
        
        stype_list = ['pre','stim','post']
        label0 = '_%.1f_'%fs
        
        PE_chanpoolmedian_list = []
        sign_chanpoollist_list = []
        PE_chanmedian_list = []
        PE_chanpool_list = []
        PE_all_list = []
        for k in range(len(ana_list)):
            PE_all_list.append([])
            PE_chanmedian_list.append([])
            PE_chanpool_list.append([])
            for kk in range(3):
                PE_chanpool_list[k].append([])
                
        alpha = 0.05
        maxy = np.zeros(len(file_list))-1000000000.0
        miny = np.zeros(len(file_list))+100000000.0
        counter = 1
        
        self.num_chan = 0
        
        sign_diffslist_list = []
        for ana_ in range(len(ana_list)):
            ana = ana_list[ana_]
            path = './analysis/'+ana+'/'
            signlist_list = []
            sign_diffslist_list.append([])
            for file0_ in range(len(file_list)):
                file0 = file_list[file0_]
                if fs==200.0:
                    label='PE'+label0+file0+'_'+ana+'.dat'
                if fs==500.0:
                    label='PE'+label0+file0+'_'+ana+'_fs500.0'+'.dat'
                data_in = np.loadtxt(label)
                self.num_chan = np.shape(data_in)[0]
                PE_all_list[ana_].append(data_in)
                
                mean_diff12 = 0.0
                mean_diff23 = 0.0
                mean_diff13 = 0.0
                collect12 = []
                collect23 = []
                collect13 = []
                for k in range(np.shape(data_in)[0]):
                    diff12 = data_in[k,1]-data_in[k,2]
                    diff23 = data_in[k,2]-data_in[k,3]
                    diff13 = data_in[k,1]-data_in[k,3]
                    collect12.append(diff12)
                    collect23.append(diff23)
                    collect13.append(diff13)
                    #print("diff12=%f  %f %f"%(diff12,diff23,diff13))
                    PE_chanpool_list[ana_][0].append(diff12)
                    PE_chanpool_list[ana_][1].append(diff23)
                    PE_chanpool_list[ana_][2].append(diff13)
                res12 = scipy.stats.wilcoxon(collect12, y=None, zero_method='wilcox', correction=False, alternative='two-sided')
                res23 = scipy.stats.wilcoxon(collect23, y=None, zero_method='wilcox', correction=False, alternative='two-sided')
                res13 = scipy.stats.wilcoxon(collect13, y=None, zero_method='wilcox', correction=False, alternative='two-sided')
                res12_23 = scipy.stats.mannwhitneyu(collect12, collect23)
                res23_13 = scipy.stats.mannwhitneyu(collect23, collect13)
                res12_13 = scipy.stats.mannwhitneyu(collect12, collect13)
                median12 = stat.median(collect12)
                median23 = stat.median(collect23)
                median13 = stat.median(collect13)
                #print("median12=%f  %f %f"%(median12,median23,median13))
                meds = [median12,median23,median13]
                PE_chanmedian_list[ana_].append(meds)
                
                sign_diffs_list = []
                if res12[1]<alpha:
                    #print("ana:%s  12  medians: %f    p=%f"%(ana,median12,res12[1]))
                    sign_diffs_list.append(1)
                if res23[1]<alpha:
                    #print("ana:%s  23  medians: %f    p=%f"%(ana,median23,res23[1]))
                    sign_diffs_list.append(2)
                if res13[1]<alpha:
                    #print("ana:%s  13  medians: %f    p=%f"%(ana,median13,res13[1]))
                    sign_diffs_list.append(3)
                sign_diffslist_list[ana_].append(sign_diffs_list)
                    
                sign_list = []
                if res12_23[1]<alpha:
                    #print("ana:%s  12-23  medians: %f  %f   p=%f"%(ana,median12,median23,res12_23[1]))
                    sign_list.append([1,2])
                if res23_13[1]<alpha:
                    #print("ana:%s  23-13  medians: %f  %f   p=%f"%(ana,median23,median13,res23_13[1]))
                    sign_list.append([2,3])
                if res12_13[1]<alpha:
                    #print("ana:%s  12-13  medians: %f  %f   p=%f"%(ana,median12,median13,res12_13[1]))
                    sign_list.append([1,3])
                signlist_list.append(sign_list)
                
                if maxy[file0_] < max(meds):
                    maxy[file0_] = max(meds) 
                if miny[file0_] > min(meds):
                    miny[file0_] = min(meds) 
                
            median1 = stat.median(PE_chanpool_list[ana_][0])
            median2 = stat.median(PE_chanpool_list[ana_][1])
            median3 = stat.median(PE_chanpool_list[ana_][2])
            PE_chanpoolmedian_list.append([median1,median2,median3])
            res1 = scipy.stats.wilcoxon(PE_chanpool_list[ana_][0], y=None, zero_method='wilcox', correction=False, alternative='two-sided')
            res2 = scipy.stats.wilcoxon(PE_chanpool_list[ana_][1], y=None, zero_method='wilcox', correction=False, alternative='two-sided')
            res3 = scipy.stats.wilcoxon(PE_chanpool_list[ana_][2], y=None, zero_method='wilcox', correction=False, alternative='two-sided')
            sign_chanpool_list = []
            if res1[1]<alpha:
                #print("ana:%s  1  median: %f    p=%f"%(ana,median1,res1[1]))
                sign_chanpool_list.append(1)
            else:
                print("ana:%s  1     p=%f"%(ana,res1[1]))
            if res2[1]<alpha:
                #print("ana:%s  2  median: %f    p=%f"%(ana,median2,res2[1]))
                sign_chanpool_list.append(2)
            else:
                print("ana:%s  2     p=%f"%(ana,res2[1]))
            if res3[1]<alpha:
                #print("ana:%s  3  median: %f    p=%f"%(ana,median3,res3[1]))
                sign_chanpool_list.append(3)
            else:
                print("ana:%s  3     p=%f"%(ana,res3[1]))
            sign_chanpoollist_list.append(sign_chanpool_list)
        
        all_cluster1_mean = []
        all_clusters1_channel = []
        all_clusters1_pe = []
        all_cluster2_mean = []
        all_clusters2_channel = []
        all_clusters2_pe = []
        all_cluster3_mean = []
        all_clusters3_channel = []
        all_clusters3_pe = []
        
        nobasecluster1_avganimal_pe_list = []
        nobasecluster2_avganimal_pe_list = []
        nobasecluster3_avganimal_pe_list = []
        
        numclusters=3
        counter = 1
        plt.figure("PE of single subjects")
        for ana_ in range(len(ana_list)):
            ana = ana_list[ana_]
            path = './analysis/'+ana+'/'
            all_cluster1_mean.append([])
            all_clusters1_channel.append([])
            all_clusters1_pe.append([])
            all_cluster2_mean.append([])
            all_clusters2_channel.append([])
            all_clusters2_pe.append([])
            all_cluster3_mean.append([])
            all_clusters3_channel.append([])
            all_clusters3_pe.append([])
            nobasecluster1_avganimal_pe_list.append([])
            nobasecluster2_avganimal_pe_list.append([])
            nobasecluster3_avganimal_pe_list.append([])
            
            for file0_ in range(len(file_list)):
        
                #ana_=1
                #file0_=2#file0='20110817B'
                #file0_=6#file0='20110913B'
                data1 = PE_all_list[ana_][file0_][:,1]
                num_cluster1, cluster_results1 = self.clustering_hierarchical(data1,numclusters)
                data2 = PE_all_list[ana_][file0_][:,2]
                num_cluster2, cluster_results2 = self.clustering_hierarchical(data2,numclusters)
                data3 = PE_all_list[ana_][file0_][:,3]
                num_cluster3, cluster_results3 = self.clustering_hierarchical(data3,numclusters)
                
                #print("1: num_cluster=%d results:"%num_cluster1,cluster_results1)
                #print("2: num_cluster=%d results:"%num_cluster2,cluster_results2)
                #print("3: num_cluster=%d results:"%num_cluster3,cluster_results3)
                
                colorlabel = ['black','red','green','yellow','blue','cyan','magenta']
                ### list clusters1_channel contains cluster members (chan) for each cluster number c
                ### list clusters1_pe contains values of PE for each cluster number c
                cluster1_mean = []
                clusters1_channel = []
                clusters1_pe = []
                for c in range(num_cluster1):
                    clusters1_channel.append([])
                    clusters1_pe.append([])
                for chan in range(self.num_chan):
                    cluster = cluster_results1[chan]-1
                    clusters1_channel[cluster].append(chan)
                    clusters1_pe[cluster].append(data1[chan])
                for c in range(num_cluster1):
                    cluster1_mean.append(np.mean(clusters1_pe[c]))
                
                cluster2_mean = []
                clusters2_channel = []
                clusters2_pe = []
                for c in range(num_cluster2):
                    clusters2_channel.append([])
                    clusters2_pe.append([])
                for chan in range(self.num_chan):
                    cluster = cluster_results2[chan]-1
                    clusters2_channel[cluster].append(chan)
                    clusters2_pe[cluster].append(data2[chan])
                for c in range(num_cluster2):
                    cluster2_mean.append(np.mean(clusters2_pe[c]))
                    
                cluster3_mean = []
                clusters3_channel = []
                clusters3_pe = []
                for c in range(num_cluster3):
                    clusters3_channel.append([])
                    clusters3_pe.append([])
                for chan in range(self.num_chan):
                    cluster = cluster_results3[chan]-1
                    clusters3_channel[cluster].append(chan)
                    clusters3_pe[cluster].append(data3[chan])
                for c in range(num_cluster3):
                    cluster3_mean.append(np.mean(clusters3_pe[c]))
                
                ################ sort cluster results according to the cluster means
                import copy
                srt_arglist1 = np.argsort(cluster1_mean,axis=0)
                #print("orig:",cluster1_mean,' arglist:',srt_arglist1)
                clusters1_channel_sorted = []
                clusters1_pe_sorted = []
                for c in range(num_cluster1):
                    clusters1_channel_sorted.append([])
                    clusters1_pe_sorted.append([])
                    #print("cond 1     cluster=%d  orig:"%c,clusters1_channel[c])
                for c in range(num_cluster1):
                    cluster=srt_arglist1[c]
                    #num = len(clusters3_channel[cluster])
                    clusters1_channel_sorted[c]=copy.deepcopy(clusters1_channel[cluster])
                    clusters1_pe_sorted[c]=copy.deepcopy(clusters1_pe[cluster])
                    #print("cluster=%d  sorted:"%cluster,clusters1_channel_sorted[c])
                cluster1_mean_sorted = np.sort(cluster1_mean,axis=0)
                
                srt_arglist2 = np.argsort(cluster2_mean,axis=0)
                #print("orig:",cluster2_mean,' arglist:',srt_arglist2)
                clusters2_channel_sorted = []
                clusters2_pe_sorted = []
                for c in range(num_cluster2):
                    clusters2_channel_sorted.append([])
                    clusters2_pe_sorted.append([])
                    #print("cond 2     cluster=%d  orig:"%c,clusters2_channel[c])
                for c in range(num_cluster2):
                    cluster=srt_arglist2[c]
                    #num = len(clusters3_channel[cluster])
                    clusters2_channel_sorted[c]=copy.deepcopy(clusters2_channel[cluster])
                    clusters2_pe_sorted[c]=copy.deepcopy(clusters2_pe[cluster])
                    #print("cluster=%d  sorted:"%cluster,clusters2_channel_sorted[cluster])
                cluster2_mean_sorted = np.sort(cluster2_mean,axis=0)
                
                srt_arglist3 = np.argsort(cluster3_mean,axis=0)
                #print("orig:",cluster3_mean,' arglist:',srt_arglist3)
                clusters3_channel_sorted = []
                clusters3_pe_sorted = []
                for c in range(num_cluster3):
                    clusters3_channel_sorted.append([])
                    clusters3_pe_sorted.append([])
                    #print("cond 3     cluster=%d  orig:"%c,clusters3_channel[c])
                for c in range(num_cluster3):
                    cluster=srt_arglist3[c]
                    #num = len(clusters3_channel[cluster])
                    clusters3_channel_sorted[c]=copy.deepcopy(clusters3_channel[cluster])
                    clusters3_pe_sorted[c]=copy.deepcopy(clusters3_pe[cluster])
                    #print("cluster=%d  sorted:"%cluster,clusters3_channel_sorted[cluster])
                cluster3_mean_sorted = np.sort(cluster3_mean,axis=0)
                
                all_cluster1_mean[ana_].append(cluster1_mean_sorted)
                all_clusters1_channel[ana_].append(clusters1_channel_sorted)
                all_clusters1_pe[ana_].append(clusters1_pe_sorted)
                all_cluster2_mean[ana_].append(cluster2_mean_sorted)
                all_clusters2_channel[ana_].append(clusters2_channel_sorted)
                all_clusters2_pe[ana_].append(clusters2_pe_sorted)
                all_cluster3_mean[ana_].append(cluster3_mean_sorted)
                all_clusters3_channel[ana_].append(clusters3_channel_sorted)
                all_clusters3_pe[ana_].append(clusters3_pe_sorted)
                
                
                nobasecluster1_pe_list = []
                nobasecluster2_pe_list = []
                nobasecluster3_pe_list = []
                cthr = 1
                for c in range(num_cluster1):
                    if c<cthr:
                        for pe in clusters1_pe_sorted[c]:
                            nobasecluster1_pe_list.append(pe)
                            nobasecluster1_avganimal_pe_list[ana_].append(pe)
                for c in range(num_cluster2):
                    if c<cthr:
                        for pe in clusters2_pe_sorted[c]:
                            nobasecluster2_pe_list.append(pe)
                            nobasecluster2_avganimal_pe_list[ana_].append(pe)
                for c in range(num_cluster3):
                    if c<cthr:
                        for pe in clusters3_pe_sorted[c]:
                            nobasecluster3_pe_list.append(pe)
                            nobasecluster3_avganimal_pe_list[ana_].append(pe)
                            #print("c=%d pe=%f"%(c,pe))
                
                #print("file0:",file0,"  nobasecluster3_pe_list:",nobasecluster3_pe_list)
                meds = []
                meds.append(stat.median(nobasecluster1_pe_list))
                meds.append(stat.median(nobasecluster2_pe_list))
                meds.append(stat.median(nobasecluster3_pe_list))
                
                pvals = []
                res = scipy.stats.mannwhitneyu(nobasecluster1_pe_list,nobasecluster2_pe_list)
                pvals.append(res[1])
                if res[1]<alpha:
                    print("ana:%d animal:%d  median1=%f  median2=%f  p=%f"%(ana_,file0_,meds[0],meds[1],res[1]))
                    
                res = scipy.stats.mannwhitneyu(nobasecluster2_pe_list,nobasecluster3_pe_list)
                pvals.append(res[1])
                if res[1]<alpha:
                    print("ana:%d animal:%d  median2=%f  median3=%f  p=%f"%(ana_,file0_,meds[1],meds[2],res[1]))
                    
                res = scipy.stats.mannwhitneyu(nobasecluster1_pe_list,nobasecluster3_pe_list)
                pvals.append(res[1])
                if res[1]<alpha:
                    print("ana:%d animal:%d  median1=%f  median3=%f  p=%f"%(ana_,file0_,meds[0],meds[2],res[1]))
                    
                ax = plt.subplot(3,7,counter)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.plot(np.arange(1,4,1),meds,'ko-')
                ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
                ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre','stim','post']))
                
                
#               for c in range(num_cluster1):
#                   plt.plot(clusters1_channel_sorted[c],clusters1_pe_sorted[c],ls='',marker='o',ms=4,color=colorlabel[c])
#                   print("c=%d"%c,clusters1_channel_sorted[c],clusters1_pe_sorted[c])
#               #print("cluster means of cond pre:",cluster1_mean_sorted)
#               for c in range(num_cluster2):
#                   plt.plot(clusters2_channel_sorted[c],clusters2_pe_sorted[c],ls='',marker='v',ms=4,color=colorlabel[c])
#               #print("cluster means of cond stim:",cluster2_mean_sorted)
#               for c in range(num_cluster3):
#                   plt.plot(clusters3_channel_sorted[c],clusters3_pe_sorted[c],ls='',marker='*',ms=4,color=colorlabel[c])
#               #print("cluster means of cond post:",cluster3_mean_sorted)
                counter += 1
                #print("ana_:%d file0:%d      nobasecluster1_avganimal_pe_list:"%(ana_,file0_),nobasecluster1_avganimal_pe_list)
        
        pe_mutual_channel = []
        pe_mutual_channel_all = []
        for file0 in range(7):
            pe_mutual_channel.append([])
            for ana_ in range(3):
                #print("cluster1:",all_clusters1_channel[ana_][file0])
                #print("cluster2:",all_clusters2_channel[ana_][file0])
                #print("cluster3:",all_clusters3_channel[ana_][file0])
                data1 = all_clusters1_channel[ana_][file0][0]
                data2 = all_clusters2_channel[ana_][file0][0]
                data3 = all_clusters3_channel[ana_][file0][0]
                list_ = []
                for k1 in range(len(data1)):
                    chan1 = data1[k1]
                    for k2 in range(len(data2)):
                        chan2 = data2[k2]
                        if chan1 == chan2:
                            for k3 in range(len(data3)):
                                chan3 = data3[k3]
                                if chan3 == chan2:
                                    list_.append(chan1)
                if len(list_)>0:
                    pe_mutual_channel[file0].append(list_)
                print("###### ana=%d file0=%d  mutual_channels:"%(ana_,file0),pe_mutual_channel[file0])
            list_ = []
            for k1 in range(len(pe_mutual_channel[file0][0])):
                chan1 = pe_mutual_channel[file0][0][k1]
                for k2 in range(len(pe_mutual_channel[file0][1])):
                    chan2 = pe_mutual_channel[file0][1][k2]
                    if chan1 == chan2:
                        for k3 in range(len(pe_mutual_channel[file0][2])):
                            chan3 = pe_mutual_channel[file0][2][k3]
                            if chan1 == chan3:
                                list_.append(chan1)
            if len(list_)>0:
                pe_mutual_channel_all.append(list_)
        quit()
        flag_plotnetwork = 1
        if flag_plotnetwork == 1:
            import networkx as nx
            # Create a complete graph with an odd number of nodes
            nnodes = self.num_chan
            connectionmatrix_list = []
            G_list = []
            pos_list = []
            for nsession in range(7):
                connectionmatrix=np.eye(nnodes)*0.0
                connectionmatrix_list.append(connectionmatrix)   
                G = nx.from_numpy_array(connectionmatrix)
                G_list.append(G)
                pos_list.append(nx.circular_layout(G))
            node_labels0_1 = {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7, 7: 8, \
                            8: 9, 9: 10, 10: 11, 11: 12, 12: 13, 13: 14, 14: 15, 15: 16,}
            node_labels0_2 = {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7, 7: 8, \
                            8: 9, 9: 10, 10: 11, 11: 12, 12: 13, 13: 14, 14: 15, 15: 16,}
            node_color = ["tab:red","tab:orange","tab:olive","tab:green","tab:blue","tab:purple","tab:cyan"]
            
            labelcolor_list = []
            labelcolor_list.append('white')
            labelcolor_list.append('black')
            labelcolor_list.append('white')
            
            node_list = []
            node_label_list = []
            
            for k in range(7):
                node_list.append([])
            #node_list0_1 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
            ana_ = 0
            for file0 in range(7):
                node_list[file0].append(all_clusters3_channel[ana_][file0][0])
                node_list[file0].append(all_clusters3_channel_rest[ana_][file0][0])
                #print("channel:",all_clusters3_channel[ana_][file0][0])
                #print("channel rest:",all_clusters3_channel_rest[ana_][file0][0])
                
            # Create a figure with 1:1 aspect ratio to preserve the circle.
            fig= plt.figure(figsize=(6, 10))
            for k in range(7):
                kk = k+1
                ax = plt.subplot(4,2,kk)
                node_opts_list = []
                node_opts_list.append({\
                    "node_size": 195, \
                    "node_color": node_color[k], \
                    "edgecolors": node_color[k], \
                    "linewidths": 1.0})
                node_opts_list.append({\
                    "node_size": 195, \
                    "node_color": "black", \
                    "edgecolors": node_color[k], \
                    "linewidths": 1.0})
                #nx.draw_networkx_nodes(G_list[k], pos_list[k], nodelist=[0, 1, 2, 3], **node_opts)
                #for l in range(len(node_list[k])):
                #print("l=%d"%l,"  ",labelcolor_list[l])
                l=0
                print("k=%d l=%d nodelist:"%(k,l),node_list[k][l])
                print("poslist:",pos_list[k])
                nx.draw_networkx_nodes(G_list[k], pos_list[k], nodelist=node_list[k][l],**node_opts_list[l])
                nx.draw_networkx_labels(G_list[k], pos_list[k], font_size=10, labels=node_labels0_1, font_color="black")#labelcolor_list[l])#labelcolor[l])
                nx.draw_networkx_edges(G_list[k], pos_list[k], width=1.0, edge_color=node_color[k], alpha=0.5)
                l=1
                print("k=%d l=%d rest nodelist:"%(k,l),node_list[k][l])
                print("rest poslist:",pos_list[k])
                nx.draw_networkx_nodes(G_list[k], pos_list[k], nodelist=node_list[k][l],**node_opts_list[l])
                nx.draw_networkx_labels(G_list[k], pos_list[k], font_size=10, labels=node_labels0_1, font_color="black")#labelcolor_list[l])#labelcolor[l])
                nx.draw_networkx_edges(G_list[k], pos_list[k], width=1.0, edge_color=node_color[k], alpha=0.5)
                #break
                ax.set_axis_off()
            fig.tight_layout()
            fig.subplots_adjust(hspace=0.06,wspace=0.1)
            plt.savefig("pe_nodes_ana0.png", dpi='figure', format='png')
            plt.show()
            quit()    
                    
                
        #quit()
        
        counter = 1
        plt.figure("PE grand average ")
        for ana_ in range(3):
            meds = []
            meds.append(stat.median(nobasecluster1_avganimal_pe_list[ana_]))
            meds.append(stat.median(nobasecluster2_avganimal_pe_list[ana_]))
            meds.append(stat.median(nobasecluster3_avganimal_pe_list[ana_]))
            
            pvals = []
            
            res = scipy.stats.mannwhitneyu(nobasecluster1_avganimal_pe_list[ana_],nobasecluster2_avganimal_pe_list[ana_])
            pvals.append(res[1])
            if res[1]<alpha:
                print("ana:%d avganimal  median1=%f  median2=%f  p=%f"%(ana_,meds[0],meds[1],res[1]))
            else:
                print("ana:%d avganimal  1-2  p=%f"%(ana_,res[1]))
            
            res = scipy.stats.mannwhitneyu(nobasecluster2_avganimal_pe_list[ana_],nobasecluster3_avganimal_pe_list[ana_])
            pvals.append(res[1])
            if res[1]<alpha:
                print("ana:%d avganimal  median2=%f  median3=%f  p=%f"%(ana_,meds[1],meds[2],res[1]))
            else:
                print("ana:%d avganimal  2-3  p=%f"%(ana_,res[1]))
                
            res = scipy.stats.mannwhitneyu(nobasecluster1_avganimal_pe_list[ana_],nobasecluster3_avganimal_pe_list[ana_])
            pvals.append(res[1])
            if res[1]<alpha:
                print("ana:%d avganimal  median1=%f  median3=%f  p=%f"%(ana_,meds[0],meds[2],res[1]))
            else:
                print("ana:%d avganimal  1-3  p=%f"%(ana_,res[1]))
                
            ax = plt.subplot(1,3,counter)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.plot(np.arange(1,4,1),meds,'ko-')
            ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre','stim','post']))
            counter += 1
        
        counter = 1
        plt.figure("PE grand average for experimental conditions")
        nobaseclusters = []
        nobaseclusters.append(nobasecluster1_avganimal_pe_list)
        nobaseclusters.append(nobasecluster2_avganimal_pe_list)
        nobaseclusters.append(nobasecluster3_avganimal_pe_list)
        meds = []
        meds.append([stat.median(nobasecluster1_avganimal_pe_list[0]),\
                    stat.median(nobasecluster1_avganimal_pe_list[1]),\
                    stat.median(nobasecluster1_avganimal_pe_list[2])])
        meds.append([stat.median(nobasecluster2_avganimal_pe_list[0]),\
                    stat.median(nobasecluster2_avganimal_pe_list[1]),\
                    stat.median(nobasecluster2_avganimal_pe_list[2])])
        meds.append([stat.median(nobasecluster3_avganimal_pe_list[0]),\
                    stat.median(nobasecluster3_avganimal_pe_list[1]),\
                    stat.median(nobasecluster3_avganimal_pe_list[2])])
        for cond in range(3):
            pvals = []
            
            res = scipy.stats.mannwhitneyu(nobaseclusters[cond][0],nobaseclusters[cond][1])
            pvals.append(res[1])
            if res[1]<alpha:
                print("expcond:%d avganimal  median1=%f  median2=%f  p=%f"%(cond,meds[cond][0],meds[cond][1],res[1]))
            else:
                print("expcond:%d avganimal  1-2  p=%f"%(cond,res[1]))
            
            res = scipy.stats.mannwhitneyu(nobaseclusters[cond][1],nobaseclusters[cond][2])
            pvals.append(res[1])
            if res[1]<alpha:
                print("expcond:%d avganimal  median1=%f  median2=%f  p=%f"%(cond,meds[cond][1],meds[cond][2],res[1]))
            else:
                print("expcond:%d avganimal  2-3  p=%f"%(cond,res[1]))
                
            res = scipy.stats.mannwhitneyu(nobaseclusters[cond][0],nobaseclusters[cond][2])
            pvals.append(res[1])
            if res[1]<alpha:
                print("expcond:%d avganimal  median1=%f  median2=%f  p=%f"%(cond,meds[cond][0],meds[cond][2],res[1]))
            else:
                print("expcond:%d avganimal  1-3  p=%f"%(cond,res[1]))
           
            ax = plt.subplot(1,3,counter)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.plot(np.arange(1,4,1),meds[cond],'ko-')
            ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(['low','medium','high']))
            counter += 1
        
        plt.show()
        
        
        
        #quit()
        
        
        
        fig00 = plt.figure("all data")
        counter = 1
        for ana_ in range(len(ana_list)):
            ana = ana_list[ana_]
            for file0_ in range(len(file_list)):
                ax = plt.subplot(3,len(file_list),counter)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.plot(range(16),PE_all_list[ana_][file0_][:,1],'ko-')
                plt.plot(range(16),PE_all_list[ana_][file0_][:,2],'ro-')
                plt.plot(range(16),PE_all_list[ana_][file0_][:,3],'bo-')
                counter += 1
        plt.show()
        quit()
                
        fig01 = plt.figure("single subject medians")
        for ana_ in range(len(ana_list)):
            ana = ana_list[ana_]
            for file0_ in range(len(file_list)):
                ax = plt.subplot(3,len(file_list),counter)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
#               plt.plot(np.arange(1,4,1),meds,'ko-',color=color[ana_])
#               ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
#               ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre','stim','post']))
#               
                plt.plot(np.arange(1,4,1),PE_chanmedian_list[ana_][file0_],'ko-')
                ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
                ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre-stim','stim-post','pre-post']))
                dscale=(maxy[file0_]-miny[file0_])*0.2
                ylimlow = miny[file0_]-dscale
                ylimhigh = maxy[file0_]+3*dscale
                plt.ylim(ylimlow,ylimhigh)
                #print("miny=%f maxy=%f"%(miny[file0_]-dscale,maxy[file0_]+dscale))
                ax.axhline(0, color='grey', linewidth=0.8)
                meds = PE_chanmedian_list[ana_][file0_]
                
                h = dscale*0.2
                y0 = ylimhigh-2*h
                #num_sign = len(signlist_list[file0_])
                num_sign = len(sign_diffslist_list[ana_][file0_])
                for kk in range(num_sign):
                    #sign = signlist_list[file0_][kk]
                    sign = sign_diffslist_list[ana_][file0_][kk]
                    print("sign_diffslist_list[%d][%d][%d]:"%(ana_,file0_,kk),sign_diffslist_list[ana_][file0_][kk])
                    x1 = sign#,sign[1]
                    y = PE_chanmedian_list[ana_][file0_][x1-1]+h
                    print("y:",y)
                    #print("kk:",kk," y:",y," y0:",y0," h:",h)
                    #print("kk=%d    y=%f  y0=%f   h=%f  maxy=%f   y+h=%f"%(kk,y,y0,h,max(meds),y+h))
                    #lab = '*'
                    #ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                    #ax.text((x1+x2)*.5, y+h,'*', ha='center', va='bottom', color='k')
                    ax.text(x1, y,'*', ha='center', va='bottom', color='k')
                counter += 1
        
        counter = 1
        fig01 = plt.figure("grand average medians")
        for ana_ in range(len(ana_list)):
            ana = ana_list[ana_]
            ax = plt.subplot(1,len(ana_list),counter)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            #for stype_ in range(len(stype_list)):
            plt.plot(np.arange(1,4,1),PE_chanpoolmedian_list[ana_],'ko-')
            ax.axhline(0, color='grey', linewidth=0.8)
            num_sign = len(sign_chanpoollist_list[ana_])
            #y = PE_chanpoolmedian_list[ana_]
            for kk in range(num_sign):
                #sign = signlist_list[file0_][kk]
                sign = sign_chanpoollist_list[ana_][kk]
                print("sign_chanpoollist_list[%d][%d]:"%(ana_,kk),sign_chanpoollist_list[ana_][kk])
                x1 = sign#,sign[1]
                print("PE_chanpoolmedian_list[ana_]:",PE_chanpoolmedian_list[ana_])
                y = PE_chanpoolmedian_list[ana_][x1-1]
                #print("kk:",kk," y:",y," y0:",y0," h:",h)
                #print("kk=%d    y=%f  y0=%f   h=%f  maxy=%f   y+h=%f"%(kk,y,y0,h,max(meds),y+h))
                #lab = '*'
                #ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                #ax.text((x1+x2)*.5, y+h,'*', ha='center', va='bottom', color='k')
                ax.text(x1, y,'*', ha='center', va='bottom', color='k')
            counter += 1
        
        plt.show()
        quit()
#                   #fig000 = plt.figure()
#                   #counter = 0
#                   pe_global_list = []
#                   pe_chanpool1_list = []
#                   pe_chanpool2_list = []
#                   pe_chanpool3_list = []
#                   for ana_ in range(3):
#                       pe_global_list.append([])
#                       pe_chanpool1_list.append([])
#                       pe_chanpool2_list.append([])
#                       pe_chanpool3_list.append([])
#                       for animal in range(7):
#                           pe = final_results[ana_][animal]
#                           pemean = []
#                           for l in range(3):
#                               pemean.append(np.mean(pe[l,:]))
#                           pe_global_list[ana_].append(pemean)
#                           for l in range(num_chan):
#                               pe_chanpool1_list[ana_].append(pe[0,l])
#                               pe_chanpool2_list[ana_].append(pe[1,l])
#                               pe_chanpool3_list[ana_].append(pe[2,l])
#                   
#                   color = ['ko-','ro-','bo-']
#                   counter = 1
#                   plt.figure("global pe")
#                   for animal in range(7):
#                       ax = plt.subplot(1,7,counter)
#                       for ana_ in range(3):
#                           list_ = pe_global_list[ana_][animal]
#                           plt.plot(range(3),list_,color[ana_])
#                       counter += 1
#                       
#                   import statistics as stat
#                   import scipy.stats
#                   counter = 1
#                   alpha = 0.05
#                   maxy=-100000000
#                   miny=100000000
#                   for ana_ in range(3):
#                       list1_ = pe_chanpool1_list[ana_]
#                       list2_ = pe_chanpool2_list[ana_]
#                       list3_ = pe_chanpool3_list[ana_]
#                       median1 = stat.median(list1_[:])
#                       median2 = stat.median(list2_[:])
#                       median3 = stat.median(list3_[:])
#                       meds = [median1,median2,median3]
#                       if max(meds)>maxy:
#                           maxy = max(meds)
#                       if min(meds)<miny:
#                           miny = min(meds)
#                   plt.figure("channel pool  pe")   
#                   for ana_ in range(3):
#                       ax = plt.subplot(1,3,counter)
#                       list1_ = pe_chanpool1_list[ana_]
#                       list2_ = pe_chanpool2_list[ana_]
#                       list3_ = pe_chanpool3_list[ana_]
#                       median1 = stat.median(list1_[:])
#                       median2 = stat.median(list2_[:])
#                       median3 = stat.median(list3_[:])
#                       meds = [median1,median2,median3]
#                       res12 = scipy.stats.mannwhitneyu(list1_, list2_)
#                       res23 = scipy.stats.mannwhitneyu(list2_, list3_)
#                       res13 = scipy.stats.mannwhitneyu(list1_, list3_)
#                       sign_list = []
#                       if res12[1]<alpha:
#                           print("ana:%d  1-2  medians: %f  %f   p=%f"%(analevel,median1,median2,res12[1]))
#                           sign_list.append([1,2])
#                       if res23[1]<alpha:
#                           print("ana:%d  2-3  medians: %f  %f   p=%f"%(analevel,median2,median3,res23[1]))
#                           sign_list.append([2,3])
#                       if res13[1]<alpha:
#                           print("ana:%d  1-3  medians: %f  %f   p=%f"%(analevel,median1,median3,res13[1]))
#                           sign_list.append([1,3])
#                       ax.spines['top'].set_visible(False)
#                       ax.spines['right'].set_visible(False)
#                       plt.plot(np.arange(1,4,1),meds,'ko-',color=color[ana_])
#                       ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
#                       ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre','stim','post']))
#                       h = 0.08*(max(meds)-min(meds))
#                       y0 = maxy*1.15#max(meds)+15*h
#                       num_sign = len(sign_list)
#                       for kk in range(num_sign):
#                           sign = sign_list[kk]
#                           x1, x2 = sign[0],sign[1]
#                           y = y0-kk*h*10
#                           #print("kk=%d    y=%f  y0=%f   h=%f  maxy=%f   y+h=%f"%(kk,y,y0,h,max(meds),y+h))
#                           #lab = '*'
#                           ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
#                           ax.text((x1+x2)*.5, y+h,'*', ha='center', va='bottom', color='k')
#                       plt.ylim(miny*0.97,maxy*1.2)   
#                       counter += 1
#                           
#                       plt.plot(range(3),meds,color[ana_])
#                       counter += 1
#                       
#                   plt.show()
#                   quit()
        
    
    
        
        
    
    
    def compute_IT_AIS(self,data_in):
        from idtxl.active_information_storage import ActiveInformationStorage
        #from idtxl.multivariate_te import MultivariateTE
        from idtxl.data import Data
        from idtxl.visualise_graph import plot_network
        import hde_fast_embedding
        
        data = Data(data_in,dim_order='spr')
        settings = {
            'cmi_estimator': 'JidtKraskovCMI',
            'n_perm_max_stat': 200,
            'n_perm_min_stat': 200,
            'max_lag': 10,
            'tau': 1
            }
            
        num_chan = np.shape(data_in)[1]
        #process = 2
        
        #processes = [0,4,5,6,7,8,9,10,11,12,13,15]# animal #1  (non-functional)
        #processes = [1,2,3,14]# animal #1  (functional)
        
        #processes = [3,6,9,10,11,12] #animal #2 (non-functional) 
        #processes = [0,1,2,4,5,7,8,13,14,15] #animal #2 (functional)
        
        #processes = [3,4,5,8,9,10,11,12]# animal #3 (non-functional)
        #processes = [0,1,2,6,7,13,14,15]# animal #3 - remaining (functional)
        
        #processes = [2,3,4,5,8,9,10,11,12,13]# animal #4 (non-functional)
        processes = [0,1,6,7,14,15]# animal #4 - remaining (functional)
        
        #processes = [1,2,3,4,5,6,9]# animal #5  (non-functional)
        #processes = [0,7,8,10,11,12,13,14,15]# animal #5  (functional)
        
        #processes = [2,3,4,5,6,9,10,11,12,13]# animal #6  (non-functional)
        #processes = [0,1,7,8,14,15]# animal #6  (functional)
        
        #processes = [1,2,3,4,13,14]# animal #7  (non-functional)
        #processes = [0,5,6,7,8,9,10,11,12,15]# animal #7 - remaining (functional)  
        
        network_analysis = ActiveInformationStorage()
        print("start AIS....")
        results = network_analysis.analyse_network(settings, data,processes)
        print(results)
        resp_list = []
        for process in processes:
            res_sp = results.get_single_process(process, fdr=False)
            resp_list.append(res_sp)
            print("process %d : "%process,res_sp)
        return resp_list
    
    def compute_IT_mTE(self,data_in):
        from idtxl.multivariate_te import MultivariateTE
        from idtxl.data import Data
        from idtxl.visualise_graph import plot_network
        
        data = Data(data_in,dim_order='spr')
        # Initialise analysis object and define settings
        network_analysis = MultivariateTE()
        settings = {'cmi_estimator': 'JidtGaussianCMI',
                'max_lag_sources': 5,
                'min_lag_sources': 1,
                'n_perm_max_stat': 200,
                'n_perm_min_stat': 200,}
        
        # Run analysis
        results = network_analysis.analyse_network(settings=settings, data=data)
        
        # Plot inferred network to console and via matplotlib
        results.print_edge_list(weights='max_te_lag', fdr=True)
        plot_network(results=results, weights='max_te_lag', fdr=True)
        plt.show()
        #quit()
    
    def analyse_InformationTheory(self,ana,file,stype,fs,duration):
        
        ########## compute multivariate Transfer Entropy
        #stype = 'pre'
        #ana = ana_list[0]
        path = './analysis/'+ana+'/'
        #file='20110913B'
        self.dt = 1.0/fs
        timewindow=duration
        if fs==200.0:
            label = path+file+'_'+ana+'_'+stype+'_notch.h5'
        if fs==500.0:
            label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
        print("datafile:",label)
        data=self.read_data_h5(label)
        if self.num_time<(timewindow/self.dt):
            print("self.num_time<(timewindow/self.dt): %d < %d    : choose duration %d"%(self.num_time,int(timewindow/self.dt),self.num_time))
            timewindow0=self.num_time*self.dt
        else:
            timewindow0 = timewindow
        num_times = int(timewindow0/self.dt)
        
        
        num_repetitions = 20
        repwindow = num_times#int(num_times//num_repetitions)
        #timewindow_final = repwindow*num_repetitions
        print("num_times:%d "%num_times)
        
        results_list = []
        label = "AIS_results_%s_%s_%s.dat"%(file,ana,stype)
        file_ = open(label,"w+")
        data_IT = np.zeros((repwindow,self.num_chan,1))
        for rep in range(num_repetitions):
            for chan in range(self.num_chan):
                for i in range(repwindow):
                    counter = np.shape(data)[0]-(num_repetitions+1)*repwindow+rep*repwindow+i
                    data_IT[i,chan,0]=data[counter,chan]
            #self.compute_IT_mTE(data_IT)
            res = self.compute_IT_AIS(data_IT)
            results_list.append(res)
            str = '%d  '%rep
            for res_ in res:
                print("%s %s %s   rep:%d AIS=%f p=%f"%(file,ana,stype,rep,res_.ais,res_.ais_pval))
                str += '   %f %f'%(res_.ais,res_.ais_pval)
            str += '\n'
            file_.write(str)
        file_.close()
    
    def compute_two_way_anova_allanimals_AIS(self,data0,file_list):
            ## balanced ANOVA 
            import pandas as pd
            import statistics as stat
            
            #data0 : [animals][cond][list_of_channels]
            print(np.shape(data0))
            num_trials = 20
            num_cond = 3
            num_animals  = 7
            num_animals_stimcond = num_animals*num_cond
            
            data_list = []
            col_animal = []
            col_cond = []
            col_chan = []
            str_cond = ['pre-low','post-low','post-high']
            str_animal = file_list
            for anim in range(num_animals):
                for n_cond in range(num_cond):
                    num_chan = len(data0[anim][n_cond])
                    print("n_cond:%d num_chan=%d"%(n_cond,num_chan))
                    str_chan = np.arange(1,num_chan+1,1)
                    for n_chan in range(num_chan):
                        print("n_cond:%d n_chan=%d"%(n_cond,n_chan))
                        for l in range(num_trials):
                            print("n_cond:%d n_chan=%d trial:%d"%(n_cond,n_chan,l))
                            data_list.append(data0[anim][n_cond][n_chan][l])
                            col_animal.append(str_animal[anim])
                            col_cond.append(str_cond[n_cond])
                            col_chan.append(str_chan[n_chan])
            df = pd.DataFrame({
                            'animal': col_animal,
                            'cond': col_cond,\
                            'channel': col_chan,\
                            'data': data_list})
            print(df[::])
            #quit()
            import statsmodels.api as sm
            from statsmodels.formula.api import ols
            #perform two-way ANOVA with interaction level*stim
            model = ols('data ~ C(animal) + C(cond) + C(animal):C(cond)', data=df).fit()
            res = sm.stats.anova_lm(model, typ=3) # typ = 1,2,3 are all equivalent for balanced ANOVA
            print("ANOVA with two factors C(animal) + C(cond) and interaction animal*cond")
            print(res)
            #perform two-way ANOVA with interaction level*animal
            #model = ols('medianPLV ~ C(a_level) + C(stim_cond) + C(animal) + C(a_level):C(animal)', data=df).fit()
            #res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
            #print("ANOVA with three factors C(a_level) + C(stim_cond) + C(animal) and interaction level*animal")
            #print(res)
            
            from statsmodels.stats.multicomp import pairwise_tukeyhsd
            ### if results depend on category 'animal', then one should specify the animal.
            print("Post hoc statistics: the stimulations")
            tukey = pairwise_tukeyhsd(endog=df['data'],groups=df['cond'], alpha=0.05)
            print(tukey)
            quit()
    
    def compute_two_way_anova_allanimals(self,data0):
        ## balanced ANOVA 
        import pandas as pd
        import statistics as stat
        num_analevel = 3
        num_stimcond = 3
        num_animals  = 7
        num_animals_stimcond = num_animals*num_stimcond
        num_animals_stimcond = num_animals*num_stimcond
        
        data_list = []
        
        col_ana = []
        col_stim = []
        str_ana = ['low','medium','high']
        str_stim = ['pre', 'stim', 'post']
        for n_ana in range(num_analevel):
            for n_stim in range(num_stimcond):
                for n_animal in range(num_animals):
                    data_list.append(data0[n_ana][n_stim][n_animal])
                    col_ana.append(str_ana[n_ana])
                    col_stim.append(str_stim[n_stim])
        df = pd.DataFrame({'a_level': col_ana,\
                        'stim_cond': col_stim,\
                        'data': data_list})
        print(df[::])
        #quit()
        import statsmodels.api as sm
        from statsmodels.formula.api import ols
        #perform two-way ANOVA with interaction level*stim
        model = ols('data ~ C(a_level) + C(stim_cond) + C(a_level):C(stim_cond)', data=df).fit()
        res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
        print("ANOVA with two factors C(a_level) + C(stim_cond) and interaction level*stim")
        print(res)
        #perform two-way ANOVA with interaction level*animal
        #model = ols('medianPLV ~ C(a_level) + C(stim_cond) + C(animal) + C(a_level):C(animal)', data=df).fit()
        #res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
        #print("ANOVA with three factors C(a_level) + C(stim_cond) + C(animal) and interaction level*animal")
        #print(res)
        
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
        ### if results depend on category 'animal', then one should specify the animal.
        print("Post hoc statistics: the stimulations")
        tukey = pairwise_tukeyhsd(endog=df['data'],groups=df['stim_cond'], alpha=0.05)
        print(tukey)
        print("Post hoc statistics: the anaesthesia levels")
        tukey = pairwise_tukeyhsd(endog=df['data'],groups=df['a_level'], alpha=0.05)
        print(tukey)


    def compute_one_way_anova_allanimals_fmin(self,data0):
        ## balanced ANOVA 
        import pandas as pd
        import statistics as stat
        num_fmin = 2
        num_analevel = 3
        num_stimcond = 3
        num_animals  = 7
        
        data_list = []
        
        col_fmin = []
        col_ana  = []
        col_stim = []
        str_fmin = ['slow','gamma']
        str_ana  = ['low','medium','high']
        str_stim = ['pre', 'stim', 'post']
        for n_fmin in range(num_fmin):
            for n_ana in range(num_analevel):
                for n_stim in range(num_stimcond):
                    for n_animal in range(num_animals):
                        data_list.append(data0[n_fmin][n_ana][n_stim][n_animal])
                        col_fmin.append(str_fmin[n_fmin])
                        col_ana.append(str_ana[n_ana])
                        col_stim.append(str_stim[n_stim])
                        #print("col_fmin:",len(col_fmin))
        print(len(col_fmin),len(col_ana),len(col_stim),len(data_list))
        df = pd.DataFrame({'band': col_fmin,\
                        'a_level': col_ana,\
                        'stim_cond': col_stim,\
                        'data': data_list})
        print(df[::])
        #quit()
        import statsmodels.api as sm
        from statsmodels.formula.api import ols
        #perform two-way ANOVA with interaction level*stim
        model = ols('data ~ C(band) ', data=df).fit()
        res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
        print("ANOVA with one factor 'band' ")
        print(res)
        #perform two-way ANOVA with interaction level*animal
        #model = ols('medianPLV ~ C(a_level) + C(stim_cond) + C(animal) + C(a_level):C(animal)', data=df).fit()
        #res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
        #print("ANOVA with three factors C(a_level) + C(stim_cond) + C(animal) and interaction level*animal")
        #print(res)
        
        from statsmodels.stats.multicomp import pairwise_tukeyhsd
        ### if results depend on category 'animal', then one should specify the animal.
        print("Post hoc statistics: the stimulations")
        tukey = pairwise_tukeyhsd(endog=df['data'],groups=df['band'], alpha=0.05)
        print(tukey)
    
    
    def compute_two_way_anova(self,data0):
        import pandas as pd
        import statistics as stat
        #create data
        #df = pd.DataFrame({'water': np.repeat(['daily', 'weekly'], 15),
        #            'sun': np.tile(np.repeat(['low', 'med', 'high'], 5), 2),
        #            'height': [6, 6, 6, 5, 6, 5, 5, 6, 4, 5,
        #                    6, 6, 7, 8, 7, 3, 4, 4, 4, 5,
        #                    4, 4, 4, 4, 4, 5, 6, 6, 7, 8]})
        #print(df)
        num_analevel = 3
        num_stimcond = 3
        num_animals  = 7
        num_animals_stimcond = num_animals*num_stimcond
        num_animals_stimcond = num_animals*num_stimcond
        ##note: we have three categorial parameters: analevel, stimcond and animal !!!            
        ## But this will be difficult since there are different numbers of elements in each 
        ## category, i.e. different number of cluster members which yields an unbalanced ANOVA.
        ## Consequently, for a balanced ANOVA (and thus stronger statistical test) one should 
        ## start with two categories analevel and stimcond and 
        ## consider the median of cluster members to gain one observable in each animal group.
        ## This is a crude approximation, but it fits to the previous animal pooling.
        data_list = []
        
        flag_balanedANOVA = 0
        
        if flag_balanedANOVA == 1:
            ## balanced ANOVA
            for n_ana in range(num_analevel):
                for n_stim in range(num_stimcond):
                    for n_animal in range(num_animals):
                        data_list.append(stat.median(data0[n_ana][n_animal][n_stim]))
            df = pd.DataFrame({'a_level': np.repeat(['low', 'medium','high'], num_animals_stimcond),\
                            'stim_cond': np.tile(np.repeat(['pre', 'stim', 'post'], num_animals),num_stimcond),\
                                'data': data_list})
            #print(df[:23])
            #quit()
            import statsmodels.api as sm
            from statsmodels.formula.api import ols
            #perform two-way ANOVA
            model = ols('data ~ C(a_level) + C(stim_cond) + C(a_level):C(stim_cond)', data=df).fit()
            res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
            print(model)
            print(res)
        
        if flag_balanedANOVA == 0:
            ## unbalanced ANOVA
            #d = {'col1': [1, 2], 'col2': [3, 4]}
            
            col_ana = []
            col_stim = []
            col_animal = []
            str_ana = ['low','medium','high']
            str_stim = ['pre', 'stim', 'post']
            for n_ana in range(num_analevel):
                for n_stim in range(num_stimcond):
                    for n_animal in range(num_animals):
                        #print("%d %d %d   %d"%(n_ana,n_stim,n_animal,len(data0[n_ana][n_animal][n_stim])))
                        for n_plv in range(len(data0[n_ana][n_animal][n_stim])):
                            data_list.append(data0[n_ana][n_animal][n_stim][n_plv])
                            col_ana.append(str_ana[n_ana])
                            col_stim.append(str_stim[n_stim])
                            col_animal.append(n_animal)
            df = pd.DataFrame({'a_level': col_ana,\
                             'stim_cond': col_stim,\
                             'animal': col_animal,\
                            'data': data_list})
            #print(df[:20])
            #quit()
            import statsmodels.api as sm
            from statsmodels.formula.api import ols
            #perform two-way ANOVA with interaction level*stim
            model = ols('data ~ C(a_level) + C(stim_cond) + C(animal) + C(a_level):C(stim_cond)', data=df).fit()
            res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
            print("ANOVA with three factors C(a_level) + C(stim_cond) + C(animal) and interaction level*stim")
            print(res)
            #perform two-way ANOVA with interaction level*animal
            #model = ols('medianPLV ~ C(a_level) + C(stim_cond) + C(animal) + C(a_level):C(animal)', data=df).fit()
            #res = sm.stats.anova_lm(model, typ=2) # typ = 1,2,3 are all equivalent for balanced ANOVA
            #print("ANOVA with three factors C(a_level) + C(stim_cond) + C(animal) and interaction level*animal")
            #print(res)
            
            from statsmodels.stats.multicomp import pairwise_tukeyhsd
            ### if results depend on category 'animal', then one should specify the animal.
            for n_animal in range(num_animals):
                print ("Now let us focus on animal #%d"%n_animal)
                
                col_ana = []
                col_stim = []
                data_list = []
                for n_ana in range(num_analevel):
                    for n_stim in range(num_stimcond):
                        for n_plv in range(len(data0[n_ana][n_animal][n_stim])):
                            data_list.append(data0[n_ana][n_animal][n_stim][n_plv])
                            col_ana.append(str_ana[n_ana])
                            col_stim.append(str_stim[n_stim])
                df = pd.DataFrame({'a_level': col_ana,\
                                'stim_cond': col_stim,\
                                'medianPLV': data_list})
                #print(df[:160])
                model = ols('medianPLV ~ C(a_level) + C(stim_cond) + C(a_level):C(stim_cond)', data=df).fit()
                res = sm.stats.anova_lm(model, typ=3) # typ = 1,2,3 are all equivalent for balanced ANOVA
                print("ANOVA with two factors C(a_level) + C(stim_cond)")
                print(res)
                
                print("Post hoc statistics: the stimulations")
                tukey = pairwise_tukeyhsd(endog=df['medianPLV'],groups=df['stim_cond'], alpha=0.05)
                print(tukey)
                print("Post hoc statistics: the anaesthesia levels")
                tukey = pairwise_tukeyhsd(endog=df['medianPLV'],groups=df['a_level'], alpha=0.05)
                print(tukey)
            ############## results for gamma ##############
            #                            sum_sq      df           F         PR(>F)
            #C(a_level)                0.080169     2.0    2.473195   8.442365e-02
            #C(stim_cond)              0.002632     2.0    0.081212   9.219997e-01
            #C(animal)                19.382448     6.0  199.315698  2.824112e-227
            #C(a_level):C(stim_cond)   0.216194     4.0    3.334790   9.799667e-03
            #Residual                 76.953180  4748.0         NaN            NaN
#           Now let us focus on animal #0
#           ANOVA with two factors C(a_level) + C(stim_cond)
#                                       sum_sq     df            F         PR(>F)
#           Intercept                18.692773    1.0  2019.513017  1.361204e-185
#           C(a_level)                0.323209    2.0    17.459295   4.468084e-08
#           C(stim_cond)              0.129507    2.0     6.995758   1.000130e-03
#           C(a_level):C(stim_cond)   0.340855    4.0     9.206254   3.333569e-07
#           Residual                  5.053819  546.0          NaN            NaN
#           Post hoc statistics: the stimulations
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               post    pre  -0.0081  0.679 -0.0315 0.0153  False
#               post   stim   0.0294 0.0177  0.0041 0.0547   True
#               pre   stim   0.0375 0.0013  0.0125 0.0625   True
#           ---------------------------------------------------
#           Post hoc statistics: the anaesthesia levels
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               high    low  -0.0052 0.8586 -0.0295 0.0192  False
#               high medium   0.0262 0.0399   0.001 0.0514   True
#               low medium   0.0314 0.0063  0.0074 0.0554   True
#           ---------------------------------------------------
#           Now let us focus on animal #1
#           ANOVA with two factors C(a_level) + C(stim_cond)
#                                       sum_sq     df            F         PR(>F)
#           Intercept                34.976754    1.0  2054.953223  4.821526e-214
#           C(a_level)                0.118353    2.0     3.476741   3.142286e-02
#           C(stim_cond)              0.087500    2.0     2.570402   7.719991e-02
#           C(a_level):C(stim_cond)   0.130557    4.0     1.917620   1.056595e-01
#           Residual                 12.374053  727.0          NaN            NaN
#           Post hoc statistics: the stimulations
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05 
#           ====================================================
#           group1 group2 meandiff p-adj   lower   upper  reject
#           ----------------------------------------------------
#               post    pre   0.0246 0.1068 -0.0039  0.0531  False
#               post   stim  -0.0311 0.0207 -0.0584 -0.0038   True
#               pre   stim  -0.0557  0.001 -0.0833 -0.0281   True
#           ----------------------------------------------------
#           Post hoc statistics: the anaesthesia levels
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               high    low  -0.0045    0.9 -0.0327 0.0237  False
#               high medium   0.0082 0.7562 -0.0202 0.0365  False
#               low medium   0.0127 0.5334 -0.0152 0.0405  False
#           ---------------------------------------------------
#           Now let us focus on animal #2
#           ANOVA with two factors C(a_level) + C(stim_cond)
#                                       sum_sq     df           F         PR(>F)
#           Intercept                19.081440    1.0  918.149647  2.998047e-121
#           C(a_level)                0.245924    2.0    5.916618   2.861068e-03
#           C(stim_cond)              0.039709    2.0    0.955338   3.852919e-01
#           C(a_level):C(stim_cond)   0.612799    4.0    7.371571   8.564179e-06
#           Residual                 11.949935  575.0         NaN            NaN
#           Post hoc statistics: the stimulations
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               post    pre  -0.0064 0.8994 -0.0418 0.0289  False
#               post   stim   0.0224 0.3155 -0.0139 0.0588  False
#               pre   stim   0.0289 0.1294 -0.0062 0.0639  False
#           ---------------------------------------------------
#           Post hoc statistics: the anaesthesia levels
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               high    low   0.0565  0.001  0.0211  0.092   True
#               high medium   0.0092 0.7872 -0.0255 0.0439  False
#               low medium  -0.0473  0.005 -0.0827 -0.012   True
#           ---------------------------------------------------
#           Now let us focus on animal #3
#           ANOVA with two factors C(a_level) + C(stim_cond)
#                                       sum_sq     df            F         PR(>F)
#           Intercept                24.435245    1.0  1034.309881  3.714475e-142
#           C(a_level)                0.230977    2.0     4.888467   7.780707e-03
#           C(stim_cond)              0.167584    2.0     3.546796   2.931311e-02
#           C(a_level):C(stim_cond)   0.117591    4.0     1.244368   2.906595e-01
#           Residual                 17.293269  732.0          NaN            NaN
#           Post hoc statistics: the stimulations
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               post    pre   0.0208 0.2855 -0.0115  0.053  False
#               post   stim   0.0197 0.3382 -0.0132 0.0525  False
#               pre   stim  -0.0011    0.9 -0.0339 0.0317  False
#           ---------------------------------------------------
#           Post hoc statistics: the anaesthesia levels
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               high    low   0.0349 0.0322  0.0024 0.0675   True
#               high medium   0.0241 0.1914 -0.0084 0.0565  False
#               low medium  -0.0109  0.694 -0.0434 0.0216  False
#           ---------------------------------------------------
#           Now let us focus on animal #4
#           ANOVA with two factors C(a_level) + C(stim_cond)
#                                       sum_sq     df            F         PR(>F)
#           Intercept                42.030088    1.0  3377.980305  2.789792e-270
#           C(a_level)                0.004155    2.0     0.166970   8.462587e-01
#           C(stim_cond)              0.005118    2.0     0.205676   8.141462e-01
#           C(a_level):C(stim_cond)   0.028558    4.0     0.573805   6.817273e-01
#           Residual                  8.722103  701.0          NaN            NaN
#           Post hoc statistics: the stimulations
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ==================================================
#           group1 group2 meandiff p-adj  lower  upper  reject
#           --------------------------------------------------
#               post    pre   0.0038   0.9 -0.0203 0.0279  False
#               post   stim   0.0019   0.9  -0.022 0.0258  False
#               pre   stim  -0.0019   0.9  -0.026 0.0223  False
#           --------------------------------------------------
#           Post hoc statistics: the anaesthesia levels
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               high    low   0.0126 0.4379 -0.0115 0.0368  False
#               high medium   0.0005    0.9 -0.0234 0.0244  False
#               low medium  -0.0122 0.4632 -0.0362 0.0119  False
#           ---------------------------------------------------
#           Now let us focus on animal #5
#           ANOVA with two factors C(a_level) + C(stim_cond)
#                                       sum_sq     df            F         PR(>F)
#           Intercept                31.048551    1.0  1778.816944  1.332586e-188
#           C(a_level)                0.109896    2.0     3.148053   4.358799e-02
#           C(stim_cond)              0.134092    2.0     3.841151   2.195562e-02
#           C(a_level):C(stim_cond)   0.166420    4.0     2.383605   5.016080e-02
#           Residual                 11.397858  653.0          NaN            NaN
#           Post hoc statistics: the stimulations
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               post    pre  -0.0097 0.7026 -0.0393 0.0199  False
#               post   stim  -0.0092 0.7269 -0.0389 0.0205  False
#               pre   stim   0.0005    0.9 -0.0291 0.0301  False
#           ---------------------------------------------------
#           Post hoc statistics: the anaesthesia levels
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ==================================================
#           group1 group2 meandiff p-adj  lower  upper  reject
#           --------------------------------------------------
#               high    low   0.0035   0.9 -0.0262 0.0332  False
#               high medium   0.0002   0.9 -0.0292 0.0296  False
#               low medium  -0.0033   0.9 -0.0332 0.0266  False
#           --------------------------------------------------
#           Now let us focus on animal #6
#           ANOVA with two factors C(a_level) + C(stim_cond)
#                                       sum_sq     df            F         PR(>F)
#           Intercept                35.670948    1.0  3806.615437  2.069808e-299
#           C(a_level)                0.198941    2.0    10.614975   2.835976e-05
#           C(stim_cond)              0.256344    2.0    13.677864   1.456765e-06
#           C(a_level):C(stim_cond)   0.216469    4.0     5.775106   1.395554e-04
#           Residual                  7.178016  766.0          NaN            NaN
#           Post hoc statistics: the stimulations
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05 
#           ====================================================
#           group1 group2 meandiff p-adj   lower   upper  reject
#           ----------------------------------------------------
#               post    pre  -0.0349  0.001 -0.0554 -0.0145   True
#               post   stim  -0.0247 0.0149 -0.0455 -0.0039   True
#               pre   stim   0.0102 0.4612 -0.0099  0.0304  False
#           ----------------------------------------------------
#           Post hoc statistics: the anaesthesia levels
#           Multiple Comparison of Means - Tukey HSD, FWER=0.05
#           ===================================================
#           group1 group2 meandiff p-adj   lower  upper  reject
#           ---------------------------------------------------
#               high    low  -0.0195 0.0692 -0.0402 0.0012  False
#               high medium  -0.0354  0.001 -0.0557 -0.015   True
#               low medium  -0.0159 0.1599 -0.0362 0.0045  False
#           ---------------------------------------------------
#           
                    
#           for n_stimcond in range(num_stimcond):
#               print ("Now let us focus on stimulation condition %d"%n_stimcond)
#               
#               col_ana = []
#               col_animal = []
#               data_list = []
#               for n_ana in range(num_analevel):
#                   for n_animal in range(num_animals):
#                       for n_plv in range(len(data0[n_ana][n_animal][n_stim])):
#                           data_list.append(data0[n_ana][n_animal][n_stim][n_plv])
#                           col_ana.append(str_ana[n_ana])
#                           col_animal.append(n_animal)
#               df = pd.DataFrame({'a_level': col_ana,\
#                               'animal': col_animal,\
#                               'medianPLV': data_list})
#               #print(df[:160])
#               model = ols('medianPLV ~ C(a_level) + C(animal) + C(a_level):C(animal)', data=df).fit()
#               res = sm.stats.anova_lm(model, typ=3) # typ = 1,2,3 are all equivalent for balanced ANOVA
#               print("ANOVA with two factors C(a_level) + C(animal)")
#               print(res)
#               
#               print("Post hoc statistics: the animals")
#               tukey = pairwise_tukeyhsd(endog=df['medianPLV'],groups=df['animal'], alpha=0.05)
#               print(tukey)
    
            
        #quit()

    def compute_connectivity_statistics_all(self,ana_list,fmin,fmax,fs,flags):
        from scipy.stats import circstd
        from matplotlib import ticker
        
        flag_plots,flag_inference,flag_twofmin = flags[0],flags[1],flags[2]
        
        ############## summarize projections on network interactions
        counter = 1
        #plt.figure("projections for all animals and levels and stimulations")
        data_foranova_nn = []
        data_foranova_fn = []
        for ana in range(3):
            data_foranova_nn.append([])
            data_foranova_fn.append([])
            for animal in range(7):
                res_nn = self.projections_list[ana][animal][0]#nearest neighbor
                res_fn = self.projections_list[ana][animal][0]#rest
                #res_foranova = self.projections_list[ana][animal][1][0]## first projection mode
                #print(np.shape(self.projections_list[ana][animal][1]))
                #print(res_foranova)
                #quit()
                data_foranova_nn[ana].append(res_nn)
                data_foranova_fn[ana].append(res_fn)
                #ax = plt.subplot(3,7,counter)
                #plt.plot(range(3),res[0],'ko-',label='proj1')
                #plt.plot(range(3),res[1],'ro-',label='proj2')
                #plt.plot(range(3),res[2],'bo-',label='proj3')
                #counter += 1
        #plt.show()
        #print(np.shape(data_foranova))
        #print(np.shape(data_foranova[2][6][2]))
        #print(data_foranova[2][6][2])
        #quit()
        print("####################### nearest neighbor ANOVA ######################")
        #self.compute_two_way_anova(data_foranova_nn)
        print("####################### far neighbors ANOVA ######################")
        #self.compute_two_way_anova(data_foranova_fn)
        ############## end of : summarize projections on network interactions
        
        ############## summarize symmetry check ########################
        if flag_twofmin == 1:
            fmin_index = -1
            symmetry_all_list = []
            cumul = []
            for fminindex in range(2):
                #print("sizes:",np.size(self.symmetry_list[fminindex]))
                symmetry_all_list.append([])
                l_list=[]
                for ana in range(3):
                    symmetry_all_list[fminindex].append([])
                    for stim in range(3):
                        l_ = []
                        for animals in range(7):
                            #print("fminindex:%d  ana:%d stim=%d animal:%d"%(fminindex,ana,stim,animals))
                            dat = self.symmetry_list[fminindex][ana][animals]
                            #print("ana:%d stim=%d animal:%d"%(ana,stim,animals),dat)
                            l_.append(dat[stim])
                        symmetry_all_list[fminindex][ana].append(l_)
                        l_list = l_list + l_
                cumul.append(l_list)
            import scipy.stats
            res = scipy.stats.mannwhitneyu(cumul[0],cumul[1])
            print("p-value: ",res[1])
            plt.figure("all ratios")
            ax = plt.subplot(1,1,1)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.boxplot(cumul)
            plt.show()
            quit()
            
            
        if flag_twofmin == 0:
            symmetry_all_list = []
            cumul = []
            for ana in range(3):
                symmetry_all_list.append([])
                for stim in range(3):
                    l_ = []
                    for animals in range(7):
                        #print("fminindex:%d  ana:%d stim=%d animal:%d"%(fminindex,ana,stim,animals))
                        dat = self.symmetry_list[0][ana][animals]
                        #print("ana:%d stim=%d animal:%d"%(ana,stim,animals),dat)
                        l_.append(dat[stim])
                    symmetry_all_list[ana].append(l_)
            self.compute_two_way_anova_allanimals(symmetry_all_list)
            
        #print(np.shape(symmetry_all_list))
        #self.compute_two_way_anova_allanimals(symmetry_all_list)
        
        #self.compute_one_way_anova_allanimals_fmin(symmetry_all_list)
        #quit()
        ############## end of : summarize symmetry check ########################
        
        
        
        num1 = len(self.topcluster_meanphase1_sessionpool_list)
        num2 = len(self.topcluster_meanphase2_sessionpool_list)
        num3 = len(self.topcluster_meanphase3_sessionpool_list)
        phasemean1 = np.mean(self.topcluster_meanphase1_sessionpool_list)
        phasemean2 = np.mean(self.topcluster_meanphase2_sessionpool_list)
        phasemean3 = np.mean(self.topcluster_meanphase3_sessionpool_list)
        phasestd1 = circstd(self.topcluster_meanphase1_sessionpool_list)
        phasestd2 = circstd(self.topcluster_meanphase2_sessionpool_list)
        phasestd3 = circstd(self.topcluster_meanphase3_sessionpool_list)
        #print("cond 1 : mean:",phasemean1," stddev:",phasestd1," num:",num1)
        #print("cond 2 : mean:",phasemean2," stddev:",phasestd2," num:",num2)
        #print("cond 3 : mean:",phasemean3," stddev:",phasestd3," num:",num3)
        

        counter = 0
        h1 = np.zeros(7)
        yerr1 = np.zeros(7)
        h2 = np.zeros(7)
        yerr2 = np.zeros(7)
        h3 = np.zeros(7)
        yerr3 = np.zeros(7)
        for s in range(7):
            for ana_ in range(3):
                h1[s]+=self.topcluster_meanphase1_singlesession_list[ana_][s]/3.0
                yerr1[s]+=self.topcluster_stdphase1_singlesession_list[ana_][s]/3.0
                h2[s]+=self.topcluster_meanphase2_singlesession_list[ana_][s]/3.0
                yerr2[s]+=self.topcluster_stdphase2_singlesession_list[ana_][s]/3.0
                h3[s]+=self.topcluster_meanphase3_singlesession_list[ana_][s]/3.0
                yerr3[s]+=self.topcluster_stdphase3_singlesession_list[ana_][s]/3.0
                #print("ana=%d yerr=%f"%(ana_,self.topcluster_stdphase1_singlesession_list[ana_][s]) )
            #print(yerr1,yerr2,yerr3)        

        labels = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7']
        x = np.arange(len(labels))  # the label locations
        width = 0.25  # the width of the bars
        
        #fig, axs = plt.subplots(1, 3, figsize=(9, 3), sharey=True)
        #axs[0].bar(names, values)
        #axs[1].scatter(names, values)
        #axs[2].plot(names, values)
        if flag_plots == 2:
            fig, ax = plt.subplots(3,1)
            #ax[0].axis('off')
            #ax[1].axis('off')
            #ax[2].axis('off')
            ax[0].spines['top'].set_visible(False)
            ax[0].spines['right'].set_visible(False)
            ax[1].spines['top'].set_visible(False)
            ax[1].spines['right'].set_visible(False)
            ax[2].spines['top'].set_visible(False)
            ax[2].spines['right'].set_visible(False)
            
            h0=self.topcluster_meanphase1_singlesession_list[0]
            yerr0=self.topcluster_stdphase1_singlesession_list[0]
            h1=self.topcluster_meanphase2_singlesession_list[0]
            yerr1=self.topcluster_stdphase2_singlesession_list[0]
            h2=self.topcluster_meanphase3_singlesession_list[0]
            yerr2=self.topcluster_stdphase3_singlesession_list[0]
            max1 = np.max(self.topcluster_stdphase1_singlesession_list[0])
            max2 = np.max(self.topcluster_stdphase2_singlesession_list[0])
            max3 = np.max(self.topcluster_stdphase3_singlesession_list[0])
            maxmax=max(max1,max2,max3)
            rects01 = ax[0].bar(x - width, h0, width, align='center', yerr=yerr0, ecolor='lightgrey',label='pre')
            rects02 = ax[0].bar(x        , h1, width, align='center', yerr=yerr1, ecolor='lightgrey', label='stim')
            rects03 = ax[0].bar(x + width, h2, width, align='center', yerr=yerr2, ecolor='lightgrey', label='post')
            ax[0].axhline(0, color='grey', linewidth=0.8)
            ax[0].xaxis.set_major_locator(ticker.FixedLocator(x))
            ax[0].xaxis.set_major_formatter(ticker.FixedFormatter([]))
            ax[0].yaxis.set_major_locator(ticker.FixedLocator([-maxmax,0,maxmax]))
            #ax[0].yaxis.set_major_formatter(ticker.FixedFormatter([]))
            #ax[0].legend(frameon='off',mode='expand')
            
            h0=self.topcluster_meanphase1_singlesession_list[1]
            yerr0=self.topcluster_stdphase1_singlesession_list[1]
            h1=self.topcluster_meanphase2_singlesession_list[1]
            yerr1=self.topcluster_stdphase2_singlesession_list[1]
            h2=self.topcluster_meanphase3_singlesession_list[1]
            yerr2=self.topcluster_stdphase3_singlesession_list[1]
            max1 = np.max(self.topcluster_stdphase1_singlesession_list[1])
            max2 = np.max(self.topcluster_stdphase2_singlesession_list[1])
            max3 = np.max(self.topcluster_stdphase3_singlesession_list[1])
            maxmax=max(max1,max2,max3)
            rects11 = ax[1].bar(x - width, h0, width, align='center', yerr=yerr0, ecolor='lightgrey')
            rects12 = ax[1].bar(x        , h1, width, align='center', yerr=yerr1, ecolor='lightgrey')
            rects13 = ax[1].bar(x + width, h2, width, align='center', yerr=yerr2, ecolor='lightgrey')
            ax[1].axhline(0, color='grey', linewidth=0.8)
            ax[1].xaxis.set_major_locator(ticker.FixedLocator(x))
            ax[1].xaxis.set_major_formatter(ticker.FixedFormatter([]))
            ax[1].yaxis.set_major_locator(ticker.FixedLocator([-maxmax,0,maxmax]))
            #ax[1].yaxis.set_major_formatter(ticker.FixedFormatter([]))
            #ax[1].legend()
            
            h0=self.topcluster_meanphase1_singlesession_list[2]
            yerr0=self.topcluster_stdphase1_singlesession_list[2]
            h1=self.topcluster_meanphase2_singlesession_list[2]
            yerr1=self.topcluster_stdphase2_singlesession_list[2]
            h2=self.topcluster_meanphase3_singlesession_list[2]
            yerr2=self.topcluster_stdphase3_singlesession_list[2]
            #min1 = np.min(self.topcluster_stdphase1_singlesession_list[2])
            #min2 = np.min(self.topcluster_stdphase2_singlesession_list[2])
            #min3 = np.min(self.topcluster_stdphase3_singlesession_list[2])
            #minmin = min(min1,min2,min3)
            max1 = np.max(self.topcluster_stdphase1_singlesession_list[2])
            max2 = np.max(self.topcluster_stdphase2_singlesession_list[2])
            max3 = np.max(self.topcluster_stdphase3_singlesession_list[2])
            maxmax=max(max1,max2,max3)
            rects21 = ax[2].bar(x - width, h0, width, align='center', yerr=yerr0, ecolor='lightgrey')
            rects32 = ax[2].bar(x        , h1, width, align='center', yerr=yerr1, ecolor='lightgrey')
            rects33 = ax[2].bar(x + width, h2, width, align='center', yerr=yerr2, ecolor='lightgrey')
            ax[2].axhline(0, color='grey', linewidth=0.8)
            ax[2].xaxis.set_major_locator(ticker.FixedLocator(x))
            ax[2].xaxis.set_major_formatter(ticker.FixedFormatter([]))
            ax[2].yaxis.set_major_locator(ticker.FixedLocator([-maxmax,0,maxmax]))
            #ax[2].yaxis.set_major_formatter(ticker.FixedFormatter([]))
            #ax[2].legend()
            fig.tight_layout()
            plt.show()
        
        ######### Mann-Whitney Test results
        import statistics as stat
        import scipy.stats
        #all_plvs = np.zeros((7,3,3))
        fig00 = plt.figure("Mann-Whitney tests")
        alpha = 0.05
        counter = 1
        plv1_animalavg_list = []
        plv2_animalavg_list = []
        plv3_animalavg_list = []
        for k in range(3): ## loop over ana levels
            plv1_animalavg_list.append([])
            plv2_animalavg_list.append([])
            plv3_animalavg_list.append([])
        num_C1 = np.zeros(3,dtype=int)
        num_C2 = np.zeros(3,dtype=int)
        num_C3 = np.zeros(3,dtype=int)
        color = ['black','red','green','yellow','blue','cyan','magenta']
        maxy = -100000.0
        miny = 100000.0
        
        ##### choice of cluster to be analyzsed statistically
        ##### flag_topclusterdata = 0 : mutual_connections_not_data_list
        ##### flag_topclusterdata = 1 : MWTest_conds_list, e.g. nobasecluster
        flag_topclusterdata = 1 
        
        all_data_for_inference = []
        if flag_topclusterdata == 1:
            all_data_for_inference = self.MWTest_conds_list
        else:
            all_data_for_inference = self.mutual_connections_not_data_list
        self.compute_two_way_anova(all_data_for_inference)
        
        
        for analevel in range(3):
            for animal in range(7):
                if flag_topclusterdata == 1:
                    list_ = self.MWTest_conds_list[analevel][animal]
                else:
                    list_ = self.mutual_connections_not_data_list[analevel][animal]
                for k in range(3):
                    median = stat.median(list_[k])
                    if median>maxy:
                        maxy = median
                    if median<miny:
                        miny = median
                        
        
        
        for analevel in range(3):
            for animal in range(7):
                #list_ = self.MWTest_conds_list[analevel][animal]
                #list_ = self.mutual_connections_data_list[analevel][animal]
                #list_ = self.mutual_connections_not_data_list[analevel][animal]
                if flag_topclusterdata == 1:
                    list_ = self.MWTest_conds_list[analevel][animal]
                else:
                    list_ = self.mutual_connections_not_data_list[analevel][animal]
                num_C1[analevel] += np.size(list_[0])
                num_C2[analevel] += np.size(list_[1])
                num_C3[analevel] += np.size(list_[2])
                plv1_animalavg_list[analevel].append(list_[0])
                plv2_animalavg_list[analevel].append(list_[1])
                plv3_animalavg_list[analevel].append(list_[2])
                median1 = stat.median(list_[0])
                median2 = stat.median(list_[1])
                median3 = stat.median(list_[2])
                list_0 = [median1,median2,median3]
                #print("median1:"median1)
                res12 = scipy.stats.mannwhitneyu(list_[0], list_[1])
                res23 = scipy.stats.mannwhitneyu(list_[1], list_[2])
                res13 = scipy.stats.mannwhitneyu(list_[0], list_[2])
                sign_list = []
                if res12[1]<alpha:
                    print("ana:%d  animal:#%d    1-2  medians: %f  %f   p=%f"%(analevel,animal,median1,median2,res12[1]))
                    sign_list.append([1,2])
                if res23[1]<alpha:
                    print("ana:%d  animal:#%d    2-3  medians: %f  %f   p=%f"%(analevel,animal,median2,median3,res23[1]))
                    sign_list.append([2,3])
                if res13[1]<alpha:
                    print("ana:%d  animal:#%d    1-3  medians: %f  %f   p=%f"%(analevel,animal,median1,median3,res13[1]))
                    sign_list.append([1,3])
            
                ax = plt.subplot(3,7,counter)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                #plt.boxplot(list_,notch=True)
                plt.plot(np.arange(1,4,1),list_0,'ko-',color=color[animal])
                ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
                ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre','stim','post']))
                meds = list_0
                h = 0.08*(max(meds)-min(meds))
                y0 = maxy*1.15#max(meds)+15*h
                num_sign = len(sign_list)
                for kk in range(num_sign):
                    sign = sign_list[kk]
                    x1, x2 = sign[0],sign[1]
                    y = y0-kk*h*10
                    #print("kk=%d    y=%f  y0=%f   h=%f  maxy=%f   y+h=%f"%(kk,y,y0,h,max(meds),y+h))
                    #lab = '*'
                    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                    ax.text((x1+x2)*.5, y+h,'*', ha='center', va='bottom', color='k')
                plt.ylim(miny*0.97,maxy*1.2)   
                counter += 1
        fig00.subplots_adjust(hspace=0.3,wspace=0.3)
        
        data1_all_list = []#np.zeros((3,num_C1[0]))
        data2_all_list = []#np.zeros((3,num_C2[0]))
        data3_all_list = []#np.zeros((3,num_C3[0]))
        for k in range(3):
            data1_all_list.append([])
            data2_all_list.append([])
            data3_all_list.append([])
            
        fig01 = plt.figure("Inference avg animals")
        medians_list = []
        for analevel in range(3):
            for list_ in plv1_animalavg_list[analevel]: 
                for l in range(len(list_)):
                    data1_all_list[analevel].append(list_[l]) 
            counter = 0
            for list_ in plv2_animalavg_list[analevel]: 
                for l in range(len(list_)):
                    data2_all_list[analevel].append(list_[l]) 
            counter = 0
            for list_ in plv3_animalavg_list[analevel]: 
                for l in range(len(list_)):
                    data3_all_list[analevel].append(list_[l]) 
            median1 = stat.median(data1_all_list[analevel])
            median2 = stat.median(data2_all_list[analevel])
            median3 = stat.median(data3_all_list[analevel]) 
            medians_list.append([median1,median2,median3])
        maxy = max(max(medians_list))   
        miny = min(min(medians_list))   
        for analevel in range(3):
            counter = 0
            median1 = medians_list[analevel][0]
            median2 = medians_list[analevel][1]
            median3 = medians_list[analevel][2]
            res12 = scipy.stats.mannwhitneyu(data1_all_list[analevel],data2_all_list[analevel])
            res13 = scipy.stats.mannwhitneyu(data1_all_list[analevel],data3_all_list[analevel])
            res23 = scipy.stats.mannwhitneyu(data2_all_list[analevel],data3_all_list[analevel])
            sign_list = []
            if res12[1]<alpha:
                print("animalavg   ana:%d    1-2  medians: %f  %f   p=%f"%(\
                        analevel,median1,median2,res12[1]\
                    ))
                sign_list.append([1,2])
            else:
                print("animalavg   ana:%d    1-2   p=%f"%(\
                        analevel,res12[1]\
                    ))
            if res23[1]<alpha:
                print("animalavg   ana:%d    2-3  medians: %f  %f   p=%f"%(\
                        analevel,median2,median3,res23[1]\
                    ))
                sign_list.append([2,3])
            else:
                print("animalavg   ana:%d    2-3  p=%f"%(\
                        analevel,res23[1]\
                    ))
            if res13[1]<alpha:
                print("animalavg   ana:%d    1-3  medians: %f  %f   p=%f"%(\
                        analevel,median1,median3,res13[1]\
                    ))
                sign_list.append([1,3])
            else:
                print("animalavg   ana:%d    1-3   p=%f"%(\
                        analevel,res13[1]\
                    ))
            dat_list = [data1_all_list[analevel],data2_all_list[analevel],data3_all_list[analevel]]
            
            ax = plt.subplot(2,3,analevel+1)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.boxplot(dat_list,notch=True)
            ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre','stim','post']))
            if fmin==0.3 and fs==500.0:
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.4,0.7,1.0]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(['','','']))
            if fmin==30.0 and fs==500.0:
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.4,0.6,0.8]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(['','','']))
            #ax.xaxis.set_major_formatter(ticker.FixedFormatter(['','','']))
            #ax.yaxis.set_ticks([])
            #ax.xaxis.set_ticks([])
            
            ax = plt.subplot(2,3,analevel+1+3)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.plot(np.arange(1,4,1),[median1,median2,median3],'ko-')
            ax.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(['pre','stim','post']))
            if fmin==0.3 and fs==500.0:
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.8,0.85,0.9]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(['','','']))
            if fmin==30.0 and fs==500.0:
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.6,0.62,0.64]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(['','','']))
            #ax.yaxis.set_ticks([])
            #ax.xaxis.set_ticks([])
            
            
            ####### for <1Hz-band
            if fmin==0.1 and fs==200.0:
                plt.ylim(miny*0.95,maxy*1.08)
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.84,0.86,0.88]))
            ####### for delta-band
            if fmin==0.5 and fs==200.0:
                plt.ylim(0.84,0.90)
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.84,0.86,0.88]))
            ###### for gamma-band
            if fmin==30.0 and fs==200.0:
                plt.ylim(0.57,0.615)
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.57,0.59,0.61]))
            
            meds = [median1,median2,median3]
            h = 0.05*(max(meds)-min(meds))
            y0 = maxy+15*h
            num_sign = len(sign_list)
            for kk in range(num_sign):
                sign = sign_list[kk]
                x1, x2 = sign[0],sign[1]
                y = y0-kk*h*6
                print("kk=%d    y=%f  y0=%f   h=%f  maxy=%f   y+h=%f"%(kk,y,y0,h,max(meds),y+h))
                #lab = '*'
                ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax.text((x1+x2)*.5, y+h*0.1,'*', ha='center', va='bottom', color='k')
            plt.ylim(miny*0.98,maxy+18*h)   
            #ax.yaxis.set_major_formatter(ticker.FixedFormatter([0.84,'stim','post']))
        fig01.subplots_adjust(hspace=0.3,wspace=0.3)
    
        sign_list = []
        datalistlist = [data1_all_list,data2_all_list,data3_all_list]
        fig02 = plt.figure("avganimals stat on analevel")
        maxy = max(max(medians_list))
        miny = min(min(medians_list))
        for k in range(3):
            sign_list = []
            l_list = [medians_list[0][k],medians_list[1][k],medians_list[2][k]]
            datalist = datalistlist[k]
            res01 = scipy.stats.mannwhitneyu(datalist[0],datalist[1])
            res02 = scipy.stats.mannwhitneyu(datalist[0],datalist[2])
            res12 = scipy.stats.mannwhitneyu(datalist[1],datalist[2])
            if res01[1]<alpha:
                print("animalavg   cond:%d    low-medium  medians: %f  %f   p=%f"%(\
                        k,l_list[0],l_list[1],res01[1]\
                    ))
                sign_list.append([0,1])
            else:
                print("animalavg   cond:%d    low-medium   p=%f"%(\
                        k,res01[1]\
                    ))
            if res02[1]<alpha:
                print("animalavg   cond:%d    low-high  medians: %f  %f   p=%f"%(\
                        k,l_list[0],l_list[2],res02[1]\
                    ))
                sign_list.append([0,2])
            else:
                print("animalavg   cond:%d    low-high   p=%f"%(\
                        k,res02[1]\
                    ))
            if res12[1]<alpha:
                print("animalavg   cond:%d    medium-high  medians: %f  %f   p=%f"%(\
                        k,l_list[1],l_list[2],res12[1]\
                    ))
                sign_list.append([1,2])
            else:
                print("animalavg   cond:%d    medium-high   p=%f"%(\
                        k,res12[1]\
                    ))
                    
            ax = plt.subplot(2,3,k+1)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.plot(np.arange(0,3,1),l_list,'ko-')
            ax.xaxis.set_major_locator(ticker.FixedLocator([0,1,2]))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(['low','medium','high']))
            #ax.yaxis.set_ticks([])
            if fmin==0.3 and fs==500.0:
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.80,0.85,0.90]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(['','','']))
            if fmin==30.0 and fs==500.0:
                ax.yaxis.set_major_locator(ticker.FixedLocator([0.59,0.62,0.64]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(['','','']))
                
            meds = l_list
            h = 0.05*(max(meds)-min(meds))
            y0 = maxy*1.05#max(meds)+15*h
            num_sign = len(sign_list)
            for kk in range(num_sign):
                sign = sign_list[kk]
                x1, x2 = sign[0],sign[1]
                y = y0-kk*h*6
                print("kk=%d    y=%f  y0=%f   h=%f  maxy=%f   y+h=%f"%(kk,y,y0,h,max(meds),y+h))
                #lab = '*'
                ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax.text((x1+x2)*.5, y+h*0.1,'*', ha='center', va='bottom', color='k')
            plt.ylim(miny*0.95,maxy*1.08)
            
            #ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(min(meds)*0.9,y0*1.1,0.05)))
            #ax.xaxis.set_major_formatter(ticker.FixedFormatter(['low','medium','high']))
            
        fig02.subplots_adjust(hspace=0.3,wspace=0.3)
        
        #plt.show()

        for ana_ in range(3):
            num_sessions = len(self.connection_mutual_list[ana_])
            print("num_sessions = %d"%num_sessions)
            ## average top connections over all sessions
            for nsession in range(num_sessions):
                for k in range(self.num_chan):
                    for l in range(self.num_chan):
                        if self.connection_mutual_list[ana_][nsession][k,l]>0:
                                self.mutual_connections_allssesions_list[ana_][k,l] +=1
            for k in range(self.num_chan):
                for l in range(self.num_chan):
                    if self.mutual_connections_allssesions_list[ana_][k,l] == num_sessions:
                        self.mutual_connections_allssesions_list[ana_][k,l] = 1
                    else:
                        self.mutual_connections_allssesions_list[ana_][k,l] = 0
        
        mlabel = "./mutuals_%.1f-%.1f.dat"%(fmin,fmax)
        f = open(mlabel,"w+")
        str = ''
        for nsession in range(7):
            ## average top connections over ana conditions
            str = ''
            for ana_ in range(3):
                for k in range(self.num_chan):
                    for l in range(self.num_chan):
                        if self.connection_mutual_list[ana_][nsession][k,l]>0:
                            if flag_inference == 1:
                                pvals = self.connection_mutual_sig_list[ana_][nsession][1]
                                if pvals[0][k,l]==0 and pvals[1][k,l]==0 and pvals[2][k,l]==0:
                                    self.mutual_connections_allanaconds_list[nsession][k,l] +=1
                            else:
                                self.mutual_connections_allanaconds_list[nsession][k,l] +=1
            for k in range(self.num_chan):
                for l in range(self.num_chan):
                    if self.mutual_connections_allanaconds_list[nsession][k,l] == 3:
                        self.mutual_connections_allanaconds_list[nsession][k,l] = 1
                        if l>k:
                            str += ' %d-%d '%(k+1,l+1)
                    else:
                        self.mutual_connections_allanaconds_list[nsession][k,l] = 0
            str += '\n'
            f.write(str)
        f.close()
        
        if flag_plots==1:
            label = '_0'
            fig0 = plt.figure(label)
            for nsession0 in range(num_sessions):
                nsession = nsession0+1
                ax = plt.subplot(1,num_sessions,nsession)
                plt.imshow(self.mutual_connections_allanaconds_list[nsession0], interpolation='nearest', cmap='Reds',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.xticks([])
                plt.yticks([])
                
            label = '_1'
            fig0 = plt.figure(label)
            for ana_ in range(3):
                for nsession0 in range(num_sessions):
                    nsession = nsession0+1+ana_*num_sessions
                    ax = plt.subplot(3,num_sessions,nsession)
                    plt.imshow(self.connection_mutual_list[ana_][nsession0], interpolation='nearest', cmap='Reds',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                    plt.xticks([])
                    plt.yticks([])
                    #fig0.subplots_adjust(hspace=0.1,wspace=0.3)
                    
            label = '_2'
            fig1 = plt.figure(label)
            for ana_0 in range(3):
                ana_ = ana_0+1
                ax = plt.subplot(3,1,ana_)
                plt.imshow(self.mutual_connections_allssesions_list[ana_0],interpolation='nearest', cmap='Reds',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.xticks([])
                plt.yticks([])
                plt.colorbar()
        
            
            ## connect_mutual_num : [cond1-cond2,transition]
            ## cond1-cond2=0 : 1-2
            ## cond1-cond2=1 : 1-3
            ## cond1-cond2=2 : 2-3
            ## transition=0  : 1  1st is set, next is not set
            ## transition=1  : 2  1st is not set, next is set
            ## transition=2  : 3  1st and next are set
            connections = np.zeros((7,3,3),dtype=int)
            for animal in range(7):
                for analevel in range(3):
                    list_ = self.connection_mutual_transitions_list[analevel]
                    connection = list_[animal]
                    for stimcond in range(3):
                        ####### if connections > 0: next is more synchronized
                        ####### if connections < 0: 1st is more synchronized
                        connections[animal,analevel,stimcond]=connection[stimcond,1]-connection[stimcond,0]
            color = ['black','red','green','yellow','blue','cyan','magenta']
            marker = ['o','s','*']
            label = ['pre-stim','pre-post','stim-post']
            fig1 = plt.figure("1-2")
            for k in range(7): ## loop over animals
                kk = k+1
                ax = plt.subplot(1,7,kk)
                for l in range(3): ## loop over stim condition 
                    plt.plot(range(3),connections[k,:,l],marker=marker[0],color=color[l],label=label[l])
                plt.xticks([1,2,3])
                ax.legend()
            #plt.title("1-2")
            #ax = plt.subplot(1,3,2)
            #for l in range(7):
            #    plt.plot(range(3),connections[l,:,2],marker=marker[0],color=color[l])
            #plt.title("2-3")
            #ax = plt.subplot(1,3,3)
            #for l in range(7):
            #    plt.plot(range(3),connections[l,:,1],marker=marker[0],color=color[l])
            #plt.title("1-3")
        plt.show()
        
        
     
        
        
        import networkx as nx
        # Create a complete graph 
        nnodes = self.num_chan
        connectionmatrix_list = []
        G_list = []
        pos_list = []
        for nsession in range(7):
            connectionmatrix=self.mutual_connections_allanaconds_list[nsession]-np.eye(nnodes)
            connectionmatrix_list.append(connectionmatrix)   
            G = nx.from_numpy_array(connectionmatrix)
            G_list.append(G)
            pos_list.append(nx.circular_layout(G))
        
        node_labels0_1 = {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7, 7: 8, \
                        8: 9, 9: 10, 10: 11, 11: 12, 12: 13, 13: 14, 14: 15, 15: 16,}
        node_labels0_2 = {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7, 7: 8, \
                        8: 9, 9: 10, 10: 11, 11: 12, 12: 13, 13: 14, 14: 15, 15: 16,}
        node_color = ["tab:red","tab:orange","tab:olive","tab:green","tab:blue","tab:purple","tab:cyan"]
        
        labelcolor_list = []
        labelcolor_list.append('white')
        labelcolor_list.append('black')
        labelcolor_list.append('white')
        
        node_list = []
        
        for k in range(7):
            node_list.append([])
        #node_list0_1 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        if fmin==0.1:
            node_list[0].append([1,2,3,4,5,6,7,8,9,10,11,12,13,14])
            node_list[0].append([0,15])
            node_list[1].append([3,6,9,12])
            node_list[1].append([0,1,2,4,5,7,8,10,11,13,14,15])
            node_list[2].append([2,3,4,5,6,7,8,9,10,11,12])
            node_list[2].append([1,13,14])
            node_list[2].append([0,15])
            node_list[3].append([1,2,3,4,5,8,9,10,11,12,13])
            node_list[3].append([0,6,7,14,15])
            node_list[4].append([0,2,5,6,9])
            node_list[4].append([1,3,4,7,8,10,11,12,13,14,15])
            node_list[5].append([1,2,3,4,5,6,10,11,12])
            node_list[5].append([0,7,8,9,13,14,15])
            node_list[6].append([0,2,3,13,14,15])
            node_list[6].append([1,4,5,6,7,8,9,10,11,12])
        if fmin==0.3:
            node_list[0].append([0,1,2,5,6,7,8,9,10,11,12,13,14,15])
            node_list[0].append([3,4])
            node_list[1].append([3,6,9,10,11,12])
            node_list[1].append([0,1,2,4,5,7,8,13,14,15])
            node_list[2].append([2,3,4,5,7,8,9,10,11,12])
            node_list[2].append([0,1,6,13,14,15])
            node_list[3].append([1,2,3,4,5,8,9,10,11,12,13])
            node_list[3].append([0,6,7,14,15])
            node_list[4].append([1,2,3,4,5,6,9])
            node_list[4].append([0,7,8,10,11,12,13,14,15])
            node_list[5].append([1,2,3,4,5,6,9,10,11,12,13])
            node_list[5].append([0,7,8,14,15])
            node_list[6].append([0,1,2,3,4,13,14,15])
            node_list[6].append([5,6,7,8,9,10,11,12])
        if fmin==0.5:
            node_list[0].append([1,2,3,4,5,6,7,8,9,10,11,12,13,14])
            node_list[0].append([0,15])
            node_list[1].append([3,4,5,6,8,9,10,11,12,13])
            node_list[1].append([7,14,15])
            node_list[1].append([0,1,2])
            node_list[2].append([1,2,3,4,5,6,7,8,9,10,11,12])
            node_list[2].append([13,14])
            node_list[2].append([0,15])
            node_list[3].append([1,2,3,4,5,6,8,9,10,11,12,13,14])
            node_list[3].append([0,7,15])
            node_list[4].append([1,2,3,4,5,6,9,12,14])
            node_list[4].append([0,7,8,10,11,13,15])
            node_list[5].append([1,2,3,4,5,6,9,10,11,12])
            node_list[5].append([0,7,8,13,14,15])
            node_list[6].append([0,1,2,3,4,5,6,7,11,12,13,14,15])
            node_list[6].append([8,9,10])
        if fmin==30.0:
            node_list[0].append([3,4,5,6,7,8,9,10,11,12,13])
            node_list[0].append([1,2,14])
            node_list[0].append([0,15])
            node_list[1].append([2,3,4,5,6,7,8,9,10,11,12,13])
            node_list[1].append([0,1,14,15])
            node_list[2].append([3,4,5,6,8,9,10,11,12])
            node_list[2].append([1,2,7,13,14])
            node_list[2].append([0,15])
            node_list[3].append([2,3,4,5,6,7,8,9,10,11,12,13])
            node_list[3].append([0,1,14,15])
            node_list[4].append([1,2,3,4,5,6,9,10,11,12,13,14])
            node_list[4].append([0,7,8,15])
            node_list[5].append([2,3,4,5,6,7,9,10,11,12,13])
            node_list[5].append([0,1,8,14,15])
            node_list[6].append([1,2,3,4,5,6,7,9,10,11,12,13,14])
            node_list[6].append([0,8,15])
        # Create a figure with 1:1 aspect ratio to preserve the circle.
        if flag_plots==1:
            fig= plt.figure(figsize=(6, 10))
            for k in range(7):
                kk = k+1
                ax = plt.subplot(4,2,kk)
                node_opts_list = []
                node_opts_list.append({"node_size": 195, "node_color": node_color[k], "edgecolors": node_color[k], "linewidths": 1.0})
                node_opts_list.append({"node_size": 195, "node_color": "w", "edgecolors": node_color[k], "linewidths": 1.0})
                node_opts_list.append({"node_size": 195, "node_color": node_color[k], "edgecolors": node_color[k], "linewidths": 1.0})
                #nx.draw_networkx_nodes(G_list[k], pos_list[k], nodelist=[0, 1, 2, 3], **node_opts)
                for l in range(len(node_list[k])):
                    #print("l=%d"%l,"  ",labelcolor_list[l])
                    nx.draw_networkx_nodes(G_list[k], pos_list[k], nodelist=node_list[k][l],**node_opts_list[l])
                    nx.draw_networkx_labels(G_list[k], pos_list[k], font_size=10, labels=node_labels0_1, font_color="black")#labelcolor_list[l])#labelcolor[l])
                    nx.draw_networkx_edges(G_list[k], pos_list[k], width=1.0, edge_color=node_color[k], alpha=0.5)
                    #break
                ax.set_axis_off()
            fig.tight_layout()
            fig.subplots_adjust(hspace=0.06,wspace=0.1)
            plt.savefig("mutual_networks.png", dpi='figure', format='png')
        
        ########## plot networks of connections, which are not globally mutual
        
        connect_notmutual123_global_list = []
        for file0 in range(7):
            connect_1 = self.connect_notmutual123_list[0][file0]
            connect_2 = self.connect_notmutual123_list[1][file0]
            connect_3 = self.connect_notmutual123_list[2][file0]
            #print(connect_1)
            #print(connect_2)
            #print(connect_3)
            #quit()
            connect_notmutual123_global = np.eye(self.num_chan,dtype=int)
            for k1 in range(self.num_chan):
                for k20 in range(self.num_chan-k1-1):
                    k2 = k1+k20+1
                    #print("k1=%d k2=%d   %d %d %d"%(k1,k2,connect_1[k1,k2],connect_2[k1,k2],connect_3[k1,k2]))
                    if connect_1[k1,k2]==1 and \
                        connect_2[k1,k2]==1 and \
                        connect_3[k1,k2]==1:
                            connect_notmutual123_global[k1,k2]=1
                            connect_notmutual123_global[k2,k1]=1
                            #if file0==0:
                            #    print("file0:%d k1=%d k2=%d : YES!!!   %d %d %d "%(file0,k1,k2,connect_1[k1,k2],connect_2[k1,k2],connect_3[k1,k2]))
            connect_notmutual123_global_list.append(connect_notmutual123_global)
            
            #if file0==0:
            #    print(connect_notmutual123_global)
            #    print(connect_notmutual123_global_list[0])
            #quit()
            del connect_notmutual123_global
            #print("file0=%d"%file0,connect_notmutual123_global_list[0])
        #del connectionmatrix_list
        #del G
        #del G_list
        #del pos_list
        connectionmatrix_list = []
        G_list = []
        pos_list = []
        #print("test:",connect_notmutual123_global_list[0])
        for nsession in range(7):
            #print(connect_notmutual123_global_list[nsession])
            #quit()
            connectionmatrix=connect_notmutual123_global_list[nsession]-np.eye(nnodes)
            connectionmatrix_list.append(connectionmatrix)   
            G = nx.from_numpy_array(connectionmatrix)
            G_list.append(G)
            pos_list.append(nx.circular_layout(G))
        del node_list
        node_list = []
        for k in range(7):
            node_list.append([])
        if fmin==0.3:
            node_list[0].append([1,2,3,4,5,6,7,8,9,10,11,12,13,14])
            node_list[0].append([0,15])
            node_list[1].append([0,1,2,5,9,10,13,15])
            node_list[1].append([3,4,5,6,7,8,11,12,14])
            node_list[2].append([0,2,3,4,5,7,8,9,10,11,12,15])
            node_list[2].append([1,6,13,14])
            node_list[3].append([0,1,2,3,4,5,6,8,9,10,11,12,13,14])
            node_list[3].append([7,15])
            node_list[4].append([0,1,2,3,5,6,7,8,9,11,12,13,14])
            node_list[4].append([4,10,15])
            node_list[5].append([0,1,2,3,4,5,6,8,9,10,11,12,13,14,15])
            node_list[5].append([7])
            node_list[6].append([0,1,3,4,5,7,9,10,11,12,13,14,15])
            node_list[6].append([2,6,8])
        if fmin==30.0:
            node_list[0].append([2,3,4,5,6,7,8,9,10,11,12,13,14])
            node_list[0].append([0,1,15])
            node_list[1].append([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
            node_list[2].append([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14])
            node_list[2].append([15])
            node_list[3].append([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
            node_list[4].append([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14])
            node_list[4].append([15])
            node_list[5].append([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14])
            node_list[5].append([15])
            node_list[6].append([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
        if flag_plots==1:
            fig= plt.figure(figsize=(6, 10))
            for k in range(7):
                kk = k+1
                ax = plt.subplot(4,2,kk)
                node_opts_list = []
                node_opts_list.append({"node_size": 195, "node_color": node_color[k], "edgecolors": node_color[k], "linewidths": 1.0})
                node_opts_list.append({"node_size": 195, "node_color": "w", "edgecolors": node_color[k], "linewidths": 1.0})
                node_opts_list.append({"node_size": 195, "node_color": node_color[k], "edgecolors": node_color[k], "linewidths": 1.0})
                #nx.draw_networkx_nodes(G_list[k], pos_list[k], nodelist=[0, 1, 2, 3], **node_opts)
                for l in range(len(node_list[k])):
                    #print("l=%d"%l,"  ",labelcolor_list[l])
                    nx.draw_networkx_nodes(G_list[k], pos_list[k], nodelist=node_list[k][l],**node_opts_list[l])
                    nx.draw_networkx_labels(G_list[k], pos_list[k], font_size=10, labels=node_labels0_1, font_color="black")#labelcolor_list[l])#labelcolor[l])
                    nx.draw_networkx_edges(G_list[k], pos_list[k], width=1.0, edge_color=node_color[k], alpha=0.5)
                    #break
                ax.set_axis_off()
            fig.tight_layout()
            fig.subplots_adjust(hspace=0.06,wspace=0.1)
            plt.savefig("notmutual_networks.png", dpi='figure', format='png')
            
            plt.show()
        
        plt.figure("network maps: top cluster")
        counter = 1
        for nsession in range(7):
            ax = plt.subplot(1,7,nsession+1)
            self.mutual_connections_allanaconds_list[nsession]
            plt.imshow(self.mutual_connections_allanaconds_list[nsession],interpolation='nearest', cmap='Reds',origin='upper')
        plt.colorbar()
        plt.figure("network maps: second cluster")
        counter = 1
        for nsession in range(7):
            ax = plt.subplot(1,7,nsession+1)
            self.mutual_connections_allanaconds_list[nsession]
            plt.imshow(connect_notmutual123_global_list[nsession],interpolation='nearest', cmap='Reds',origin='upper')
        plt.colorbar()
        plt.show()
        
    def symmetry_check(self,connect,flag_plots):
        dim = np.shape(connect)[0]
        dim2 = dim//2
        index = np.zeros((dim,dim),dtype=int)
        counter = 0
        sum = 0.0
        for c in range(dim):
            k=0
            while 2*k+c<dim-1:
                l = dim-1-c-k
                if connect[k,k+c]==connect[l,l+c] and connect[k,k+c]>1:
                    sum += 1
                    #print("c=%d    k: %d %d   l:%d %d     connect: %d %d"%(c,k,k+c,l,l+c,connect[k,k+c],connect[l,l+c]))
                index[k,k+c]=c+1
                index[l,l+c]=c+1
                #print("c=%d    k: %d %d   l:%d %d"%(c,k,k+c,l,l+c))
                counter += 1
                k +=1
        print("sum=%d counter=%d"%(sum,counter))
        #plt.figure("symmetry check")
        #plt.imshow(connect, interpolation='nearest', cmap='Reds',origin='upper')
        #plt.show()
        return float(sum)/float(counter)
        #quit()
    
    def projection_on_networkinteractions(self,connect,flag_plots):
        dim = np.shape(connect)[0]
        index = np.zeros((dim,dim),dtype=int)
        list_proj12 = []
        
        list_proj1 = []
        kmax = dim
        #for k in range(kmax):
        #    k1 = k
        #    k2 = k
        #    if k1<dim//2 and k2<dim//2:
        #        list_proj1.append(connect[k1,k2])
        #        index[k1,k2]=1
        #    if k1>=dim//2 and k2>=dim//2:
        #        list_proj1.append(connect[k1,k2])
        #        index[k1,k2]=1
        for k in range(kmax-1):
            k1 = k
            k2 = k+1
            if k1<dim//2 and k2<dim//2:
                list_proj1.append(connect[k1,k2])
                index[k1,k2]=1
            if k1>=dim//2 and k2>=dim//2:
                list_proj1.append(connect[k1,k2])
                index[k1,k2]=1
            k1 = k+1
            k2 = k
            if k1<dim//2 and k2<dim//2:
                list_proj1.append(connect[k1,k2])
                index[k1,k2]=1
            if k1>=dim//2 and k2>=dim//2:
                list_proj1.append(connect[k1,k2])
                index[k1,k2]=1
            
        list_proj2 = []
        kmax = dim//2-1
        for k in range(kmax):
            k1 = k
            k2 = dim-k-1
            list_proj2.append(connect[k1,k2])
            index[k1,k2]=2
            index[k2,k1]=2
            #print("%d %d"%(k1,k2))
            k2 = dim-k-2
            list_proj2.append(connect[k1,k2])
            index[k1,k2]=2
            index[k2,k1]=2
            #print("%d %d"%(k1,k2))
            k1 = k+1
            k2 = dim-k-1
            list_proj2.append(connect[k1,k2])
            index[k1,k2]=2
            index[k2,k1]=2
            #print("%d %d"%(k1,k2))
        k1 = kmax
        k2 = dim-1-kmax
        list_proj2.append(connect[k1,k2])
        index[k1,k2]=2
        index[k2,k1]=2
        
        list_proj12 = list_proj1 + list_proj2
        
        list_proj3 = []
        for k in range(dim):
            for l in range(dim):
                if index[k,l]==0:
                    list_proj3.append(connect[k1,k2])
        
        import statistics as stat
        #proj1_mean = np.mean(list_proj1)
        #proj2_mean = np.mean(list_proj2)
        #proj3_mean = np.mean(list_proj3)
        #print("projection results: %f %f %f"%(proj1_mean,proj2_mean,proj3_mean))
        proj1_median = stat.median(list_proj1)
        proj2_median = stat.median(list_proj2)
        proj3_median = stat.median(list_proj3)
        proj12_median = stat.median(list_proj12)
        print("projection results: %f %f %f     %f"%(proj1_median,proj2_median,proj3_median,proj12_median))
        
        if flag_plots==1:
            plt.figure("index")
            ax = plt.subplot(1,2,1)
            plt.imshow(connect, interpolation='nearest', cmap='Reds',origin='upper')
            ax = plt.subplot(1,2,2)
            plt.imshow(index, interpolation='nearest', cmap='Reds',origin='upper')
            plt.show()
        
        #return [proj1_median,proj2_median,proj3_median,proj12_median,list_proj1,list_proj2,list_proj3]
        #return [list_proj2,list_proj3] # not bad, some analevels are distinguish interactions
        return [list_proj12,list_proj3]
    
    def compute_connectivity_statistics(self,path,ana,stype_list,file,flags,fmin,fmax,fs):
        flag_plots = flags[0]
        flag_inference = flags[1]
        #file='20110817B'
        print("file:",file)
        ana_index=-1
        if ana=='low':
            ana_index=0
        if ana=='medium':
            ana_index=1
        if ana=='high':
            ana_index=2
                
        def remove_PLV_outliers(C1,C2,C3,phase1,phase2,phase3):
            thr_percentage = 0.8#0.95
            data_uniform1,dens1,dens_uniform1 =self.clustering_KDE(C1) 
            data_uniform2,dens2,dens_uniform2 =self.clustering_KDE(C2) 
            data_uniform3,dens3,dens_uniform3 =self.clustering_KDE(C3)
            phase_data_uniform1,phase_dens1,phase_dens_uniform1 =self.clustering_phase_KDE(phase1,0.1) 
            phase_data_uniform2,phase_dens2,phase_dens_uniform2 =self.clustering_phase_KDE(phase2,0.1) 
            phase_data_uniform3,phase_dens3,phase_dens_uniform3 =self.clustering_phase_KDE(phase3,0.1)
            globalmax1 = np.argmax(dens_uniform1)
            max1 = np.max(dens_uniform1)
            min1 = np.min(dens_uniform1)
            thr1 = max1-thr_percentage*(max1-min1)
            max1_n = np.max(dens1)
            min1_n = np.min(dens1)
            thr1_n = max1_n-thr_percentage*(max1_n-min1_n) 
            globalmax2 = np.argmax(dens_uniform2)
            max2 = np.max(dens_uniform2)
            min2 = np.min(dens_uniform2)
            thr2 = max2-thr_percentage*(max2-min2)
            max2_n = np.max(dens2)
            min2_n = np.min(dens2)
            thr2_n = max2_n-thr_percentage*(max2_n-min2_n) 
            globalmax3 = np.argmax(dens_uniform3)
            max3 = np.max(dens_uniform3)
            min3 = np.min(dens_uniform3)
            thr3 = max3-thr_percentage*(max3-min3)
            max3_n = np.max(dens3)
            min3_n = np.min(dens3)
            thr3_n = max3_n-thr_percentage*(max3_n-min3_n) 
            
            dens_uniform1_list  = []
            data_uniform1_list  = []
            C1_list             = []
            C1_list_zero        = []
            dens1_list          = []
            C1_n_list           = []
            C1_n_list_zero      = []
            phase_dens1_list    = []
            phase1_n_list       = []
            dens_uniform2_list  = []
            data_uniform2_list  = []
            C2_list             = []
            C2_list_zero        = []
            dens2_list          = []
            C2_n_list           = []
            C2_n_list_zero      = []
            phase_dens2_list    = []
            phase2_n_list       = []
            dens_uniform3_list  = []
            data_uniform3_list  = []
            C3_list             = []
            C3_list_zero        = []
            dens3_list          = []
            C3_n_list           = []
            C3_n_list_zero      = []
            phase_dens3_list    = []
            phase3_n_list       = []
            chanpairs1_list     = []
            chanpairs2_list     = []
            chanpairs3_list     = []
            absthr = 0.6
            for k in range(len(C1)):
                #if C1[k]<absthr:
                #    continue
                if dens1[k]<thr1_n and C1[k]<C1[globalmax1]:
                    continue
                dens1_list.append(dens1[k]+0.0)
                C1_n_list.append(C1[k])
                phase_dens1_list.append(phase_dens1[k]+0.0)
                phase1_n_list.append(phase1[k])
                chanpairs1_list.append(k)
                #C1_n_list_zero.append(9.0)
            for k in range(len(C2)):
                #if C2[k]<absthr:
                #    continue
                if dens2[k]<thr2_n and C2[k]<C2[globalmax2]:
                    continue
                dens2_list.append(dens2[k]+0.0)
                C2_n_list.append(C2[k])
                phase_dens2_list.append(phase_dens2[k]+0.0)
                phase2_n_list.append(phase2[k])
                chanpairs2_list.append(k)
                #C2_n_list_zero.append(10.0)
            for k in range(len(C3)):
                #if C3[k]<absthr:
                #    continue
                if dens3[k]<thr3_n and C3[k]<C3[globalmax3]:
                    continue
                dens3_list.append(dens3[k]+0.0)
                C3_n_list.append(C3[k])
                phase_dens3_list.append(phase_dens3[k]+0.0)
                phase3_n_list.append(phase3[k])
                chanpairs3_list.append(k)
                    
            #print("condition 0: original number :",len(C1),"  reduced number :",len(C1_n_list)," with threshold ",thr1_n)
            #print("condition 1: original number :",len(C2),"  reduced number :",len(C2_n_list)," with threshold ",thr2_n)
            #print("condition 2: original number :",len(C3),"  reduced number :",len(C3_n_list)," with threshold ",thr3_n)
            
            if flag_plots == 1:
                from matplotlib import ticker
                fig = plt.figure(0)
                ax = plt.subplot(2,1,1)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.plot(data_uniform1,dens_uniform1,'k--',label='pre')
                plt.plot(data_uniform2,dens_uniform2,'r--',label='stim')
                plt.plot(data_uniform3,dens_uniform3,'b--',label='post')
                #ax.legend()
                plt.plot(C1_n_list,dens1_list,'k.',C2_n_list,dens2_list,'r.',C3_n_list,dens3_list,'b.')
                #ax.xaxis.set_major_locator(ticker.FixedLocator([0,0.25,0.5,0.75,1.0]))
                #ax.xaxis.set_major_formatter(ticker.FixedFormatter([]))
                #ax.yaxis.set_major_locator(ticker.FixedLocator([0,3]))
                #ax.yaxis.set_major_formatter(ticker.FixedFormatter([]))
                
                
                ax = plt.subplot(2,1,2)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.plot(phase_data_uniform1,phase_dens_uniform1,'k--',label='pre')
                plt.plot(phase_data_uniform2,phase_dens_uniform2,'r--',label='stim')
                plt.plot(phase_data_uniform3,phase_dens_uniform3,'b--',label='post')
                ax.legend()
                plt.plot(phase1_n_list,phase_dens1_list,'k.',\
                    (phase2_n_list),phase_dens2_list,'r.',\
                    (phase3_n_list),phase_dens3_list,'b.')
                #plt.plot(np.abs(phase1_n_list),phase_dens1_list,'k.',\
                #    np.abs(phase2_n_list),phase_dens2_list,'r.',\
                #    np.abs(phase3_n_list),phase_dens3_list,'b.')
                plt.xlim(-0.4,0.4)
                
                #ax.xaxis.set_major_locator(ticker.FixedLocator([-3.14,-1.57,0,1.57,3.14]))
                #ax.xaxis.set_major_formatter(ticker.FixedFormatter([]))
                #ax.yaxis.set_major_locator(ticker.FixedLocator([0,2,4,6]))
                #ax.yaxis.set_major_formatter(ticker.FixedFormatter([]))
                
                #ax.xaxis.set_ticks([])
                
                plt.show()
            
            result_list = []
            result_list.append([C1_n_list,C2_n_list,C3_n_list])
            result_list.append([dens1_list,dens2_list,dens3_list])
            result_list.append([chanpairs1_list,chanpairs2_list,chanpairs3_list])
            return result_list
        
        flag_clustering = 1 ## 0: kernel density estimation ; 1: hierarchical clustering
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
        
        label0 = '_%.1f-%.1f_'%(fmin,fmax)
        if fs==200.0:
            label=self.filelabel_connectivity+'_n_'+label0+file+'_'+ana+'.dat'
        if fs==500.0:
            label=self.filelabel_connectivity+'_n_'+label0+file+'_'+ana+'_fs500.0'+'.dat'
        
        print("read %s"%label)
        
        #label = self.filelabel_connectivity+'_'+file+'_'+ana+'.dat'
        C0 = np.loadtxt(label)
        self.num_chan = int(C0[-1,1])+1
        
        
        res_list=remove_PLV_outliers(C0[:,2],C0[:,4],C0[:,6],C0[:,8],C0[:,10],C0[:,12],)
        C1_list=res_list[0][0]
        C2_list=res_list[0][1]
        C3_list=res_list[0][2]
        dens1_list=res_list[1][0]
        dens2_list=res_list[1][1]
        dens3_list=res_list[1][2]
        chanpairs1_list=res_list[2][0]
        chanpairs2_list=res_list[2][1]
        chanpairs3_list=res_list[2][2]
        C1  = np.zeros(len(C1_list))
        C1  = C1_list[:]
        C2  = np.zeros(len(C2_list))
        C2  = C2_list[:]
        C3  = np.zeros(len(C3_list))
        C3  = C3_list[:]
        #C1 = C0[:,2]
        #C2 = C0[:,5]
        #C3 = C0[:,8]
        del res_list
        
        #import scipy.stats
        #import statistics as stat
        #res12 = scipy.stats.mannwhitneyu(C1, C2)
        #res23 = scipy.stats.mannwhitneyu(C2, C3)
        #res13 = scipy.stats.mannwhitneyu(C1, C3)
        #median1 = stat.median(C1)
        #median2 = stat.median(C2)
        #median3 = stat.median(C3)
        #print("C1-C2: %f - %f   "%(median1,median2),res12)
        #print("C2-C3: %f - %f   "%(median2,median3),res23)
        #print("C1-C3: %f - %f   "%(median1,median3),res13)
        
        
        
        res = []
        if flag_clustering == 0:
            data_u,kde,kde_uni =self.clustering_KDE(C1) 
            res.append([data_u,kde,kde_uni])
            data_u,kde,kde_uni =self.clustering_KDE(C2) 
            res.append([data_u,kde,kde_uni])
            data_u,kde,kde_uni =self.clustering_KDE(C3) 
            res.append([data_u,kde,kde_uni])
            fig00 = plt.figure(file)
            #plt.plot(C1,np.exp(res[0]),'ko',C2,np.exp(res[1]),'ro',C2,np.exp(res[2]),'bo')
            #print(np.shape(res[0][0]))
            #print(np.shape(res[0][1]))
            #print(np.shape(res[0][2]))
            zero1 = res[0][0]*0.0
            zero2 = res[1][0]*0.0+1.0
            zero3 = res[2][0]*0.0+2.0
            #plt.plot(C1,np.exp(res[0][1]),'ko',res[0][0],np.exp(res[0][2]),'k',C1,zero1,'k*')
            #plt.plot(C2,np.exp(res[1][1]),'ro',res[1][0],np.exp(res[1][2]),'r',C2,zero2,'r*')
            #plt.plot(C3,np.exp(res[2][1]),'bo',res[2][0],np.exp(res[2][2]),'b',C3,zero3,'b*')
            plt.plot(res[0][0],(res[0][2]),'k',C1,zero1,'k*')
            plt.plot(res[1][0],(res[1][2]),'r',C2,zero2,'r*')
            plt.plot(res[2][0],(res[2][2]),'b',C3,zero3,'b*')
            plt.show()
            quit()
        
        if flag_clustering == 1:
            channels1_list = []
            channels2_list = []
            channels3_list = []
            #print(np.shape(C1))
            
            connect1=np.eye(self.num_chan)
            connect2=np.eye(self.num_chan)
            connect3=np.eye(self.num_chan)
            connectphase1=np.zeros((self.num_chan,self.num_chan))
            connectphase2=np.zeros((self.num_chan,self.num_chan))
            connectphase3=np.zeros((self.num_chan,self.num_chan)) 
            
            for k in range(len(C1)):
                chan1 = int(C0[chanpairs1_list[k],0])
                chan2 = int(C0[chanpairs1_list[k],1])
                #chan1 = int(C0[k,0])
                #chan2 = int(C0[k,1])
                channels1_list.append([chan1,chan2])
                connect1[chan1,chan2]=C1[k]
                connect1[chan2,chan1]=C1[k]
                connectphase1[chan1,chan2]=C0[k,8]
                connectphase1[chan2,chan1]=-C0[k,8]

            for k in range(len(C2)):
                chan1 = int(C0[chanpairs2_list[k],0])
                chan2 = int(C0[chanpairs2_list[k],1])
                channels2_list.append([chan1,chan2])
                connect2[chan1,chan2]=C2[k]
                connect2[chan2,chan1]=C2[k]
                connectphase2[chan1,chan2]=C0[k,10]
                connectphase2[chan2,chan1]=-C0[k,10]
                
            for k in range(len(C3)):
                chan1 = int(C0[chanpairs3_list[k],0])
                chan2 = int(C0[chanpairs3_list[k],1])
                channels3_list.append([chan1,chan2])    
                connect3[chan1,chan2]=C3[k]
                connect3[chan2,chan1]=C3[k]
                connectphase3[chan1,chan2]=C0[k,12]
                connectphase3[chan2,chan1]=-C0[k,12]
            
            
            
            ##### compute projections on nearest-neighbor networks
            projections1_list = self.projection_on_networkinteractions(connect1,flag_plots)
            projections2_list = self.projection_on_networkinteractions(connect2,flag_plots)
            projections3_list = self.projection_on_networkinteractions(connect3,flag_plots)
            results1 = []
            results1.append(projections1_list[0])
            results1.append(projections2_list[0])
            results1.append(projections3_list[0]) # nearest neighbor
            results2 = []
            results2.append(projections1_list[1])
            results2.append(projections2_list[1])
            results2.append(projections3_list[1]) # rest
            #results3 = [projections1_list[2],projections2_list[2],projections3_list[2]]
            if flag_plots==2:
                plt.figure("projections")
                ax = plt.subplot(1,1,1)
                plt.plot(range(3),results1,'ko-',label='proj 1')
                plt.plot(range(3),results2,'ro-',label='proj 2')
                #plt.plot(range(3),results3,'bo-',label='proj 3')
                ax.legend(loc=1,fontsize=8)
                plt.show()
            list_foranova = []
            #list_foranova.append([projections1_list[4],projections1_list[5],projections1_list[6]])
            #list_foranova.append([projections2_list[4],projections2_list[5],projections2_list[6]])
            #list_foranova.append([projections3_list[4],projections3_list[5],projections3_list[6]])
            list_foranova.append(results1)
            list_foranova.append(results2)
            self.projections_list[ana_index].append(list_foranova)
            #print(len(channels1_list),channels1_list)
            #quit()
            ##### end of : compute projections on nearest-neighbor networks
            
            
            num_cluster1,clusters1 =self.clustering_hierarchical(C1,3) 
            #res.append(clusters1)
            num_cluster2,clusters2 =self.clustering_hierarchical(C2,3) 
            #res.append(clusters2)
            num_cluster3,clusters3 =self.clustering_hierarchical(C3,3) 
            #res.append(clusters3)
            
            connect1C=np.eye(self.num_chan,dtype=int)*0
            connect2C=np.eye(self.num_chan,dtype=int)*0
            connect3C=np.eye(self.num_chan,dtype=int)*0
            cluster_mean1sync = np.zeros(num_cluster1)
            cluster_mean2sync = np.zeros(num_cluster2)
            cluster_mean3sync = np.zeros(num_cluster3)
            cluster_meanphase1 = np.zeros(num_cluster1)
            cluster_meanphase2 = np.zeros(num_cluster2)
            cluster_meanphase3 = np.zeros(num_cluster3)
            cluster_minmax1sync = np.zeros((num_cluster1,2))
            cluster_minmax1sync[:,0]=100000.0
            cluster_minmax1sync[:,1]=-100000.0
            cluster_minmaxphase1 = np.zeros((num_cluster1,2))
            cluster_minmaxphase1[:,0]=100000.0
            cluster_minmaxphase1[:,1]=-100000.0
            cluster_minmax2sync = np.zeros((num_cluster2,2))
            cluster_minmax2sync[:,0]=100000.0
            cluster_minmax2sync[:,1]=-100000.0
            cluster_minmaxphase2 = np.zeros((num_cluster2,2))
            cluster_minmaxphase2[:,0]=100000.0
            cluster_minmaxphase2[:,1]=-100000.0
            cluster_minmax3sync = np.zeros((num_cluster3,2))
            cluster_minmax3sync[:,0]=100000.0
            cluster_minmax3sync[:,1]=-100000.0
            cluster_minmaxphase3 = np.zeros((num_cluster3,2))
            cluster_minmaxphase3[:,0]=100000.0
            cluster_minmaxphase3[:,1]=-100000.0
            counter_clusters1 = np.zeros(num_cluster1,dtype=int)
            counter_clusters2 = np.zeros(num_cluster2,dtype=int)
            counter_clusters3 = np.zeros(num_cluster3,dtype=int)
            
            for k in range(len(C1)):
                chan1 = channels1_list[k][0]
                chan2 = channels1_list[k][1]
                connect1C[chan1,chan2]=clusters1[k]
                connect1C[chan2,chan1]=clusters1[k]
                counter_clusters1[clusters1[k]-1] +=1
                cluster_mean1sync[clusters1[k]-1] += C1[k]
                cluster_meanphase1[clusters1[k]-1] += connectphase1[chan1,chan2]
                if C1[k]<cluster_minmax1sync[clusters1[k]-1,0]:
                    cluster_minmax1sync[clusters1[k]-1,0]=C1[k]
                if C1[k]>=cluster_minmax1sync[clusters1[k]-1,1]:
                    cluster_minmax1sync[clusters1[k]-1,1]=C1[k]
                if connectphase1[chan1,chan2]<cluster_minmaxphase1[clusters1[k]-1,0]:
                    cluster_minmaxphase1[clusters1[k]-1,0]=connectphase1[chan1,chan2]
                if connectphase1[chan1,chan2]>=cluster_minmaxphase1[clusters1[k]-1,1]:
                    cluster_minmaxphase1[clusters1[k]-1,1]=connectphase1[chan1,chan2]
            for l in range(num_cluster1):
                cluster_mean1sync[l] /= float(counter_clusters1[l])
                cluster_meanphase1[l] /= float(counter_clusters1[l])
                
            for k in range(len(C2)):
                chan1 = channels2_list[k][0]
                chan2 = channels2_list[k][1]
                connect2C[chan1,chan2]=clusters2[k]
                connect2C[chan2,chan1]=clusters2[k]
                counter_clusters2[clusters2[k]-1] +=1
                cluster_mean2sync[clusters2[k]-1] += C2[k]
                cluster_meanphase2[clusters2[k]-1] += connectphase2[chan1,chan2]
                if C2[k]<cluster_minmax2sync[clusters2[k]-1,0]:
                    cluster_minmax2sync[clusters2[k]-1,0]=C2[k]
                if C2[k]>=cluster_minmax2sync[clusters2[k]-1,1]:
                    cluster_minmax2sync[clusters2[k]-1,1]=C2[k]
                if connectphase2[chan1,chan2]<cluster_minmaxphase2[clusters2[k]-1,0]:
                    cluster_minmaxphase2[clusters2[k]-1,0]=connectphase2[chan1,chan2]
                if connectphase2[chan1,chan2]>=cluster_minmaxphase2[clusters2[k]-1,1]:
                    cluster_minmaxphase2[clusters2[k]-1,1]=connectphase2[chan1,chan2]
            for l in range(num_cluster2):
                cluster_mean2sync[l] /= float(counter_clusters2[l])
                cluster_meanphase2[l] /= float(counter_clusters2[l])
                
            for k in range(len(C3)):
                chan1 = channels3_list[k][0]
                chan2 = channels3_list[k][1]
                connect3C[chan1,chan2]=clusters3[k]
                connect3C[chan2,chan1]=clusters3[k]
                counter_clusters3[clusters3[k]-1] +=1
                cluster_mean3sync[clusters3[k]-1] += C3[k]
                cluster_meanphase3[clusters3[k]-1] += connectphase3[chan1,chan2]
                if C3[k]<cluster_minmax3sync[clusters3[k]-1,0]:
                    cluster_minmax3sync[clusters3[k]-1,0]=C3[k]
                if C3[k]>=cluster_minmax3sync[clusters3[k]-1,1]:
                    cluster_minmax3sync[clusters3[k]-1,1]=C3[k]
                if connectphase3[chan1,chan2]<cluster_minmaxphase3[clusters3[k]-1,0]:
                    cluster_minmaxphase3[clusters3[k]-1,0]=connectphase3[chan1,chan2]
                if connectphase3[chan1,chan2]>=cluster_minmaxphase3[clusters3[k]-1,1]:
                    cluster_minmaxphase3[clusters3[k]-1,1]=connectphase3[chan1,chan2]                    
            for l in range(num_cluster3):
                cluster_mean3sync[l] /= float(counter_clusters3[l])
                cluster_meanphase3[l] /= float(counter_clusters3[l])
                
              
            final_clusters1_arglist = np.argsort(cluster_mean1sync,axis=0)# provides the ascending list of clusters
            final_clusters2_arglist = np.argsort(cluster_mean2sync,axis=0)# provides the ascending list of clusters
            final_clusters3_arglist = np.argsort(cluster_mean3sync,axis=0)# provides the ascending list of clusters
            #print(cluster_meansync)
            #print("final_clusters_arglist:",final_clusters_arglist)
            ## final_clusters_arglist[i,j]: i is condition, new sorted cluster j is cluster final_clusters_arglist[i,j]
            ## final_clusters_arglist_inverse[i,j]: old cluster j has the new place final_clusters_arglist_inverse[i,j]
            final_clusters1_arglist_inverse = final_clusters1_arglist*0
            final_clusters2_arglist_inverse = final_clusters2_arglist*0
            final_clusters3_arglist_inverse = final_clusters3_arglist*0
            for cc in range(num_cluster1):
                for s in range(num_cluster1):
                    if final_clusters1_arglist[s]==cc:
                        final_clusters1_arglist_inverse[cc]=int(s)
            for cc in range(num_cluster2):
                for s in range(num_cluster2):
                    if final_clusters2_arglist[s]==cc:
                        final_clusters2_arglist_inverse[cc]=int(s)
            for cc in range(num_cluster3):
                for s in range(num_cluster3):
                    if final_clusters3_arglist[s]==cc:
                        final_clusters3_arglist_inverse[cc]=int(s)
                        
            #print("final_clusters_arglist_inverse:",final_clusters_arglist_inverse)
            #quit()
            cluster_mean1sync_sorted = np.zeros((num_cluster1))
            cluster_mean2sync_sorted = np.zeros((num_cluster2))
            cluster_mean3sync_sorted = np.zeros((num_cluster3))
            cluster_minmax1sync_sorted = np.zeros((num_cluster1,2))
            cluster_minmax2sync_sorted = np.zeros((num_cluster2,2))
            cluster_minmax3sync_sorted = np.zeros((num_cluster3,2))
            cluster_meanphase1_sorted = np.zeros((num_cluster1))
            cluster_meanphase2_sorted = np.zeros((num_cluster2))
            cluster_meanphase3_sorted = np.zeros((num_cluster3))
            cluster_minmaxphase1_sorted = np.zeros((num_cluster1,2))
            cluster_minmaxphase2_sorted = np.zeros((num_cluster2,2))
            cluster_minmaxphase3_sorted = np.zeros((num_cluster3,2))
            connect1C_sorted        = np.eye(self.num_chan,dtype=int)*0
            connect2C_sorted        = np.eye(self.num_chan,dtype=int)*0
            connect3C_sorted        = np.eye(self.num_chan,dtype=int)*0
            #channels_cluster_all    = []
            connections_cluster_all = []
            cluster1_channels = np.zeros((num_cluster1,len(C1)),dtype=int)
            cluster2_channels = np.zeros((num_cluster2,len(C2)),dtype=int)
            cluster3_channels = np.zeros((num_cluster3,len(C3)),dtype=int)
            
            #if condition > 0:
            #    break
            for cc in range(num_cluster1):
                c = cc+1
                cluster_mean1sync_sorted[cc]        = cluster_mean1sync[final_clusters1_arglist[cc]]
                cluster_minmax1sync_sorted[cc,0]    = cluster_minmax1sync[final_clusters1_arglist[cc],0]
                cluster_minmax1sync_sorted[cc,1]    = cluster_minmax1sync[final_clusters1_arglist[cc],1]
                cluster_meanphase1_sorted[cc]       = cluster_meanphase1[final_clusters1_arglist[cc]]
                cluster_minmaxphase1_sorted[cc,0]   = cluster_minmaxphase1[final_clusters1_arglist[cc],0]
                cluster_minmaxphase1_sorted[cc,1]   = cluster_minmaxphase1[final_clusters1_arglist[cc],1]
                
            for cc in range(num_cluster2):
                c = cc+1
                cluster_mean2sync_sorted[cc]=cluster_mean2sync[final_clusters2_arglist[cc]]
                cluster_minmax2sync_sorted[cc,0]=cluster_minmax2sync[final_clusters2_arglist[cc],0]
                cluster_minmax2sync_sorted[cc,1]=cluster_minmax2sync[final_clusters2_arglist[cc],1]
                cluster_meanphase2_sorted[cc]       = cluster_meanphase2[final_clusters2_arglist[cc]]
                cluster_minmaxphase2_sorted[cc,0]   = cluster_minmaxphase2[final_clusters2_arglist[cc],0]
                cluster_minmaxphase2_sorted[cc,1]   = cluster_minmaxphase2[final_clusters2_arglist[cc],1]
            for cc in range(num_cluster3):
                c = cc+1
                cluster_mean3sync_sorted[cc]=cluster_mean3sync[final_clusters3_arglist[cc]]
                cluster_minmax3sync_sorted[cc,0]=cluster_minmax3sync[final_clusters3_arglist[cc],0]
                cluster_minmax3sync_sorted[cc,1]=cluster_minmax3sync[final_clusters3_arglist[cc],1]
                cluster_meanphase3_sorted[cc]       = cluster_meanphase3[final_clusters3_arglist[cc]]
                cluster_minmaxphase3_sorted[cc,0]   = cluster_minmaxphase3[final_clusters3_arglist[cc],0]
                cluster_minmaxphase3_sorted[cc,1]   = cluster_minmaxphase3[final_clusters3_arglist[cc],1]
            
            if flag_plots == 1:
                for c in range(num_cluster1):
                    help = counter_clusters1[final_clusters1_arglist[c]]
                    #print("condition 0    #%d : mean sync=%f     min sync=%f max sync=%f    number=%d"%(c,cluster_mean1sync_sorted[c],cluster_minmax1sync_sorted[c,0],cluster_minmax1sync_sorted[c,1],help))
                for c in range(num_cluster2):
                    help = counter_clusters2[final_clusters2_arglist[c]]
                    #print("condition 1    #%d : mean sync=%f     min sync=%f max sync=%f    number=%d"%(c,cluster_mean2sync_sorted[c],cluster_minmax2sync_sorted[c,0],cluster_minmax2sync_sorted[c,1],help))
                for c in range(num_cluster3):
                    help = counter_clusters3[final_clusters3_arglist[c]]
                    #print("condition 2    #%d : mean sync=%f     min sync=%f max sync=%f    number=%d"%(c,cluster_mean3sync_sorted[c],cluster_minmax3sync_sorted[c,0],cluster_minmax3sync_sorted[c,1],help))    
            
                for c in range(num_cluster1):
                    help = counter_clusters1[final_clusters1_arglist[c]]
                    #print("condition 0    #%d : mean phase=%f     min phase=%f max phase=%f    number=%d"%(c,cluster_meanphase1_sorted[c],cluster_minmaxphase1_sorted[c,0],cluster_minmaxphase1_sorted[c,1],help))
                for c in range(num_cluster2):
                    help = counter_clusters2[final_clusters2_arglist[c]]
                    #print("condition 1    #%d : mean phase=%f     min phase=%f max phase=%f    number=%d"%(c,cluster_meanphase2_sorted[c],cluster_minmaxphase2_sorted[c,0],cluster_minmaxphase2_sorted[c,1],help))
                for c in range(num_cluster3):
                    help = counter_clusters3[final_clusters3_arglist[c]]
                    #print("condition 2    #%d : mean phase=%f     min phase=%f max phase=%f    number=%d"%(c,cluster_meanphase3_sorted[c],cluster_minmaxphase3_sorted[c,0],cluster_minmaxphase3_sorted[c,1],help))    
                    
            if flag_plots == 1:
                        
                #------------------
                from matplotlib import ticker
                fig, ax = plt.subplots(3,2)
                #ax[0].axis('off')
                #ax[1].axis('off')
                #ax[2].axis('off')
                for k in range(3):
                    for l in range(2):
                        ax[k,l].spines['top'].set_visible(False)
                        ax[k,l].spines['right'].set_visible(False)
                #ax[1].spines['top'].set_visible(False)
                #ax[1].spines['right'].set_visible(False)
                #ax[2].spines['top'].set_visible(False)
                #ax[2].spines['right'].set_visible(False)
                
                #cluster_mean1sync_sorted[cc]=cluster_mean1sync[final_clusters1_arglist[cc]]
                #cluster_minmax1sync_sorted[cc,0]=cluster_minmax1sync[final_clusters1_arglist[cc],0]
                #cluster_minmax1sync_sorted[cc,1]
                
                h0=cluster_mean1sync_sorted
                yerr0=np.zeros((2,3))
                yerr0[0,:]=h0-cluster_minmax1sync_sorted[:,0]
                yerr0[1,:]=cluster_minmax1sync_sorted[:,1]-h0
                x = np.linspace(1,3,3)
                labels = ['cluster 1','cluster 2','cluster 3']
                width = 0.5
                ax[0,0].bar(x, h0, width, align='center', yerr=yerr0, color=['c','y','brown'], ecolor='lightgrey')
                #ax[0,0].axhline(0, color='grey', linewidth=0.8)
                ax[0,0].xaxis.set_major_locator(ticker.FixedLocator(x))
                ax[0,0].xaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[0,0].yaxis.set_major_locator(ticker.FixedLocator([0,0.5,1.0]))
                #ax[0,0].yaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[0,0].xaxis.set_ticks([])
                
                h0=cluster_meanphase1_sorted
                #yerr0=np.zeros((2,3))
                #yerr0[0,:]=h0-cluster_minmaxphase1_sorted[:,0]
                #yerr0[1,:]=cluster_minmaxphase1_sorted[:,1]-h0
                yerr0[0,:]=np.sqrt(-2.0*np.log(cluster_mean1sync_sorted))
                yerr0[1,:]=np.sqrt(-2.0*np.log(cluster_mean1sync_sorted))
                #minmin=np.min(cluster_minmaxphase1_sorted[:,0])
                #maxmax=np.max(cluster_minmaxphase1_sorted[:,1])
                maxmax=np.max(yerr0[1,:]+h0)
                ax[0,1].bar(x, h0, width, align='center', yerr=yerr0, color=['c','y','brown'], ecolor='lightgrey')
                ax[0,1].axhline(0, color='grey', linewidth=0.8)
                ax[0,1].xaxis.set_major_locator(ticker.FixedLocator(x))
                ax[0,1].xaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[0,1].yaxis.set_major_locator(ticker.FixedLocator([-maxmax,0.0,maxmax]))
                #ax[0,1].yaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[0,1].xaxis.set_ticks([])
                #print("maxmax=%f"%maxmax)
                #print("yerr0",yerr0[0,:],"  ",yerr0[1,:])
                #print("h0",h0)
                
                h0=cluster_mean2sync_sorted
                yerr0[0,:]=h0-cluster_minmax2sync_sorted[:,0]
                yerr0[1,:]=cluster_minmax2sync_sorted[:,1]-h0
                ax[1,0].bar(x, h0, width, align='center', yerr=yerr0, color=['c','y','brown'], ecolor='lightgrey')
                ax[1,0].xaxis.set_major_locator(ticker.FixedLocator(x))
                ax[1,0].xaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[1,0].yaxis.set_major_locator(ticker.FixedLocator([0,0.5,1.0]))
                #ax[1,0].yaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[1,0].xaxis.set_ticks([])
                
                h0=cluster_meanphase2_sorted
                #yerr0[0,:]=h0-cluster_minmaxphase2_sorted[:,0]
                #yerr0[1,:]=cluster_minmaxphase2_sorted[:,1]-h0
                yerr0[0,:]=np.sqrt(-2.0*np.log(cluster_mean2sync_sorted))
                yerr0[1,:]=np.sqrt(-2.0*np.log(cluster_mean2sync_sorted))
                maxmax=np.max(yerr0[1,:]+h0)
                minmin=-maxmax
                #minmin=np.min(cluster_minmaxphase1_sorted[:,0])
                #maxmax=np.max(cluster_minmaxphase1_sorted[:,1])
                #minmin=np.min(cluster_minmaxphase2_sorted[:,0])
                #maxmax=np.max(cluster_minmaxphase2_sorted[:,1])
                ax[1,1].bar(x, h0, width, align='center', yerr=yerr0, color=['c','y','brown'], ecolor='lightgrey')
                ax[1,1].axhline(0, color='grey', linewidth=0.8)
                ax[1,1].xaxis.set_major_locator(ticker.FixedLocator(x))
                ax[1,1].xaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[1,1].yaxis.set_major_locator(ticker.FixedLocator([minmin,0.0,maxmax]))
                #ax[1,1].yaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[1,1].xaxis.set_ticks([])
                
        
                h0=cluster_mean3sync_sorted
                yerr0[0,:]=h0-cluster_minmax3sync_sorted[:,0]
                yerr0[1,:]=cluster_minmax3sync_sorted[:,1]-h0
                ax[2,0].bar(x, h0, width, align='center', yerr=yerr0, color=['c','y','brown'], ecolor='lightgrey')
                ax[2,0].xaxis.set_major_locator(ticker.FixedLocator(x))
                ax[2,0].xaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[2,0].yaxis.set_major_locator(ticker.FixedLocator([0,0.5,1.0]))
                #ax[2,0].yaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[2,0].xaxis.set_ticks([])
                
                h0=cluster_meanphase3_sorted
                #yerr0[0,:]=h0-cluster_minmaxphase3_sorted[:,0]
                #yerr0[1,:]=cluster_minmaxphase3_sorted[:,1]-h0
                yerr0[0,:]=np.sqrt(-2.0*np.log(cluster_mean3sync_sorted))
                yerr0[1,:]=np.sqrt(-2.0*np.log(cluster_mean3sync_sorted))
                maxmax=np.max(yerr0[1,:]+h0)
                minmin=-maxmax
                #minmin=np.min(cluster_minmaxphase3_sorted[:,0])
                #maxmax=np.max(cluster_minmaxphase3_sorted[:,1])
                ax[2,1].bar(x, h0, width, align='center', yerr=yerr0, color=['c','y','brown'], ecolor='lightgrey')
                ax[2,1].axhline(0, color='grey', linewidth=0.8)
                ax[2,1].xaxis.set_major_locator(ticker.FixedLocator(x))
                ax[2,1].xaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[2,1].yaxis.set_major_locator(ticker.FixedLocator([minmin,0.0,maxmax]))
                #ax[2,1].yaxis.set_major_formatter(ticker.FixedFormatter([]))
                ax[2,1].xaxis.set_ticks([])
                
                fig.subplots_adjust(hspace=0.3,wspace=0.4)
    
                plt.show()
                #quit()
                #------------------
            
            from scipy.stats import circstd
            connections_cluster1 = []
            connections_cluster2 = []
            connections_cluster3 = []
            C1_nobase = []
            C2_nobase = []
            C3_nobase = []
            C1_top = []
            C2_top = []
            C3_top = []
            C1_medium = []
            C2_medium = []
            C3_medium = []
            
            #connect1C_sorted_nobase = np.eye(self.num_chan,dtype=int)
            num_topclusterchannels1 = 0
            topcluster_meanphase1 = 0.0
            phases_list = []
            for k in range(len(C1)):
                chan1 = channels1_list[k][0]
                chan2 = channels1_list[k][1]
                cluster = final_clusters1_arglist_inverse[int(connect1C[chan1,chan2])-1]+1
                connect1C_sorted[chan1,chan2]=cluster
                connect1C_sorted[chan2,chan1]=connect1C_sorted[chan1,chan2]
                if cluster>1:
                    #connect1C_sorted_nobase[chan1,chan2]=1 
                    C1_nobase.append(C1[k])
                    #if C1[k]<0.5:
                    #    print("*********************** ana=%d  animal: %s  C1[%d]=%f   cluster=%d"%(ana_index,file,k,C1[k],cluster))
                if cluster==num_cluster1:
                    C1_top.append(C1[k])
                if cluster==num_cluster1-1:
                    C1_medium.append(C1[k])
                cluster1_channels[cluster-1,k]=1
                connections_cluster1.append(cluster)
                if cluster == num_cluster1:
                    topcluster_meanphase1 += connectphase1[chan1,chan2]
                    num_topclusterchannels1 += 1
                    self.topcluster_meanphase1_sessionpool_list.append(connectphase1[chan1,chan2])
                    phases_list.append(connectphase1[chan1,chan2])
            topcluster_meanphase1 /= float(num_topclusterchannels1)
            self.topcluster_meanphase1_singlesession_list[ana_index].append(topcluster_meanphase1)
            topcluster_stdphase1 = circstd(phases_list)
            #print("topcluster_stdphase1",topcluster_stdphase1,"  phases_list:",phases_list)
            self.topcluster_stdphase1_singlesession_list[ana_index].append(topcluster_stdphase1)
            del phases_list
            
            #connect2C_sorted_nobase = np.eye(self.num_chan,dtype=int)
            num_topclusterchannels2 = 0
            topcluster_meanphase2 = 0.0
            phases_list = []
            for k in range(len(C2)):
                chan1 = channels2_list[k][0]
                chan2 = channels2_list[k][1]
                cluster = final_clusters2_arglist_inverse[int(connect2C[chan1,chan2])-1]+1
                connect2C_sorted[chan1,chan2]=cluster
                connect2C_sorted[chan2,chan1]=connect2C_sorted[chan1,chan2]
                if cluster>1:
                    #connect2C_sorted_nobase[chan1,chan2]=1 
                    C2_nobase.append(C2[k])
                if cluster==num_cluster2:
                    C2_top.append(C2[k])
                if cluster==num_cluster2-1:
                    C2_medium.append(C2[k])
                cluster2_channels[cluster-1,k]=1
                connections_cluster2.append(cluster)
                if cluster == num_cluster2:
                    topcluster_meanphase2 += connectphase2[chan1,chan2]
                    num_topclusterchannels2 += 1
                    self.topcluster_meanphase2_sessionpool_list.append(connectphase2[chan1,chan2])
                    phases_list.append(connectphase2[chan1,chan2])
            topcluster_meanphase2 /= float(num_topclusterchannels2)
            self.topcluster_meanphase2_singlesession_list[ana_index].append(topcluster_meanphase2)
            topcluster_stdphase2 = circstd(phases_list)
            self.topcluster_stdphase2_singlesession_list[ana_index].append(topcluster_stdphase2)
            del phases_list
            
            #connect3C_sorted_nobase = np.eye(self.num_chan,dtype=int)
            num_topclusterchannels3 = 0
            topcluster_meanphase3 = 0.0
            phases_list = []
            for k in range(len(C3)):
                chan1 = channels3_list[k][0]
                chan2 = channels3_list[k][1]
                cluster = final_clusters3_arglist_inverse[int(connect3C[chan1,chan2])-1]+1
                connect3C_sorted[chan1,chan2]=cluster
                connect3C_sorted[chan2,chan1]=connect3C_sorted[chan1,chan2]
                if cluster>1:
                    #connect3C_sorted_nobase[chan1,chan2]=1
                    C3_nobase.append(C3[k])
                if cluster==num_cluster3:
                    C3_top.append(C3[k])
                if cluster==num_cluster3-1:
                    C3_medium.append(C3[k])
                cluster3_channels[cluster-1,k]=1
                connections_cluster3.append(cluster)
                if cluster == num_cluster3:
                    topcluster_meanphase3 += connectphase3[chan1,chan2]
                    num_topclusterchannels3 += 1
                    self.topcluster_meanphase3_sessionpool_list.append(connectphase3[chan1,chan2])
                    phases_list.append(connectphase3[chan1,chan2])
            topcluster_meanphase3 /= float(num_topclusterchannels3)
            self.topcluster_meanphase3_singlesession_list[ana_index].append(topcluster_meanphase3)
            topcluster_stdphase3 = circstd(phases_list)
            self.topcluster_stdphase3_singlesession_list[ana_index].append(topcluster_stdphase3)
            del phases_list
            
            #print("file:%s  topcluster 1: meanphase=%.3f(%.3f) num_channels=%d"%(file,\
            #        topcluster_meanphase1,topcluster_stdphase1,num_topclusterchannels1\
            #    ))
            #print("file:%s  topcluster 2: meanphase=%.3f(%.3f) num_channels=%d"%(file,\
            #        topcluster_meanphase2,topcluster_stdphase2,num_topclusterchannels2\
            #    ))
            #print("file:%s  topcluster 3: meanphase=%.3f(%.3f) num_channels=%d"%(file,\
            #        topcluster_meanphase3,topcluster_stdphase3,num_topclusterchannels3\
            #    ))
            
            ##### check on symmetry in coinnectivity map
            symmetry1 = self.symmetry_check(connect1C,flag_plots)
            symmetry2 = self.symmetry_check(connect2C,flag_plots)
            symmetry3 = self.symmetry_check(connect3C,flag_plots)
            #print("ratio1=%f  ratio2=%f  ratio3=%f"%(symmetry1,symmetry2,symmetry3))
            #print("before:",np.size(self.symmetry_list[0][ana_index]))
            ar = [symmetry1,symmetry2,symmetry3]
            if fmin == 0.3:
                #print("@@@@@@@@@@@@@@@@@@@@@@@@@@ fmin=0.3@@@@@@@@@@@@@@@@@@")
                self.symmetry_list[0][ana_index].append(ar)
                #print("len(ar):%d"%len(ar),len(self.symmetry_list[0][ana_index]))  
                #print("ar:",self.symmetry_list[0][ana_index])
            if fmin == 30.0:
                #print("@@@@@@@@@@@@@@@@@@@@@@@@@@ fmin=30.0@@@@@@@@@@@@@@@@@@")
                self.symmetry_list[1][ana_index].append(ar)
                #print(np.size(self.symmetry_list[1][ana_index]))    
            connections_cluster_all.append(connections_cluster1)
            connections_cluster_all.append(connections_cluster2)
            connections_cluster_all.append(connections_cluster3)
            
            ###### globalize the PLVs with upper clusters removing the base cluster
#           for k in range(len(C1)):
#               chan1 = channels1_list[k][0]
#               chan2 = channels1_list[k][1]
#               if connect1C_sorted_nobase[chan1,chan2]==1:
#                   C1_nobase.append(C1[k])
#           for k in range(len(C2)):
#               chan1 = channels2_list[k][0]
#               chan2 = channels2_list[k][1]
#               if connect2C_sorted_nobase[chan1,chan2]==1:
#                   C2_nobase.append(C2[k])
#           for k in range(len(C3)):
#               chan1 = channels3_list[k][0]
#               chan2 = channels3_list[k][1]
#               if connect3C_sorted_nobase[chan1,chan2]==1:
#                   C3_nobase.append(C3[k])
            self.MWTest_conds_list[ana_index].append([C1_nobase,C2_nobase,C3_nobase])
            #self.MWTest_conds_list[ana_index].append([C1_medium,C2_medium,C3_medium])
            #self.MWTest_conds_list[ana_index].append([C1_top,C2_top,C3_top])
            #self.MWTest_conds_list[ana_index].append([C1,C2,C3])
            
            
            
            ############## determine all channels that participate in a cluster
            cluster1_channels_allclusters_list = [] # contains network nodes of each cluster
            cluster2_channels_allclusters_list = [] # contains network nodes of each cluster
            cluster3_channels_allclusters_list = [] # contains network nodes of each cluster
            for c in range(num_cluster1):
                cluster1_channels_list = []
                for kk in range(len(C1)):
                    if cluster1_channels[c,kk]==1:
                        cluster1_channels_list.append(channels1_list[kk][0])
                        cluster1_channels_list.append(channels1_list[kk][1])
                list_duplrmv1 = [*set(cluster1_channels_list)]
                #print(cluster_channels_list)
                #print(list_duplrmv)
                cluster1_channels_allclusters_list.append(list_duplrmv1)
            for c in range(num_cluster2):
                cluster2_channels_list = []
                for kk in range(len(C2)):
                    if cluster2_channels[c,kk]==1:
                        cluster2_channels_list.append(channels2_list[kk][0])
                        cluster2_channels_list.append(channels2_list[kk][1])
                list_duplrmv2 = [*set(cluster2_channels_list)]
                #print(cluster_channels_list)
                #print(list_duplrmv)
                cluster2_channels_allclusters_list.append(list_duplrmv2)
            for c in range(num_cluster3):
                cluster3_channels_list = []
                for kk in range(len(C3)):
                    if cluster3_channels[c,kk]==1:
                        cluster3_channels_list.append(channels3_list[kk][0])
                        cluster3_channels_list.append(channels3_list[kk][1])
                list_duplrmv3 = [*set(cluster3_channels_list)]
                #print(cluster_channels_list)
                #print(list_duplrmv)
                cluster3_channels_allclusters_list.append(list_duplrmv3)
                
            #print(cluster_mean1sync_sorted)
            #print(cluster_mean2sync_sorted)
            #print(cluster_mean3sync_sorted)
            
            ##### comparison of networks in different conditions ##########
            connect_mutual = np.eye(self.num_chan,dtype=int)
            connect12_mutual = np.eye(self.num_chan,dtype=int)*3
            connect13_mutual = np.eye(self.num_chan,dtype=int)*3
            connect23_mutual = np.eye(self.num_chan,dtype=int)*3
            
            ##### coding of condition comparison:
            ##### if 1 is set and 2 is set     : 3
            ##### if 1 is set and 2 is not set : 1
            ##### if 2 is set and 1 is not set : 2
            ##### 
            ##### if 1 is set and 3 is set     : 3
            ##### if 1 is set and 3 is not set : 1
            ##### if 3 is set and 1 is not set : 2
            ##### 
            ##### if 2 is set and 3 is set     : 3
            ##### if 2 is set and 3 is not set : 1
            ##### if 3 is set and 2 is not set : 2
            ##### 
            ##### 
            ##### 
            for kk in range(len(C1)):
                chan1 = channels1_list[kk][0]
                chan2 = channels1_list[kk][1]
                if connect1C_sorted[chan1,chan2] == num_cluster1: ## 1 is set
                    
                    for kkk in range(len(C2)):
                        if chan1 == channels2_list[kkk][0] and chan2 == channels2_list[kkk][1]:
                            if connect2C_sorted[chan1,chan2] == num_cluster2: ## 2 is set
                                
                                #print("1-2: chan1=%d chan2=%d "%(chan1,chan2))
                                connect12_mutual[chan1,chan2]=3
                                connect12_mutual[chan2,chan1]=3
                                #print("%d %d   1-2"%(chan1,chan2))
                                for kkkk in range(len(C3)):
                                    if chan1 == channels3_list[kkkk][0] and chan2 == channels3_list[kkkk][1]:
                                        #print("1-2-3: chan1=%d chan2=%d    %d (%d)"%(chan1,chan2,connect3C_sorted[chan1,chan2],num_cluster3))
                                        if connect3C_sorted[chan1,chan2] == num_cluster3: # 3 is set 
                                            connect_mutual[chan1,chan2]=1
                                            connect_mutual[chan2,chan1]=1
                                            #print("chan1=%d chan2=%d mutual=%d"%(chan1,chan2,connect_mutual[chan1,chan2]))
                            else: ## 2 is not set
                                connect12_mutual[chan1,chan2]=1
                                connect12_mutual[chan2,chan1]=1
                                
                    for kkk in range(len(C3)):
                        if chan1 == channels3_list[kkk][0] and chan2 == channels3_list[kkk][1]:
                            if connect3C_sorted[chan1,chan2] == num_cluster3: ## 3 is set
                                connect13_mutual[chan1,chan2]=3
                                connect13_mutual[chan2,chan1]=3
                            else:
                                connect13_mutual[chan1,chan2]=1
                                connect13_mutual[chan2,chan1]=1
                
            for kk in range(len(C2)):
                chan1 = channels2_list[kk][0]
                chan2 = channels2_list[kk][1]
                if connect2C_sorted[chan1,chan2] == num_cluster2 : ## 2 is set
                    for kkk in range(len(C1)): ## test on 1
                        if chan1 == channels1_list[kkk][0] and chan2 == channels1_list[kkk][1]:
                            if connect1C_sorted[chan1,chan2] != num_cluster1: ## 1 is not set
                                connect12_mutual[chan1,chan2]=2
                                connect12_mutual[chan2,chan1]=2
                                #print("%d %d   1-2"%(chan1,chan2))
                    
                    for kkk in range(len(C3)): ## test on 3
                        if chan1 == channels3_list[kkk][0] and chan2 == channels3_list[kkk][1]:
                            if connect3C_sorted[chan1,chan2] == num_cluster3: # 3 is set
                                connect23_mutual[chan1,chan2]=3
                                connect23_mutual[chan2,chan1]=3
                            else:
                                connect23_mutual[chan1,chan2]=1
                                connect23_mutual[chan2,chan1]=1
            
            for kk in range(len(C3)):
                chan1 = channels3_list[kk][0]
                chan2 = channels3_list[kk][1]
                if connect3C_sorted[chan1,chan2] == num_cluster3 : ## 3 is set
                    for kkk in range(len(C1)): ## test on 1
                        if chan1 == channels1_list[kkk][0] and chan2 == channels1_list[kkk][1]:
                            if connect1C_sorted[chan1,chan2] != num_cluster1: ## 1 is not set
                                connect13_mutual[chan1,chan2]=2
                                connect13_mutual[chan2,chan1]=2
                                #print("%d %d   1-2"%(chan1,chan2))
                            
                    for kkk in range(len(C2)): ## test on 2
                        if chan1 == channels2_list[kkk][0] and chan2 == channels2_list[kkk][1]:
                            if connect2C_sorted[chan1,chan2] != num_cluster2: # 2 is not set
                                connect23_mutual[chan1,chan2]=2
                                connect23_mutual[chan2,chan1]=2
            
            label = "mutual_channels_%.1f-%.1f_%s_%s.dat"%(fmin,fmax,ana,file)
            f = open(label,"w")
            chanlist = []
            for k in range(self.num_chan):
                for l in range(self.num_chan):
                    if connect_mutual[k,l]==1 and l>k:
                        chanlist.append([k,l])
                        str = "%d %d\n"%(k,l)
                        f.write(str)
            f.close()
            print("The number of channel pairs under investigation is %d"%len(chanlist))
            
            connect_notmutual123 = np.eye(self.num_chan,dtype=int)*0
            mutual_1_list = []
            mutual_1_not_list = []
            connect_notmutual1 = np.eye(self.num_chan,dtype=int)*0
            for kk in range(len(C1)):
                chan1 = channels1_list[kk][0]
                chan2 = channels1_list[kk][1]
                ######## consider those channel pairs, which are top clusters in all conditions
                if connect_mutual[chan1,chan2]==1 and chan2>chan1:
                    mutual_1_list.append(C1[kk])
                ######## consider those channel pairs, which are not top clusters in all conditions
                ######## but which are no in the base cluster
                if connect_mutual[chan1,chan2]==0 and chan2>chan1 and connect1C_sorted[chan1,chan2]>1:
                    mutual_1_not_list.append(C1[kk])
                    connect_notmutual1[chan1,chan2]=1
            mutual_2_list = []
            mutual_2_not_list = []
            connect_notmutual2 = np.eye(self.num_chan,dtype=int)*0
            for kk in range(len(C2)):
                chan1 = channels2_list[kk][0]
                chan2 = channels2_list[kk][1]
                if connect_mutual[chan1,chan2]==1 and chan2>chan1:
                    mutual_2_list.append(C2[kk])
                if connect_mutual[chan1,chan2]==0 and chan2>chan1 and connect2C_sorted[chan1,chan2]>1:
                    mutual_2_not_list.append(C2[kk])
                    connect_notmutual2[chan1,chan2]=1
            mutual_3_list = []
            mutual_3_not_list = []
            connect_notmutual3 = np.eye(self.num_chan,dtype=int)*0
            for kk in range(len(C3)):
                chan1 = channels3_list[kk][0]
                chan2 = channels3_list[kk][1]
                if connect_mutual[chan1,chan2]==1 and chan2>chan1:
                    mutual_3_list.append(C3[kk])
                if connect_mutual[chan1,chan2]==0 and chan2>chan1 and connect3C_sorted[chan1,chan2]>1:
                    mutual_3_not_list.append(C3[kk])
                    connect_notmutual3[chan1,chan2]=1
            
            for kk in range(len(C1)):
                chan1 = channels1_list[kk][0]
                chan2 = channels1_list[kk][1]
                if connect_notmutual1[chan1,chan2]==1 and \
                    connect_notmutual2[chan1,chan2]==1 and \
                    connect_notmutual3[chan1,chan2]==1 :
                        connect_notmutual123[chan1,chan2]=1
                
            self.mutual_connections_data_list[ana_index].append([mutual_1_list,mutual_2_list,mutual_3_list])
            self.mutual_connections_not_data_list[ana_index].append([mutual_1_not_list,mutual_2_not_list,mutual_3_not_list])
            self.connect_notmutual123_list[ana_index].append(connect_notmutual123)
            
            
            ##### coding of condition comparison in connection??_mutual:
            ##### connection12_mutual:
            ##### if 1 is set and 2 is set     : 3
            ##### if 1 is set and 2 is not set : 1
            ##### if 2 is set and 1 is not set : 2
            ##### 
            ##### connection13_mutual:
            ##### if 1 is set and 3 is set     : 3
            ##### if 1 is set and 3 is not set : 1
            ##### if 3 is set and 1 is not set : 2
            ##### 
            ##### connection23_mutual:
            ##### if 2 is set and 3 is set     : 3
            ##### if 2 is set and 3 is not set : 1
            ##### if 3 is set and 2 is not set : 2
            
            ## connect_mutual_num : [cond1-cond2,transition]
            ## cond1-cond2=0 : 1-2
            ## cond1-cond2=1 : 1-3
            ## cond1-cond2=2 : 2-3
            ## transition=0  : 1  1st is set, next is not set
            ## transition=1  : 2  1st is not set, next is set
            ## transition=2  : 3  1st and next are set
            connect_mutual_num = np.zeros((3,3),dtype=int)
            for k in range(self.num_chan):
                for l in range(self.num_chan):
                    if connect12_mutual[k,l]>0 and l>k:
                        connect_mutual_num[0,connect12_mutual[k,l]-1] += 1
                    if connect13_mutual[k,l]>0 and l>k:
                        connect_mutual_num[1,connect13_mutual[k,l]-1] += 1
                    if connect23_mutual[k,l]>0 and l>k:
                        connect_mutual_num[2,connect23_mutual[k,l]-1] += 1
            self.connection_mutual_transitions_list[ana_index].append(connect_mutual_num)
            
#           if flag_inference == 1:
#               ########## compute ANOVA in mutual channels
#               
#               ######### anova experiments have shown that p=0
#               #res_anova = self.compute_topcluster_anova(path,ana,stype_list,file,fmin,fmax,chanlist)
#               ######### all p-values are p=0 what does not yield anywhere 
#               #res_pvalues = self.compute_phase_inference(path,ana,stype_list,file,flag_plots,fmin,fmax,chanlist)
#       
#               connect_mutual_sign = np.zeros((self.num_chan,self.num_chan))
#               for k in range(self.num_chan):
#                   for l in range(self.num_chan):
#                       if connect_mutual[k,l]>0:
#                           if res_pvalues[0][k,l]==0 and res_pvalues[1][k,l]==0 and res_pvalues[2][k,l]==0:
#                               connect_mutual_sign[k,l]=1
#                               connect_mutual_sign[l,k]=1
#               
#               fig0=plt.figure("specific mutual")
#               ax = fig0.add_subplot(211)
#               plt.title("1-2-3 insignificant")
#               plt.imshow(connect_mutual_sign, interpolation='nearest', cmap='Reds',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
#               plt.colorbar()
#               forceAspect(ax,aspect=1)
#               ax.xaxis.set_ticks([])
#               ax.yaxis.set_ticks([])
#               plt.show()
#               
#               if ana == 'low':
#                   self.connection_mutual_sig_list[0].append([res_anova,res_pvalues])
#               if ana == 'medium':
#                   self.connection_mutual_sig_list[1].append([res_anova,res_pvalues])
#               if ana == 'high':
#                   self.connection_mutual_sig_list[2].append([res_anova,res_pvalues])
                    
            connect_mutual_zerophase = np.zeros((self.num_chan,self.num_chan))
            if fmin==0.1 or fmin==0.3:
                thr = 2*np.pi*(1.0/fs)/1000.0
            if fmin==0.5:
                thr = 2*np.pi*(1.0/fs)/500.0
            if fmin==30.0:
                thr = 2*np.pi*(1.0/fs)/25.0
            for k in range(self.num_chan):
                for l in range(self.num_chan):
                    if connect_mutual[k,l]>0:
                        if np.abs(connectphase1[k,l])<thr and \
                            np.abs(connectphase2[k,l])<thr and \
                            np.abs(connectphase3[k,l])<thr:
                            connect_mutual_zerophase[k,l]=1
                            connect_mutual[k,l]=2
            
            ########## If we are looking for strictly in-phase couples,
            ########## then one should consider connect_mutual_zerophase.
            ########## If we relax this condition and allow phase fluctuations
            ########## beyond zero-phase, then we consider connect_mutual.
            ########## Preference: not requesting zero-phase relation since
            ########## it is natural that distant electrodes are not strictly 
            ########## in-phase, unless this is a neural function.
            if ana == 'low':
                self.connection_mutual_list[0].append(connect_mutual)
                #self.connection_mutual_list[0].append(connect_mutual_zerophase)
            if ana == 'medium':
                self.connection_mutual_list[1].append(connect_mutual)
                #self.connection_mutual_list[1].append(connect_mutual_zerophase)
            if ana == 'high':
                self.connection_mutual_list[2].append(connect_mutual)
                #self.connection_mutual_list[2].append(connect_mutual_zerophase)
                    
                    
            if flag_plots == 1:
                
                #fig0 = plt.figure("zero phase mutuals")
                #ax = fig0.add_subplot(111)
                #plt.title("1-2-3 zero phase")
                #plt.imshow(connect_mutual_zerophase, interpolation='nearest', cmap='Reds',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                #plt.colorbar()
                #forceAspect(ax,aspect=1)
                #ax.xaxis.set_ticks([])
                #ax.yaxis.set_ticks([])
                
                label = file+'_1'
                fig0=plt.figure(label)
                ax = fig0.add_subplot(221)
                plt.title("1-2-3")
                plt.imshow(connect_mutual, interpolation='nearest', cmap='Reds',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax = fig0.add_subplot(222)
                plt.title("1-2")
                plt.imshow(connect12_mutual, interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax = fig0.add_subplot(223)
                plt.title("1-3")
                plt.imshow(connect13_mutual, interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax = fig0.add_subplot(224)
                plt.title("2-3")
                plt.imshow(connect23_mutual, interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                fig0.subplots_adjust(hspace=0.3,wspace=0.4)
                
                
                
                if fmin==0.3:
                    label = file+'_'+ana+'_slow'+'_2'
                if fmin==30.0:
                    label = file+'_'+ana+'_gamma'+'_2'
                if fmin!=0.3 and fmin!=30.0:
                    label = file+'_'+ana+'__'+'_2'
                fign=plt.figure(label)
                ax = fign.add_subplot(331)
                plt.imshow(connect1, interpolation='nearest', cmap=plt.cm.jet,origin='upper',vmin=0.0,vmax=1.0)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                #            ax = fign.add_subplot(332)
                #            plt.imshow(connect1C, interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                #            plt.colorbar()
                #            forceAspect(ax,aspect=1)
                ax = fign.add_subplot(332)
                plt.imshow(connect1C_sorted, interpolation='nearest', cmap=plt.cm.jet,origin='upper',vmin=0,vmax=num_cluster1)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax = fign.add_subplot(333)
                minphase=np.min(np.min(connectphase1))
                maxphase=np.max(np.max(connectphase1))
                minmaxphase=max(np.abs(minphase),np.abs(maxphase))
                plt.imshow(connectphase1, interpolation='nearest', cmap='seismic',origin='upper',vmin=-0.2,#minmaxphase,
                    vmax=0.2#minmaxphase
                )#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                #print("cond1: min=%f max=%f"%(minphase,maxphase))
                
                ax = fign.add_subplot(334)
                plt.imshow(connect2, interpolation='nearest', cmap=plt.cm.jet,origin='upper',vmin=0.0,vmax=1.0)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                #ax = fign.add_subplot(324)
                #plt.imshow(connect2C, interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                #plt.colorbar()
                #forceAspect(ax,aspect=1)
                ax = fign.add_subplot(335)
                plt.imshow(connect2C_sorted, interpolation='nearest', cmap=plt.cm.jet,origin='upper',vmin=0,vmax=num_cluster2)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                
                ax = fign.add_subplot(336)
                minphase=np.min(np.min(connectphase2))
                maxphase=np.max(np.max(connectphase2))
                minmaxphase=max(np.abs(minphase),np.abs(maxphase))
                plt.imshow(connectphase2, interpolation='nearest', cmap='seismic',origin='upper',\
                    vmin=-0.2,#minmaxphase
                    vmax=0.2#minmaxphase
                )#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                #print("cond2: min=%f max=%f"%(minphase,maxphase))
                
                ax = fign.add_subplot(337)
                plt.imshow(connect3, interpolation='nearest', cmap=plt.cm.jet,origin='upper',vmin=0.0,vmax=1.0)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                #ax = fign.add_subplot(338)
                #plt.imshow(connect3C, interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                #plt.colorbar()
                #forceAspect(ax,aspect=1)
                
                ax = fign.add_subplot(338)
                cmap = mpl.cm.jet
                #cmap = plt.cm.jet  # define the colormap
                # extract all colors from the .jet map
                #cmaplist = [cmap(i) for i in range(cmap.N)]
                # force the first color entry to be grey
                #cmaplist[0] = (.5, .5, .5, 1.0)
                # create the new map
                #cmap = mpl.colors.LinearSegmentedColormap.from_list(
                #    'Custom cmap', cmaplist, cmap.N)
                # define the bins and normalize
                #bounds = np.linspace(0, 3, 4)
                #norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                
                pylt.imshow(connect3C_sorted, interpolation='nearest', cmap=cmap,origin='upper',vmin=0,vmax=num_cluster3)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                #plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ticks=bounds)
                #ax2 = fig.add_axes([0.95, 0.1, 0.03, 0.8])
                #matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
                #    spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
                
                
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                
                ax = fign.add_subplot(339)
                minphase=np.min(np.min(connectphase3))
                maxphase=np.max(np.max(connectphase3))
                minmaxphase=max(np.abs(minphase),np.abs(maxphase))
                plt.imshow(connectphase3, interpolation='nearest', cmap='seismic',origin='upper',\
                    vmin=-0.2,#minmaxphase,
                    vmax=0.2#minmaxphase
                )#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
                plt.colorbar()
                forceAspect(ax,aspect=1)
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                #print("cond3: min=%f max=%f"%(minphase,maxphase))
                fign.subplots_adjust(hspace=0.3,wspace=0.4)
    
                plt.show()
        
    def compute_phase_inference(self,path,ana,stype_list,file,flag_plots,fmin,fmax,fs,chanlist):
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
        timewindow = 600#600#400#120 # in seconds
        self.FLAG_THR = 0
        label_list = []
        final_results = []
        data_list = []
        list_all = []
        timewindow_list = []
        
        for stype in stype_list:
            print(stype)
            if fs==200.0:
                label = path+file+'_'+ana+'_'+stype+'_notch.h5'
            if fs==500.0:
                label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
            print("datafile:",label)
            label_list.append(stype)
            data=self.read_data_h5(label)
            data_list.append(data)
            print("shape of data:",np.shape(data))
            if self.num_time<(timewindow/self.dt):
                print("self.num_time<(timewindow/self.dt): %d < %d    : choose duration %d"%(self.num_time,int(timewindow/self.dt),self.num_time))
                timewindow0=self.num_time*self.dt
            else:
                timewindow0 = timewindow
            timewindow_list.append(timewindow0)
        
#       print("chanlist:",chanlist)
#       chanlist0 = []
#       counter = 0
#       for l in chanlist:
#           chanlist0.append(l)
#           counter +=1
#           if counter>3:
#               break 
#       del chanlist
#       chanlist = chanlist0
#       print("chanlist:",chanlist)
        
        [chan11,results11] = self.Compute_phaseinference_band_intra(\
                                data_list[0],fmin,fmax,timewindow_list[0],chanlist)
        #[chan12,results12] = self.Compute_phaseinference_band(\
        #                        data_list[0],data_list[1],\
        #                        fmin,fmax,\
        #                        timewindow_list[0],timewindow_list[1],\
        #                        chanlist)
        #[chan13,results13] = self.Compute_phaseinference_band(\
        #                        data_list[0],data_list[2],\
        #                        fmin,fmax,\
        #                        timewindow_list[0],timewindow_list[2],\
        #                        chanlist)
        #[chan23,results23] = self.Compute_phaseinference_band(\
        #                        data_list[1],data_list[2],\
        #                        fmin,fmax,\
        #                        timewindow_list[1],timewindow_list[2],\
        #                        chanlist)
        
        inference_list = []
        
        ## question: should one apply a Bonferroni correction ?
        alpha = 0.05/len(chanlist)
        
        label = "pvalues11_%.1f-%.1f_%s_%s.dat"%(fmin,fmax,ana,file)
        f = open(label,"w")
        inference11 = np.eye(self.num_chan,dtype=int)*0
        #print(np.shape(chan11))
        #print(chan11)
        #print(type(chan11))
        for k in range(len(chan11)):
            #print("k=%d "%k)
            #print("chan11[k]:",chan11[k])
            chanpair = chan11[k] 
            chan1 = chanlist[chanpair[0]][0]
            chan2 = chanlist[chanpair[0]][1]
            chan_11 = chanlist[chanpair[1]][0]
            chan_22 = chanlist[chanpair[1]][1]
            if results11[k]<alpha:
                str = "%d %d     %d %d  %g\n"%(chan1,chan2,chan_11,chan_22,results11[k])
                #print("p-values:",str)
                f.write(str)
        f.close()
        
        
#       label = "pvalues12_%.1f-%.1f_%s_%s.dat"%(fmin,fmax,ana,file)
#       f = open(label,"w")
#       inference12 = np.eye(self.num_chan,dtype=int)*0
#       for k in range(len(chan12)):
#           chan1 = chan12[k][0]
#           chan2 = chan12[k][1]
#           if results12[k]<alpha:
#               inference12[chan1,chan2]=1
#               str = "%d %d  %g\n"%(chan1,chan2,results12[k])
#               f.write(str)
#       f.close()
#       inference_list.append(inference12)
#       
#       label = "pvalues13_%.1f-%.1f_%s_%s.dat"%(fmin,fmax,ana,file)
#       f = open(label,"w")
#       inference13 = np.eye(self.num_chan,dtype=int)*0
#       for k in range(len(chan13)):
#           chan1 = chan12[k][0]
#           chan2 = chan12[k][1]
#           if results13[k]<alpha:
#               inference13[chan1,chan2]=1
#               str = "%d %d  %g\n"%(chan1,chan2,results13[k])
#               f.write(str)
#       f.close()
#       inference_list.append(inference13)
#       
#       label = "pvalues23_%.1f-%.1f_%s_%s.dat"%(fmin,fmax,ana,file)
#       f = open(label,"w")
#       inference23 = np.eye(self.num_chan,dtype=int)*0
#       for k in range(len(chan23)):
#           chan1 = chan23[k][0]
#           chan2 = chan23[k][1]
#           if results23[k]<alpha:
#               inference23[chan1,chan2]=1
#               str = "%d %d  %g\n"%(chan1,chan2,results23[k])
#               f.write(str)
#       f.close()
#       inference_list.append(inference23)
        
#       if flag_plots == 1:
#           fig = plt.figure(5)
#           ax = plt.subplot(3,1,1)
#           plt.imshow(inference12, interpolation='nearest', cmap="Reds",origin='upper',vmin=0.0,vmax=1.0)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
#           plt.colorbar()
#           self.forceAspect(ax,aspect=3)
#           ax = plt.subplot(3,1,2)
#           plt.imshow(inference13, interpolation='nearest', cmap="Reds",origin='upper',vmin=0.0,vmax=1.0)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
#           plt.colorbar()
#           self.forceAspect(ax,aspect=3)
#           ax = plt.subplot(3,1,3)
#           plt.imshow(inference23, interpolation='nearest', cmap="Reds",origin='upper',vmin=0.0,vmax=1.0)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
#           plt.colorbar()
#           self.forceAspect(ax,aspect=3)
#           plt.show()
        
        return inference_list
    
    def compute_topcluster_anova(self,path,ana,stype_list,file,fmin,fmax,fs,chanlist):
        
        timewindow = 600#600#400#120 # in seconds
        label_list = []
        final_results = []
        
        anova_list = []
        for stype in stype_list:
            if fs==200.0:
                label = path+file+'_'+ana+'_'+stype+'_notch.h5'
            if fs==500.0:
                label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
            label_list.append(stype)
            data=self.read_data_h5(label)
            if self.num_time<(timewindow/self.dt):
                print("self.num_time<(timewindow/self.dt): %d < %d    : choose duration %d"%(self.num_time,int(timewindow/self.dt),self.num_time))
                timewindow0=self.num_time*self.dt
            else:
                timewindow0 = timewindow
            
            phase_list= self.compute_phases(data,fmin,fmax,timewindow0,stype,chanlist)
            final_results.append(phase_list)
        
            ####### write out data for R script
            num_chanpairs = len(phase_list)
            f = open("datain_anova.dat","w")
            for k  in range(num_chanpairs):
                group = k+1
                ph = phase_list[k]
                num_times = np.size(ph)
                for i in range(num_times):
                    str = '%f  %d\n'%(ph[i],group)
                    #print(str)
                    f.write(str)
            f.close()
            import subprocess
            res = subprocess.call("Rscript ./PhaseAnova.R", shell=True)
            f = open("out.dat","r")
            r=float(f.read())
            anova_list.append(r)
            print("this was type : %s with p-value:%f"%(stype,r))
            res = subprocess.call("rm ./datain_anova.dat", shell=True)
        return anova_list
                
                
            
    def compute_connectivity(self,path,ana,stype_list,file_list,fmin,fmax,fs):
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
        ####### delta band
        #fmin = 0.5
        #fmax = 4.5#5.0
        
        ####### gamma band
        #fmin = 30.0
        #fmax = 60.0
        timewindow = 600#600#400#120 # in seconds
        self.FLAG_THR = 0
        label_list = []
        final_results = []
        power_delta_list = []
        #label = self.filelabel_connectivity+'_'+file+'_'+ana+'.dat'
        #print(label)
        #quit()

        #label0 = '_%.3f-%.3f_'%(fmin,fmax)
        #label=self.filelabel_connectivity+label0+file+'_'+ana+'.dat'
        #print("label:",label)
        #quit()

        for file in file_list:
            label_list = []
            final_results = []
            power_delta_list = []
            list_all = []
            for stype in stype_list:
                print(stype)
                if fs==200.0:
                    label = path+file+'_'+ana+'_'+stype+'_notch.h5'
                if fs==500.0:
                    label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
                print("datafile:",label)
                label_list.append(stype)
                data=self.read_data_h5(label)
                #print("shape of data:",np.shape(data))
                #data_filtered = self.Compute_filter_butterworth(data,fmin,fmax)
                if self.num_time<(timewindow/self.dt):
                    print("self.num_time<(timewindow/self.dt): %d < %d    : choose duration %d"%(self.num_time,int(timewindow/self.dt),self.num_time))
                    timewindow0=self.num_time*self.dt
                else:
                    timewindow0 = timewindow
                results = self.Compute_phasecoherence_band(data,fmin,fmax,timewindow0,stype)
                final_results.append(results)
            
            label0 = '_%.1f-%.1f_'%(fmin,fmax)
            #label=self.filelabel_connectivity+'_n_'+label0+file+'_'+ana+'.dat'
            if fs==200.0:
                label=self.filelabel_connectivity+'_n_'+label0+file+'_'+ana+'.dat'
            if fs==500.0:
                label=self.filelabel_connectivity+'_n_'+label0+file+'_'+ana+'_fs500.0'+'.dat'
            
            f = open(label,"w")
            counter = 0
            #connectivity_diff = np.zeros((3,self.num_chan,self.num_chan))
            for chan1 in range(self.num_chan-1):
                chan01 = chan1+1
                for chan02 in range(self.num_chan-chan01):
                    chan2 = chan01+chan02
                    
                    cond_list = final_results[0] 
                    m1      = cond_list[2][counter][0]
                    n1      = cond_list[3]
                    PM1     = cond_list[2][counter][1]
                    PM_sd1  = cond_list[2][counter][2]
                    del cond_list
                    cond_list = final_results[1] 
                    m2      = cond_list[2][counter][0]
                    n2      = cond_list[3]
                    PM2     = cond_list[2][counter][1]
                    PM_sd2  = cond_list[2][counter][2]
                    del cond_list
                    cond_list = final_results[2] 
                    m3      = cond_list[2][counter][0]
                    n3      = cond_list[3]
                    PM3     = cond_list[2][counter][1]
                    PM_sd3  = cond_list[2][counter][2]
                    del cond_list
                    
                    print("%d %d "%(chan1,chan2),\
                        "   m1:%.3f"%m1," n1:",n1,\
                        "   m2:%.3f"%m2," n2:",n2,\
                        "   m3:%.3f"%m3," n3:",n3,\
                        "   PM1:%.3f"%PM1," Pstdv1:%.3f"%PM_sd1,\
                        "   PM2:%.3f"%PM2," Pstdv2:%.3f"%PM_sd2,\
                        "   PM3:%.3f"%PM3," Pstdv3:%.3f"%PM_sd3)
                    str = '%d %d   %f %d  %f %d   %f %d   %f %f  %f %f  %f %f\n'%(chan1,chan2,m1,n1,m2,n2,m3,n3,\
                            PM1,PM_sd1,PM2,PM_sd2,PM3,PM_sd3)
                    f.write(str)
                    counter += 1
    #              s,p=ss.ttest_ind_from_stats(m1,sdv1,n1,m2,sdv2,n2,equal_var=False)
    #              if p<0.05:
    #                   connectivity_diff[0,chans[0],chans[1]]=m2-m1
    #                   connectivity_diff[0,chans[1],chans[0]]=m2-m1
    #              s,p=ss.ttest_ind_from_stats(m1,sdv1,n1,m3,sdv3,n3,equal_var=False)
    #              if p<0.05:
    #                   connectivity_diff[1,chans[0],chans[1]]=m2-m1
    #                   connectivity_diff[1,chans[1],chans[0]]=m2-m1
    #              s,p=ss.ttest_ind_from_stats(m2,sdv2,n2,m3,sdv3,n3,equal_var=False)
    #              if p<0.05:
    #                   connectivity_diff[2,chans[0],chans[1]]=m2-m1
    #                   connectivity_diff[2,chans[1],chans[0]]=m2-m1
    #       
    #       self.clustering()
    #       
            f.close()
            
            fign=plt.figure(1)
            ax = fign.add_subplot(321)
            plt.imshow(final_results[0][0], interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
            plt.colorbar()
            forceAspect(ax,aspect=2)
            
            ax = fign.add_subplot(322)
            max_ = np.max(np.max(final_results[0][1]))
            min_ = np.min(np.min(final_results[0][1]))
            minmax = max(max_,min_)
            plt.imshow(final_results[0][1], interpolation='nearest', cmap='seismic',origin='upper',vmin=-minmax,vmax=minmax)#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
            plt.colorbar()
            forceAspect(ax,aspect=2)
            
            ax = fign.add_subplot(323)
            plt.imshow(final_results[1][0], interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
            plt.colorbar()
            forceAspect(ax,aspect=2)
            
            ax = fign.add_subplot(324)
            max_ = np.max(np.max(final_results[1][1]))
            min_ = np.min(np.min(final_results[1][1]))
            minmax = max(max_,min_)
            plt.imshow(final_results[1][1], interpolation='nearest', cmap='seismic',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
            plt.colorbar()
            forceAspect(ax,aspect=2)   
            
            ax = fign.add_subplot(325)
            plt.imshow(final_results[2][0], interpolation='nearest', cmap=plt.cm.jet,origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
            plt.colorbar()
            forceAspect(ax,aspect=2)
            
            ax = fign.add_subplot(326)
            max_ = np.max(np.max(final_results[2][1]))
            min_ = np.min(np.min(final_results[2][1]))
            minmax = max(max_,min_)
            plt.imshow(final_results[2][1], interpolation='nearest', cmap='seismic',origin='upper')#, #extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
            plt.colorbar()
            forceAspect(ax,aspect=2)
            
            del label_list
            del final_results
            del power_delta_list
            del list_all
            
            plt.show()
    
    def extract_set_timeseries(self,path_,parameters,file,channum):
        #file_ = ['20110913B']
        #stype_list = ['pre','post']

        fig0=plt.figure("time series")
        len_ = len(parameters)
        counter = 0
        for counter in range(len_):
            ana = parameters[counter][0]
            stype = parameters[counter][1] 
            path = path_+ana+'/'
            signal = self.show_timeseries(path,ana,stype,file,channum)
            orig_signal = signal[0] # 0: stype=pre  0: orig_data
            slow_signal = signal[1] # 0: stype=pre  0: orig_data
            gamma_signal = signal[2] # 0: stype=pre  0: orig_data
            ax=plt.subplot(len_,3,(counter)*3+1)
            plt.plot(range(len(orig_signal)),orig_signal,'k')
            #labeltitle = '%s - %s   all'%(ana,stype)
            #plt.title(labeltitle)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.ylim(-4000,4000)
            plt.yticks(color='w')
            if counter<3:
                plt.xticks(color='w')
            
            ax=plt.subplot(len_,3,(counter)*3+2)
            plt.plot(range(len(orig_signal)),slow_signal,'b')
            #labeltitle = '%s - %s   slow'%(ana,stype)
            #plt.title(labeltitle)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.ylim(-4000,4000)
            plt.yticks(color='w')
            if counter<3:
                plt.xticks(color='w')
            
            ax=plt.subplot(len_,3,(counter)*3+3)
            plt.plot(range(len(orig_signal)),gamma_signal,'r')
            #labeltitle = '%s - %s   gamma'%(ana,stype)
            #plt.title(labeltitle)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.ylim(-500,500)
            plt.yticks(color='w')
            if counter<3:
                plt.xticks(color='w')
        plt.show()
        
    def compute_amplitudes(self,path,ana,stype,file):
        #f0min = 0.1
        #f0max = 2.5
        f0min = 0.1
        f0max = 80.0
        fmin = 0.3
        fmax = 2.5
        f1min = 30.0
        f1max = 70.0
        self.fs = 500.0
        fs = self.fs
        counter = 0
        
        from scipy.signal import hilbert
                    
        print(stype)
        if fs==200.0:
            label = path+file+'_'+ana+'_'+stype+'_notch.h5'
        if fs==500.0:
            label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
        print("datafile:",label)
        data=self.read_data(label)
        
        signals_list = []
        signals_list.append([])
        signals_list.append([])
        power_list = []
        power_list.append([])
        power_list.append([])
        
        #signals_list.append([])
        fig0 = plt.figure("Phase-amplitude coupling ?")
        for kk in range(self.num_chan):
            k = kk+1
            print("channel %d(%d)"%(kk,self.num_chan))
            #data_filtered = self.Compute_filter_butterworth(data[:,kk],fmin,fmax,4)
            data_filtered1 = self.Compute_filter_butterworth(data[:,kk],f1min,f1max,4)
            #signals_list[0].append(data[:,kk])
            #signals_list[1].append(data_filtered)
            signals_list[0].append(data_filtered1)
            signals_list[1].append(np.abs(hilbert(data_filtered1)))
            res = self.Compute_powerspectrum_band_freqs(signals_list[0],f1min,f1max) # gamma range
            power_list[0].append(res)
            res1 = self.Compute_powerspectrum_band_freqs(signals_list[1],f0min,f0max) # full frequency range
            power_list[1].append(res1)
            plt.subplot(4,4,k)
            plt.plot(res[0],res[1],'b')
            plt.plot(res1[0],res1[1],'r')
        plt.show()
 
    
    
    def show_timeseries(self,path,ana,stype,file,channum):
        
        #f0min = 0.1
        #f0max = 2.5
        f0min = 0.1
        f0max = 80.0
        fmin = 0.3
        fmax = 2.5
        f1min = 30.0
        f1max = 70.0
        self.fs = 500.0
        fs = self.fs
        counter = 0
        
        print(stype)
        if fs==200.0:
            label = path+file+'_'+ana+'_'+stype+'_notch.h5'
        if fs==500.0:
            label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
        print("datafile:",label)
        data=self.read_data(label)
        
        signals_list = []
        signals_list.append([])
        signals_list.append([])
        signals_list.append([])
        for kk in range(self.num_chan):
            print("channel %d(%d)"%(kk,self.num_chan))
            data_filtered = self.Compute_filter_butterworth(data[:,kk],fmin,fmax,4)
            data_filtered1 = self.Compute_filter_butterworth(data[:,kk],f1min,f1max,4)
            signals_list[0].append(data[:,kk])
            signals_list[1].append(data_filtered)
            signals_list[2].append(data_filtered1)
            #res = self.Compute_powerspectrum(data[:,kk])
        
        duration = len(signals_list[0][0])
        print(np.shape(signals_list[0]))
        #figlabel = '%s'%stype
        #fig2 = plt.figure(figlabel,figsize=(10,7))
        #for chan in range(self.num_chan):
        #    chan0 = chan+1
        #    ax = plt.subplot(4,4,chan0)
        #    plt.plot(range(duration),signals_list[0][chan],'r',range(duration),signals_list[1][chan],'k',range(duration),signals_list[2][chan])
    
        list_ = []
        for l in range(3):# loop over orig data, slow data and gamma data
            list_.append(signals_list[l][channum][-5000:])## append time series of the last 1000 points
        #plt.show()
        return list_
            
 
    def powerbands_all(self,path_,ana_list,stype_list,file_list):
        
        import statistics as stat 
        import scipy.stats as ss
        
        num_stype = 3
        counter = 0
        counter_all = 0
        power_list = []
        power_filtered_list = []
        bands_list = []
        label_list = []
        freqs_list = []
        
        #f0min = 0.1
        #f0max = 2.5
        f0min = 0.1
        f0max = 80.0
        
        
        fmin = 0.5
        fmax = 2.5
        f1min = 30.0
        f1max = 70.0
        
        self.fs = 500.0
        fs = self.fs
        
        Dpower_list = []
        Dpower_chans_list = []
        p_list = []
        num_list = len(file_list)
        
        from scipy.stats import wilcoxon
        import statistics as stat
        for k in range(num_list):# loop over animals
            file = file_list[k]
            Dpower_list.append([])
            Dpower_chans_list.append([])
            p_list.append([])
            
            for ana in ana_list:# loop over ana_levels
                path = path_+ana+'/'
                dpower_ = []
                
                stype=stype_list[0]
                print(stype)
                if fs==200.0:
                    label = path+file+'_'+ana+'_'+stype+'_notch.h5'
                if fs==500.0:
                    label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
                print("datafile:",label)
                data=self.read_data(label)
                
                #power_slow_avg1 = 0.0
                #power_gamma_avg1 = 0.0
                
                pop1_list = []
                pop1_list.append([])
                pop1_list.append([])
                for kk in range(self.num_chan):
                    #print("channel %d(%d)"%(kk,self.num_chan))
                    power_slow = self.Compute_powerspectrum_band(data[:,kk],fmin,fmax)
                    power_gamma = self.Compute_powerspectrum_band(data[:,kk],f1min,f1max)
                    #power_slow_avg1 += power_slow/float(self.num_chan)
                    #power_gamma_avg1 += power_gamma/float(self.num_chan)
                    pop1_list[0].append(power_slow)
                    pop1_list[1].append(power_gamma)
                power_slow_avg1 = stat.median(pop1_list[0])
                power_gamma_avg1 = stat.median(pop1_list[1])
                print("avg: slow power=%f   gamma power=%f"%(power_slow_avg1,power_gamma_avg1))
                
                stype=stype_list[1]
                print(stype)
                if fs==200.0:
                    label = path+file+'_'+ana+'_'+stype+'_notch.h5'
                if fs==500.0:
                    label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
                print("datafile:",label)
                data=self.read_data(label)
                #power_slow_avg2 = 0.0
                #power_gamma_avg2 = 0.0
                pop2_list = []
                pop2_list.append([])
                pop2_list.append([])
                for kk in range(self.num_chan):
                    #print("channel %d(%d)"%(kk,self.num_chan))
                    power_slow = self.Compute_powerspectrum_band(data[:,kk],fmin,fmax)
                    power_gamma = self.Compute_powerspectrum_band(data[:,kk],f1min,f1max)
                    #power_slow_avg2 += power_slow/float(self.num_chan)
                    #power_gamma_avg2 += power_gamma/float(self.num_chan)
                    pop2_list[0].append(power_slow)
                    pop2_list[1].append(power_gamma)
                power_slow_avg2 = stat.median(pop2_list[0])
                power_gamma_avg2 = stat.median(pop2_list[1])
                print("avg: slow power=%f   gamma power=%f"%(power_slow_avg2,power_gamma_avg2))
        
                #print("means:  ",np.mean(pop1_list[0]),' ',np.mean(pop2_list[0]))
                res_slow = wilcoxon(pop1_list[0], y=pop2_list[0], zero_method='wilcox', correction=False, alternative='two-sided')
                res_gamma = wilcoxon(pop1_list[1], y=pop2_list[1], zero_method='wilcox', correction=False, alternative='two-sided')
                print(res_slow[1])
                print(res_gamma[1])
                p_list[k].append([res_slow[1],res_gamma[1]])
        
                Dpower_slow = (power_slow_avg2-power_slow_avg1)/power_slow_avg1
                Dpower_gamma = (power_gamma_avg2-power_gamma_avg1)/power_gamma_avg1
                Dpower_list[k].append([Dpower_slow,Dpower_gamma])
        
        
        labels = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7']
        #labels = ['A1', 'A2']
        x = np.arange(len(labels))  # the label locations
        width = 0.25  # the width of the bars
        fig, ax = plt.subplots(2,1)
        #ax[0].axis('off')
        #ax[1].axis('off')
        #ax[2].axis('off')
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)
        #ax[2].spines['top'].set_visible(False)
        #ax[2].spines['right'].set_visible(False)
        
        from matplotlib import ticker
        ####### slow
        h0 = []
        for k in range(num_list):
            p = Dpower_list[k][0][0]
            #if p<0:
            #    p = - np.log(np.abs(p))
            #else:
            #    p = np.log(p)
            h0.append(p)# ana=low
        #print(np.shape(x),np.shape(h0))
        #print(h0)
        
        rects01 = ax[0].bar(x - width, h0, width, align='center', ecolor='lightgrey',label='low')
        #ax[0].boxplot(h0,positions = x, notch=True)
        h1 = []
        for k in range(num_list):
            p = Dpower_list[k][1][0]
            #if p<0:
            #    p = - np.log(np.abs(p))
            #else:
            #    p = np.log(p)
            h1.append(p)#ana=high
        rects02 = ax[0].bar(x        , h1, width, align='center', ecolor='lightgrey', label='high')
        ax[0].axhline(0, color='grey', linewidth=0.8)
        ax[0].xaxis.set_major_locator(ticker.FixedLocator(x))
        ax[0].xaxis.set_major_formatter(ticker.FixedFormatter([]))
        #ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        #ax[0].yaxis.set_major_locator(ticker.FixedLocator([-maxmax,0,maxmax]))
        alpha = 0.05
        maxx0 = max(h0)
        maxx1 = max(h1)
        for k in range(num_list):
            if p_list[k][0][0]<alpha:
                val = h0[k]*1.1
                if val<0:
                    val = 0.1*maxx0
                ax[0].text(x[k]-width, val,'*', ha='center', va='bottom', color='k') 
            if p_list[k][1][0]<alpha:
                val = h1[k]*1.1
                if val<0:
                    val = 0.1*maxx1
                ax[0].text(x[k], val,'*', ha='center', va='bottom', color='k')
        
        ####### gamma
        h0 = []
        for k in range(num_list):
            p = Dpower_list[k][0][1]
            #if p<0:
            #    p = -np.log(np.abs(p))
            #else:
            #    p = np.log(p)
            h0.append(p)
        
        rects01 = ax[1].bar(x - width, h0, width, align='center', ecolor='lightgrey',label='low')
        h1 = []
        for k in range(num_list):
            p = Dpower_list[k][1][1]
            #if p<0:
            #    p = - np.log(np.abs(p))
            #else:
            #    p = np.log(p)
            h1.append(p)
        rects02 = ax[1].bar(x        , h1, width, align='center', ecolor='lightgrey', label='high')
        ax[1].axhline(0, color='grey', linewidth=0.8)
        ax[1].xaxis.set_major_locator(ticker.FixedLocator(x))
        ax[1].xaxis.set_major_formatter(ticker.FixedFormatter([]))
        #ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        #ax[0].yaxis.set_major_locator(ticker.FixedLocator([-maxmax,0,maxmax]))
        maxx0 = max(h0)
        maxx1 = max(h1)
        for k in range(num_list):
            if p_list[k][0][1]<alpha:
                val = h0[k]*1.1
                if val<0:
                    val = 0.1*maxx0
                ax[1].text(x[k]-width, val,'*', ha='center', va='bottom', color='k') 
            if p_list[k][1][1]<alpha:
                val = h1[k]*1.1
                if val<0:
                    val = 0.1*maxx1
                ax[1].text(x[k], val,'*', ha='center', va='bottom', color='k')
                
                
        
        plt.show()
        
        
    def average_powerbands(self,path,ana,stype_list,file):
        
        import statistics as stat 
        import scipy.stats as ss
        
        num_stype = 3
        counter = 0
        counter_all = 0
        power_list = []
        power_filtered_list = []
        bands_list = []
        label_list = []
        freqs_list = []
        
        #f0min = 0.1
        #f0max = 2.5
        f0min = 0.1
        f0max = 80.0
        
        
        fmin = 0.5
        fmax = 2.5
        f1min = 30.0
        f1max = 70.0
        
        self.fs = 500.0
        fs = self.fs
        for stype in stype_list:
            print(stype)
            if fs==200.0:
                label = path+file+'_'+ana+'_'+stype+'_notch.h5'
            if fs==500.0:
                label = path+file+'_'+ana+'_'+stype+'_fs500.0'+'_notch.h5'
            print("datafile:",label)
            label_list.append(stype)
            data=self.read_data(label)
            
            
            power_bands = []
            power_spectra = []
            power_spectra_filtered = []
            signals_list = []
            signals_list.append([])
            signals_list.append([])
            signals_list.append([])
            for kk in range(self.num_chan):
                print("channel %d(%d)"%(kk,self.num_chan))
                data_filtered = self.Compute_filter_butterworth(data[:,kk],fmin,fmax)
                data_filtered1 = self.Compute_filter_butterworth(data[:,kk],f1min,f1max)
                signals_list[0].append(data_filtered)
                signals_list[1].append(data_filtered1)
                res = self.Compute_powerspectrum(data[:,kk])
                res_filtered = self.Compute_powerspectrum(data_filtered)
                res_filtered1 = self.Compute_powerspectrum(data_filtered1)
                if kk == 0:
                    freqs = res[1]
                poweravg = 1.0#np.mean(res[2])
                power_rel = res[2]/poweravg
                power_bands.append(res[0])
                power_spectra.append(power_rel)
                power_spectra_filtered.append([res_filtered[2],res_filtered1[2]])
            print("length of power_spectra_filtered:",len(power_spectra_filtered),np.shape(power_spectra_filtered))
            
            freqs_list.append(freqs)    
            power_list.append(power_spectra)
            power_filtered_list.append(power_spectra_filtered)
            bands_list.append(power_bands)
            
            counter += 1
            counter_all += 1
        
        duration = len(signals_list[0][0])
        fig2 = plt.figure(file,figsize=(10,7))
        for chan in range(self.num_chan):
            chan0 = chan+1
            ax = plt.subplot(4,4,chan0)
            plt.plot(range(duration),signals_list[0][chan])
        plt.show()
        quit()
        
        print(len(power_filtered_list))
        res_filtered_all = power_filtered_list[0]
        print(len(res_filtered))
        for k in range(len(res_filtered)):
            print(np.shape(res_filtered_all[k]))
        #quit()
        
        df = freqs[1]-freqs[0]
        num_f0min = int((f0min-freqs[0])/df)
        num_f0max = int((f0max-freqs[0])/df)
        freqs_n = freqs[num_f0min:num_f0max]
        
        print("num of freqs:",len(freqs_list[0]))
        fig1 = plt.figure(file,figsize=(10,7))
        for chan in range(self.num_chan):
            chan0 = chan+1
            #plt.title(label)
            ax = plt.subplot(4,4,chan0)
            power_n = (power_list[0][chan][:])[num_f0min:num_f0max]
            #plt.plot(freqs_n,power_n,color='k',label='pre')
            plt.plot(freqs_n,np.log(power_n),color='k',label='pre')
            #plt.plot(freqs_list[0],np.log(power_n),color='k',label='pre')
            #plt.plot(np.log(freqs_list[0]),np.log(power_n),color='k',label='pre')
            #plt.plot(freqs_list[0],np.log(power_filtered_list[0][chan][0][:]),'k--')
            #plt.plot(freqs_list[0],np.log(power_filtered_list[0][chan][1][:]),'k--')
            #plt.ylim(0.8,16.0)
            if chan0<=12:
                plt.xticks(color='w')
                #ax.xaxis.set_ticks([])
            plt.xticks(color='w')
            plt.yticks(color='w')
            #ax.yaxis.set_ticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            #plt.plot(np.log(freqs_list[0]),np.log(power_list[0][chan][:]),color='k',label='pre')
            #print("shape:",np.shape(bands_list[0][chan][8]))
            #limit = (bands_list[0][chan][8])
            #low_limit0= np.log(power_list[0][chan][:] - limit)
            #high_limit0 = np.log(power_list[0][chan][:] + limit)
            #print("limit:",limit)
            #plt.fill_between(freqs_list[0], low_limit0, high_limit0, color='black', alpha=0.2)
            #plt.fill_between(np.log(freqs_list[0]), low_limit0, high_limit0, \
            #                    color='black', alpha=0.2)
            
            power_n = (power_list[1][chan][:])[num_f0min:num_f0max]
            #plt.plot(freqs_n,power_n,color='r',label='stim')
            plt.plot(freqs_n,np.log(power_n),color='r',label='stim')
            #plt.plot(freqs_list[1],np.log(power_list[1][chan][:]),color='r',label='stim')
            #plt.plot(np.log(freqs_list[1]),np.log(power_list[1][chan][:]),color='r',label='stim')
            #plt.ylim(0.8,16.0)
            #plt.plot(np.log(freqs_list[1]),np.log(power_list[1][chan][:]),color='r',label='stim')
            #limit = (bands_list[1][chan][8])
            #low_limit1= np.log(power_list[1][chan][:] - limit)
            #high_limit1 = np.log(power_list[1][chan][:] + limit)
            #plt.fill_between(freqs_list[1], low_limit1, high_limit1,color='red', alpha=0.2)
            #plt.fill_between(np.log(freqs_list[1]), low_limit1, high_limit1, \
            #                    color='red', alpha=0.2)
            
            power_n = (power_list[2][chan][:])[num_f0min:num_f0max]
            #plt.plot(freqs_n,power_n,color='b',label='post')
            plt.plot(freqs_n,np.log(power_n),color='b',label='post')
            #plt.plot(freqs_list[2],np.log(power_list[2][chan][:]),color='b',label='post')
            #plt.plot(np.log(freqs_list[2]),np.log(power_list[2][chan][:]),color='b',label='post')
            #plt.plot(np.log(freqs_list[2]),np.log(power_list[2][chan][:]),color='b',label='post')
            #plt.ylim(0.8,16.0)
            #limit = (bands_list[2][chan][8])
            #low_limit2= np.log(power_list[2][chan][:] - limit)
            #high_limit2 = np.log(power_list[2][chan][:] + limit)
            #plt.fill_between(freqs_list[2], low_limit2, high_limit2,color='blue', alpha=0.2)
            #plt.fill_between(np.log(freqs_list[2]), low_limit2, high_limit2, \
            #                    color='blue', alpha=0.2)
            if chan==0:
                ax.legend(loc=1,fontsize=8)
        fig1.subplots_adjust(hspace=0.5,wspace=0.5)
        plt.show()
        label = './spectra_'+file+'.png'
        print("write ",label)
        plt.savefig(label,format='png',bbox_inches='tight')
        #quit()
        #return freqs,power_all
    

    def plot_data(self):
        
        bands_channels = []
        
        fig = plt.figure(1)
        
        power_all = np.zeros(self.num_chan)
        numchannels_chosen = self.num_chan#4
        #for kk in range(self.num_chan):
        for kk in range(numchannels_chosen):
            #print("kk=%d"%kk)
            k = kk+1
            ax = plt.subplot(numchannels_chosen,2,2*k-1)            
            plt.plot(self.times,self.data[kk,:],'k*')
            str = "chn %d"%self.channel_indices[kk]
            plt.ylabel(str)
            if kk==self.num_chan-1:
                plt.xlabel("time [s]")
            if kk==0:
                plt.title("time series")
                
            ax = plt.subplot(numchannels_chosen,2,2*k)            
            res = self.Compute_powerspectrum(self.data[kk,:])
            bands_channels.append(res[0])
            power_f = res[1]
            lmax = len(power_f)
            for l in range(len(power_f)):
                if power_f[l]>50.0:
                    lmax=l-1
                    break
            power_p = res[2]
            power_all[kk]=np.sum(power_p)
            
            plt.plot((power_f[:lmax]),np.log(power_p[:lmax]))
            #plt.plot(np.log(power_f[:lmax]),np.log(power_p[:lmax]))
            #plt.plot(power_f[:lmax],(power_p[:lmax]))
            #plt.plot(self.times_stored,self.data_stored[kk,:],'r')
            if kk==self.num_chan-1:
                plt.xlabel("frequency [Hz]")
            if kk==0:
                plt.title("power spectrum")
                
        fig.subplots_adjust(hspace=0.5)
        
#       fig = plt.figure(0)
#       plt.plot(self.times_stored,self.signal)
#       
#       fig = plt.figure(2)
#       plt.plot(range(self.num_chan),power_all,'*')
#       
#       fig = plt.figure(3)
#       num_bands=4
#       channels_bands = np.zeros(self.num_chan)
#       band_labels = ['delta','theta','alpha','beta']
#       
#       #for kk in range(self.num_chan):
#       for kk in range(num_bands):
#           k = kk+1
#           ax = plt.subplot(num_bands,1,k)            
#           
#           print("shape channels_band:",np.shape(channels_bands),"  bands_channels:",np.shape(bands_channels))
#           for l in range(numchannels_chosen):
#               channels_bands[l] = (bands_channels[l])[2*kk]
#           plt.plot(range(self.num_chan),channels_bands)
#           str = "channels"
#           plt.xlabel(str)
#           str = "%s"%band_labels[kk]
#           plt.ylabel(str)
#       fig.subplots_adjust(hspace=0.5)
        
        plt.show()
        
    def write_data(self):
        f = open(self.outfile,"w+")
        for i in range(self.num_time):
            str = '%f'%(i*self.dt)
            for k in range(self.num_chan):
                str = str + '\t%f'%self.data[k,i]
            str = str + '\n'
            f.write(str)
        f.close()
        print("written %s"%self.outfile)
        
    def convert_dat_to_hdf5(self,dir0,dir1,label,analevel,labeln):
        import h5py
        
        path = './'
        outpath = './analysis/'+analevel+'/'
        
        dir=dir0+dir1+'/'
        label0 = dir0+dir1+'_'
        #label_file_in = path+dir+dir0+'-'+label+'_artred'
        label_file_in = path+dir+dir0+'-'+label+'_fs500.0_artred'
        label_file_out = outpath+label0+labeln
        print("dir0:%s  dir1:%s label=%s  analevel=%s "%(dir0,dir1,label,analevel))
        print("label0:%s  label_file_out:%s "%(label0,label_file_out))
        
        infile = label_file_in+'.dat' ## artifact reduced data
        #outfile = label_file_out+'.h5' ## artifact reduced data
        outfile = label_file_out+'_fs500.0_'+'.h5'
        data0 = np.loadtxt(infile)
        #self.data = np.zeros((np.shape(data0)[0],np.shape(data0)[1]-1))
        #self.data=data0[:,1:]
        #self.num_time = np.shape(self.data)[0]
        #self.num_chan = np.shape(self.data)[1]
        #self.dt = data0[1,0]-data0[0,0]
        #self.times = data0[:,0]
        #del data0
        
        #a = np.random.random(size=(100,20))
        h5f = h5py.File(outfile, 'w')
        label_data = 'LFP_anaesthesia_level:'+analevel
        ret = h5f.create_dataset(label_data, data=data0)
        print(ret)
        h5f.close()
        print("converted %s to %s"%(infile,outfile))
        
    def plot_barplots(self,AIS_1,AIS_2,p1,p2):
        
        from matplotlib import ticker
        
        width = 0.25  # the width of the bars
        fig, ax = plt.subplots(2,1)
        #ax[0].axis('off')
        #ax[1].axis('off')
        #ax[2].axis('off')
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['right'].set_visible(False)
        ax[1].spines['top'].set_visible(False)
        ax[1].spines['right'].set_visible(False)
        #ax[2].spines['top'].set_visible(False)
        #ax[2].spines['right'].set_visible(False)
        
        
#       AIS_1 = np.zeros((3,7))# animal #3
#       AIS_1[0,0]=2.0690464598985487
#       AIS_1[0,1]=2.0089225815447023
#       AIS_1[0,2]=2.4390186966492093
#       AIS_1[0,3]=2.2622683869967535
#       AIS_1[0,4]=2.3358639167782655
#       AIS_1[0,5]=2.261648808897635
#       AIS_1[0,6]=2.3219136478192626
#       AIS_1[1,0]=2.1391026866004337
#       AIS_1[1,1]=2.204971855414203
#       AIS_1[1,2]=2.1058514030148707
#       AIS_1[1,3]=2.2118131282069875
#       AIS_1[1,4]=2.397188654382849
#       AIS_1[1,5]=2.380077164208342
#       AIS_1[1,6]=2.506495599702075
#       AIS_1[2,0]=2.3640029667885685
#       AIS_1[2,1]=2.4583910729267098
#       AIS_1[2,2]=2.5884813189427796
#       AIS_1[2,3]=2.5551843341530915
#       AIS_1[2,4]=2.4139472541729603
#       AIS_1[2,5]=2.388160093505393
#       AIS_1[2,6]=2.4989214288667325
        
#       AIS_2 = np.zeros((3,6))# animal #7
#       AIS_2[0,0]=2.389689870252096
#       AIS_2[0,1]=2.3868976848403367
#       AIS_2[0,2]=2.440840514783541
#       AIS_2[0,3]=2.475381133737063
#       AIS_2[0,4]=2.4391648802247765
#       AIS_2[0,5]=2.3999313326329705
#       AIS_2[1,0]=2.2669783340405063
#       AIS_2[1,1]=2.258289986546287
#       AIS_2[1,2]=2.3352610009110792
#       AIS_2[1,3]=2.2805480781750314
#       AIS_2[1,4]=2.1797899891623036
#       AIS_2[1,5]=2.210518320080003
#       AIS_2[2,0]=2.200777187080667
#       AIS_2[2,1]=2.5037081576515243
#       AIS_2[2,2]=2.547505698365228
#       AIS_2[2,3]=2.468714160213641
#       AIS_2[2,4]=2.504895321948341
#       AIS_2[2,5]=2.314519868528584
        
        
        #labels_1 = ['#1','#5', '#6','#7', '#8','#9', '#10','#11', '#12','#13','#14','#16']# animal #1  (non-functional)
        #labels_1 = ['#2', '#3','#4', '#15']# animal #1  (non-functional)
        #labels_2 = ['#4', '#7', '#10', '#11', '#12', '#13'] #animal #2 - non-functional
        #labels_2 = ['#1', '#2', '#3', '#5','#6','#8', '#9','#14','#15','16'] #animal #2 - functional
        
        #labels_1 = ['#4', '#5', '#6', '#9', '#10', '#11', '#12', '#13'] #animal #3 - non-functional
        #labels_1 = ['#1', '#2', '#3', '#7', '#8', '#14', '#15','#16'] #animal #3 - functional
        #labels_2 = ['#2', '#3', '#4', '#5', '#14', '#15'] #animal #7 - non-functional
        #labels_2 = ['#1', '#6', '#7', '#8','#9', '#10', '#11', '#12', '#13','#16'] #animal #7 - functional
        
        #labels_1 = ['#2', '#3','#4', '#5','#6', '#7','#10']# animal #5  (non-functional)
        #labels_1 = ['#1', '#8','#9', '#11','#12', '#13','#14','#15','#16']# animal #5  (non-functional)
        #labels_2 = ['#3', '#4', '#5', '#6', '#7', '#10', '#11', '#12', '#13', '#14'] #animal #6 - non-functional
        #labels_2 = ['#1', '#2', '#8', '#9','#15','#16'] #animal #6 - functional
        
        #labels_1 = ['3','4','5','6','9','10','11','12','13','14']# animal #4 (non-functional)
        labels_1 = ['1','2','7','8','15','16']# animal #4 - remaining (functional)
        #labels_2 = ['3','4','5','6','9','10','11','12','13','14']# animal #4 (non-functional)
        labels_2 = ['1','2','7','8','15','16']# animal #4 - remaining (functional)
        
        
        
        ylabels = [1.0,2.0]
        
        x = np.arange(len(labels_1))  # the electrodes        
        h0 = AIS_1[0,:]
        h1 = AIS_1[1,:]
        h2 = AIS_1[2,:]
        rects01 = ax[0].bar(x - width, h0, width, align='center', ecolor='lightgrey', label='pre - low')
        rects02 = ax[0].bar(x        , h1, width, align='center', ecolor='lightgrey', label='post - low')
        rects03 = ax[0].bar(x + width, h2, width, align='center', ecolor='lightgrey', label='post - high')
        ax[0].axhline(0, color='grey', linewidth=0.8)
        ax[0].xaxis.set_major_locator(ticker.FixedLocator(x))
        #ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(labels_1))
        ax[0].xaxis.set_major_formatter(ticker.FixedFormatter([]))
        ax[0].yaxis.set_major_locator(ticker.FixedLocator(ylabels))
        ax[0].yaxis.set_major_formatter(ticker.FixedFormatter([]))
        maxx = max(np.max(h0),np.max(h1),np.max(h2))
        h = 0.08
        #y0 = maxx*1.15#max(meds)+15*h
        y0 = maxx*1.00#num_sign = len(sign_list)
        #for kk in range(num_sign):
        #sign = sign_list[kk]
        alpha = 0.05/float(np.shape(p1)[1])## Bonferroni
        for k in range(np.shape(p1)[1]):
            if p1[0,k]<alpha:
                x1, x2 = x[k] - width, x[k] - width*0.1
                y = y0+h*0.5#max(h0[k],h1[k])+h
                ax[0].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax[0].text((x1+x2)*.5, y-h*0.5,'*', ha='center', va='bottom', color='k')
                #                       plt.ylim(miny*0.97,maxy*1.2) 
            if p1[1,k]<alpha:
                x1, x2 = x[k] - width, x[k] + width
                y = y0+5*h#max(h0[k],h2[k])+h*0.8
                ax[0].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax[0].text((x1+x2)*.5, y-h*0.5,'*', ha='center', va='bottom', color='k')
            if p1[2,k]<alpha:
                x1, x2 = x[k]+width*0.1 , x[k] + width
                y = y0+3*h#max(h1[k],h2[k])+h*1.4
                ax[0].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax[0].text((x1+x2)*.5, y-h*0.5,'*', ha='center', va='bottom', color='k')    
        
        print("AIS_2:",np.shape(AIS_2))
        print("p2:",np.shape(p2))
        x = np.arange(len(labels_2))  # the electrodes        
        h0 = AIS_2[0,:]
        h1 = AIS_2[1,:]
        h2 = AIS_2[2,:]
        rects01 = ax[1].bar(x - width, h0, width, align='center', ecolor='lightgrey', label='pre - low')
        rects02 = ax[1].bar(x        , h1, width, align='center', ecolor='lightgrey', label='post - low')
        rects03 = ax[1].bar(x + width, h2, width, align='center', ecolor='lightgrey', label='post - high')
        ax[1].axhline(0, color='grey', linewidth=0.8)
        ax[1].xaxis.set_major_locator(ticker.FixedLocator(x))
        #ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(labels_2))
        ax[1].xaxis.set_major_formatter(ticker.FixedFormatter([]))
        ax[1].yaxis.set_major_locator(ticker.FixedLocator(ylabels))
        ax[1].yaxis.set_major_formatter(ticker.FixedFormatter([]))
        maxx = max(np.max(h0),np.max(h1),np.max(h2))
        print(maxx)
        h = 0.08
        y0 = maxx*1.00#max(meds)+15*h
        alpha = 0.05/float(np.shape(p2)[1]) # Bonferroni
        for k in range(np.shape(p2)[1]):
            if p2[0,k]<alpha:
                x1, x2 = x[k] - width, x[k] - width*0.1
                y = y0+h*0.5
                ax[1].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax[1].text((x1+x2)*.5, y-h*0.5,'*', ha='center', va='bottom', color='k')
                #                       plt.ylim(miny*0.97,maxy*1.2) 
            if p2[1,k]<alpha:
                x1, x2 = x[k] - width, x[k] + width
                y = y0+5*h#9*h
                ax[1].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax[1].text((x1+x2)*.5, y-h*0.5,'*', ha='center', va='bottom', color='k')
            if p2[2,k]<alpha:
                x1, x2 = x[k]+width*0.1 , x[k] + width
                y = y0+3*h#5*h
                ax[1].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')
                ax[1].text((x1+x2)*.5, y-h*0.5,'*', ha='center', va='bottom', color='k')    
        
        
        plt.show()
        #quit()
#       alpha = 0.05
#       maxx0 = max(h0)
#       maxx1 = max(h1)
#       for k in range(num_list):
#           if p_list[k][0][0]<alpha:
#               val = h0[k]*1.1
#               if val<0:
#                   val = 0.1*maxx0
#               ax[0].text(x[k]-width, val,'*', ha='center', va='bottom', color='k') 
#           if p_list[k][1][0]<alpha:
#               val = h1[k]*1.1
#               if val<0:
#                   val = 0.1*maxx1
#               ax[0].text(x[k], val,'*', ha='center', va='bottom', color='k')
#               
        