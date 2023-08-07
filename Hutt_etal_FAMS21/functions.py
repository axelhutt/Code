import numpy as np
import scipy.special as special
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def dglsyst_test(V,W,P):
    help1 =  P.F*S1(V,P)-P.M*S2(W,P)+P.I10
    help2 = - P.F*S2(W,P)+P.M*S1(V,P)+P.I20
    return help1,help2

def dglsyst(V,W,P):
    help1 = -V + P.F*S1(V,P)-P.M*S2(W,P)+P.I10
    help2 = -W - P.F*S2(W,P)+P.M*S1(V,P)+P.I20
    return help1,help2

def dglsyst_linear(x,y,V,W,P):
    help1 = -x + P.F*S1_deriv(V,P)*x-P.M*S2_deriv(W,P)*y
    help2 = -y - P.F*S2_deriv(W,P)*y+P.M*S1_deriv(V,P)*x
    return help1,help2

def S(V):
    return 1.0/(1.0+np.exp(-200.0*(V)))


def S1(V,P):
    #return 1.0/(1.0+np.exp(-P.c*(V-P.theta))) 
    #return P.S10*0.5*(1.0-special.erf((P.theta-(V-P.mu1))/(np.sqrt(2*P.noisevariance_u))))
    variance_1 = P.noisevariance_u_1*P.alpha
    variance_2 = P.noisevariance_u_2*P.alpha
    #print("S1 &&&&&&&&& noisevariance:  %f    %f   "%(np.sqrt(2*variance_1*P.alpha),P.noisevariance_u_1*P.alpha))
    return P.S10*0.5*( P.p1*(1.0-special.erf((P.theta-(V-P.mu1))/(np.sqrt(2*variance_1*P.alpha))))\
                      +P.p2*(1.0-special.erf((P.theta-(V-P.mu2))/(np.sqrt(2*variance_2*P.alpha)))))

#     theory[l]=0.5*(percentage*(1.0-special.erf(-(V[l]-mu1)/(np.sqrt(2*D_1/tau))))+ \
#              (1.0-percentage)*(1.0-special.erf(-(V[l]-mu2)/(np.sqrt(2*D_2/tau)))))                        
                        
def S1_deriv(V,P):
    #return 1.0/(1.0+np.exp(-P.c*(V-P.theta))) 
    print("P.noisevariance_u_1=%f P.noisevariance_u_2=%f"%(P.noisevariance_u_1*P.alpha,P.noisevariance_u_2*P.beta))
    variance_1 = P.noisevariance_u_1*P.alpha
    variance_2 = P.noisevariance_u_2*P.alpha
    #return P.S10*np.exp(-np.power(V-P.mu1-P.theta,2)/(2.0*P.noisevariance_u))/(np.sqrt(2*np.pi*P.noisevariance_u))
    return P.S10*( P.p1*np.exp(-np.power(V-P.mu1-P.theta,2)/(2.0*variance_1*P.alpha))/
                            (np.sqrt(2*np.pi*variance_1*P.alpha))\
                  +P.p2*np.exp(-np.power(V-P.mu2-P.theta,2)/(2.0*variance_2*P.alpha))/
                            (np.sqrt(2*np.pi*variance_2*P.alpha)))

def S2(V,P):
    #return 1.0/(1.0+np.exp(-P.c2*(V-P.theta))) 
    return 0.5*(1.0-special.erf((P.theta-V)/(np.sqrt(2*P.noisevariance_v*P.beta))))
def S2_deriv(V,P):
    #print("P.noisevariance_v=%f"%(P.noisevariance_v))
    return np.exp(-np.power(V-P.theta,2)/(2.0*P.noisevariance_v*P.beta))/(np.sqrt(2*np.pi*P.noisevariance_v*P.beta))

def stat(V,W,P):
    help1 = V-(P.R1*P.F*S1(V-P.z_e,P)-P.R1*P.M*S2(W,P)+P.R1*P.I10+P.p1*P.z_e)
    help2 = W-(-P.R2*P.F*S2(W,P)+P.R2*P.M*S1(V-P.z_e,P)+P.R2*P.I20)
    #help2 = W-(-P.F*S2(W,P)+P.R2*P.M*S1(V,P)+P.R2*P.I20)
    #print("**********#################### %f  %f  %f  %f"%(P.R1*P.F*S1(V-P.z_e,P),P.R1,P.F,S1(V-P.z_e,P)))
    return help1,help2

def find_stat(V,W,P):
    num_random = 100000
    deltaV = 1.0
    deltaW = 1.0
    thr1 = 0.005
    thr2 = 0.005
    V_min = V-deltaV
    V_max = V+deltaV
    W_min = W-deltaW
    W_max = W+deltaW
    flag_stop = 0
    for k in range(num_random):
        v=(V_max-V_min)*np.random.uniform()+V_min
        w=(W_max-W_min)*np.random.uniform()+W_min
        g1,g2=stat(v,w,P)
        if np.abs(g1)<thr1 and np.abs(g2)<thr2:
            flag_stop = 1
            break 
    if flag_stop == 1:
        return v,w,g1,g2
    else:
        return P.valthr,P.valthr,g1,g2
    

def find_stat_precise(V,W,P):
    num_random = 50000
    deltaV = 0.01
    deltaW = 0.01
    thr1 = 0.00
    thr2 = 0.00
    V_min = V-deltaV
    V_max = V+deltaV
    W_min = W-deltaW
    W_max = W+deltaW
    flag_stop = 0
    g_min=1000000.0
    g1_min=1000000.0
    g1_min=1000000.0
    for k in range(num_random):
        v=(V_max-V_min)*np.random.uniform()+V_min
        w=(W_max-W_min)*np.random.uniform()+W_min
        g1,g2=stat(v,w,P)
        if np.abs(g1)+np.abs(g2)<g_min : 
            v_min=v
            w_min=w
            g1_min=g1
            g2_min=g2
            g_min=np.abs(g1)+np.abs(g2)
            #print "v=%f w=%f g1=%f g2=%f g_min=%f" %(v,w,g1,g2,g_min)
    return v_min,w_min,g1_min,g2_min

def find_mult_stat(Vmin,Vmax,Wmin,Wmax,P):
    num_random = 500000
    thr1 = 0.2#0.0005#0.005
    thr2 = 0.06#0.0005#0.005
    V_min = Vmin
    V_max = Vmax
    W_min = Wmin
    W_max = Wmax
    flag_stop = 0
    stat_sols = []
    for k in range(num_random):
        v=(V_max-V_min)*np.random.uniform()+V_min
        w=(W_max-W_min)*np.random.uniform()+W_min
        g1,g2=stat(v,w,P)
        if np.abs(g1)<thr1 and np.abs(g2)<thr2:
            flag_stop = 1
            stat_sols.append([v,w,g1,g2]) 
    if flag_stop == 1:
        return stat_sols
    else:
        stat_sols.append([P.valthr,P.valthr,g1,g2])
        return stat_sols

def simulate_solution(V0,W0,P):
    V = V0#+np.random.normal(0.0,1.0)*0.1
    W = W0#+np.random.normal(0.0,1.0)*0.1
    
    S1prime = S1_deriv(V,P)
    S2prime = S2_deriv(W,P)
    a = P.alpha*(-1.0+P.F*S1prime)
    b = P.alpha*P.M*S2prime
    c = P.beta*P.M*S1prime
    d = P.beta*(1.0+P.F*S2prime)
    gamma = np.power(a+d,2)-4.0*b*c
    lambda_r = a-d
    print("variance_u_1=%f  variance_u_2=%f  lambda_r=%f gamma=%f a=%f b=%f c=%f d=%f     S1=%f S2=%f"%(P.noisevariance_u_1,P.noisevariance_u_2,lambda_r,gamma,a,b,c,d,S1prime,S2prime))
    
    num_time = 200000#2000000
    dt = 0.0001#0.0005
    flag_powerspectrum = 1 ###### no power spectrum plotted
    if flag_powerspectrum == 0:
        D = 0.0#1.e-4
        free_period=1.0
        num_free_period = int(free_period/dt)
    else:
        D = 1.e-6
    signal = np.zeros(num_time)
    f = open('./signal.dat','w+')
    for i in range(num_time):
        if flag_powerspectrum == 0 and i%num_free_period == 0:
            V=V+0.0
            W=W+0.0
        noise = np.random.normal(0.0,1.0)
        help1 = -V+(P.R1*P.F*S1(V,P)-P.R1*P.M*S2(W,P)+P.R1*P.I10)
        help2 = -W+(-P.R2*P.F*S2(W,P)+P.R2*P.M*S1(V,P)+P.R2*P.I20)
        #print("i=%d %f %f       %f %f      %f %f %f %f"%(i,help1,help2,V,W,-V,P.R1*P.F*S1(V,P),-P.R1*P.M*S2(W,P),P.R1*P.I10))
        #print("i=%d %f %f       %f %f      %f %f %f %f     -W=%f %f %f %f      F=%f M=%f"%(i,help1,help2,V,W,-V,P.R1*P.F*S1(V,P),-P.R1*P.M*S2(W,P),P.R1*P.I10,-W,-P.F*S2(W,P),P.M*S1(V,P),P.R2*P.I20,P.F,P.M))
        #print("Help1=%f Help2=%f"%(P.alpha*dt*help1,P.beta*dt*help2),P.beta)
        V=V+P.alpha*(dt*help1+np.sqrt(2.*D*dt)*noise)
        W=W+P.beta*(dt*help2)
        signal[i]=V
        str = '%f %f\n'%(dt*float(i),V)
        f.write(str)
        #if i>1000:
        #    break
    f.close()
    print("written signal.dat !!!!!!!")
    if flag_powerspectrum == 1:
        df = 0.5
        nfft = int(1.0/(df*dt))
        if nfft>num_time:
            print ("time segment is too short !!! nfft=%d duration=%d"%(nfft,num_time))
            quit()
        fs = 1.0/dt
        [power0_n,freqs_n]=mlab.psd(signal, NFFT=nfft,Fs=fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=nfft*0.8, pad_to=None, sides='onesided', scale_by_freq=True)		
        num_freqs_n=np.shape(freqs_n)[0]
        df_n=freqs_n[2]-freqs_n[1]
        #print("freqs_n.shape:",np.shape(freqs_n))
        print("min freq=%f      max freq=%f   df=%f"%(freqs_n[0],freqs_n[num_freqs_n-1],df_n))
        fmin_n=15.0
        fmax_n = 60.0
        fmin_num_n=int(fmin_n/df_n)
        fmax_num_n=int(fmax_n/df_n)
        num_n=fmax_num_n-fmin_num_n
        maxpower = np.max(power0_n[fmin_num_n:fmax_num_n])
        maxpower_index = np.argmax(power0_n[fmin_num_n:fmax_num_n])
        maxpower_freq = freqs_n[maxpower_index+fmin_num_n]
        plt.plot(freqs_n[fmin_num_n:fmax_num_n],power0_n[fmin_num_n:fmax_num_n])
        plt.show()
        return maxpower_freq 
    else:
        return -1

def find_multiple_stationary_solutions(V0,W0,Vmin,Vmax,Wmin,Wmax,Imin,Imax,I2min,I2max,deltamu_init,deltamu_final,noisevar0,noisevar1,p0,p1, R20,R21,inputrate_min,inputrate_max,num,P):
    #I10_0=P.I10
    #Imin=P.I10#-0.5
    #Imax=P.I10#+0.5
    #Vmin=-0.05#V0-2.5
    #Vmax=1.0#V0+6.5
    #Wmin=0.9#W0-2.5
    #Wmax=3.7#W0+6.5
    
    #num = 200 # number of grid points in search direction I10
    
    f_stat_stable = open('./stat_sols_stable.dat','w')
    f_stat_unstable = open('./stat_sols_unstable.dat','w')
    f_stat_stable_Hopf = open('./stat_sols_stable_Hopf.dat','w')
    f_stat_unstable_Hopf = open('./stat_sols_unstable_Hopf.dat','w')
    f_stat_stable_node = open('./stat_sols_stable_node.dat','w')
    f_stat_unstable_node = open('./stat_sols_unstable_node.dat','w')
    
    Vavg=-100000.0
    Wavg=-100000.0
    
    flag_simulation = 0
    
    for k in range(num+1):
        P.p = p0+(p1-p0)*float(k)/float(num)
        
        P.R2 = R20+(R21-R20)*float(k)/float(num)
        
        P.I20=I2min+(I2max-I2min)*float(k)/float(num)
        
        P.I10=Imin+(Imax-Imin)*float(k)/float(num)
        
        rate = (inputrate_min+(inputrate_max-inputrate_min)*float(k)/float(num))
        input_mu = P.inputstrength*P.inputtau*rate
        P.mu1 = deltamu_init+(deltamu_final-deltamu_init)*float(k)/float(num)\
                + input_mu
        P.mu2= 0.0
        
        input_var = 0.5*(P.inputstrength**2)*P.inputtau*rate
        P.noisevariance_u_1=noisevar0+(noisevar1-noisevar0)*float(k)/float(num)\
                            + input_var/P.alpha
        noisevariance_u_1_scale = P.noisevariance_u_1*P.alpha
         
        P.z_e = input_mu#/P.alpha
        print("ze=%f  inputstrength=%f  inputtau=%f  rate=%f inputvar=%f"%(P.z_e,P.inputstrength,P.inputtau,rate,input_var))
        
        print("input_var ::::: %f :::  input_mu=%f :::: z_e=%f"%(input_var,input_mu,P.z_e))
        print("variance ::::: %f  "%(noisevariance_u_1_scale))
        
        P.R2 *= P.p
        
        P.R1 = 1.0
        P.R2 = 1.0
        
        maxpowerfreq=-1
        
        if k%5 == 0:
             print( "%d(%d) noisevariance_u_1=%f I10=%f I20=%.3f mu=%.4f p=%.3f" %(k,num,P.noisevariance_u_1*P.alpha,P.I10,P.I20,P.mu1,P.p))
        mult_statsols=find_mult_stat(Vmin,Vmax,Wmin,Wmax,P)
        numsols = len(mult_statsols)
        if k==0:
            #             print mult_statsols
            Vmin_=1000000.0
            Wmin_=1000000.0
            #if (mult_statsols[0])[0]>P.valthr :
            for ll in range(numsols):
                if (mult_statsols[ll])[0]<Vmin_:
                    Vmin_=(mult_statsols[ll])[0]
                if (mult_statsols[ll])[1]<Wmin_:
                    Wmin_=(mult_statsols[ll])[1]
                print( "%d Vmin_=%f" %(ll,Vmin_))
            #print "numsols=%d" %numsols
        for l in range(numsols):
            V,W,g1,g2=mult_statsols[l]
            print("check:"," k=%d V=%f W=%f : noisevar=%f"%(k,V,W,P.noisevariance_u_1),stat(V,W,P))
            if V>P.valthr :
                print("k=%d l=%d"%(k,l))
                #if k==num-2 and l==0:
                if k>int(num/2) and l==0 and flag_simulation==0:
                    #maxpowerfreq = simulate_solution(V,W,P)
                    print("############ maximum power frequency: %f\n"%maxpowerfreq)
                    flag_simulation = 1
                lambda_r,gamma=eigenvalues(V,W,P)
                freq = 0.0
                if gamma>=0:
                    lambda_=0.5*lambda_r+0.5*np.sqrt(gamma)
                    lambda_1=0.5*lambda_r-0.5*np.sqrt(gamma)
                    print( "found: %d (%d)   noisevar=%f rate=%.3f I=%.3f I2=%.3f mu=%.4f p=%.3f R2=%.3f V=%f W=%f g1=%f g2=%f  lambda_=%f lambda_1=%f   gamma=%f      lambda_r=%f" %(l,numsols,P.noisevariance_u_1*P.alpha,rate,P.I10,P.I20,P.mu1,P.p,P.R2,V,W,g1,g2,lambda_,lambda_1,gamma,lambda_r))
                else:
                    lambda_=0.5*lambda_r
                    lambda_1=0.5*lambda_r
                    freq = 0.5*np.sqrt(-gamma)/(2.0*np.pi)
                    print( "found: %d (%d)   noisevar=%f rate=%.3f  I=%.3f mu=%.4f p=%.3f R2=%.3f V=%f W=%f g1=%f g2=%f  lambda_=%f gamma=%f    freq=%f   " %(l,numsols,P.noisevariance_u_1*P.alpha,rate,P.I10,P.mu1,P.p,P.R2,V,W,g1,g2,lambda_,gamma,freq))
                #print( "found: %d (%d)   noisevar=%f V=%f W=%f g1=%f g2=%f  lambda=%f lambda_1=%f   gamma=%f" %(l,numsols,P.noisevariance_u,V,W,g1,g2,lambda_,lambda_1,gamma))
                str = '%f %f %f %f %f  %f %f  %f \n' %\
                            (P.noisevariance_u_1*P.alpha,rate,P.p,V,W,lambda_,lambda_1,gamma)
                if lambda_>0.0:
                    f_stat_unstable.write(str)
                else:
                    f_stat_stable.write(str)
                
                if gamma>0 and lambda_>0.0:
                    f_stat_unstable_node.write(str)
                if gamma>0 and lambda_<=0.0:
                    f_stat_stable_node.write(str)
                if gamma<=0 and lambda_>0.0:
                    str = '%f %f %f  %f %f  %f %f  \n' %\
                            (P.noisevariance_u_1*P.alpha,rate,P.p,V,W,lambda_,freq)
                    f_stat_unstable_Hopf.write(str)
                if gamma<=0 and lambda_<=0.0:
                    str = '%f %f  %f   %f %f  %f %f  \n' %\
                            (P.noisevariance_u_1*P.alpha,rate,P.p,V,W,lambda_,freq)
                    f_stat_stable_Hopf.write(str)
                         
    f_stat_stable.close()
    f_stat_unstable.close()
    f_stat_stable_Hopf.close()
    f_stat_unstable_Hopf.close()
    f_stat_stable_node.close()
    f_stat_unstable_node.close()

    return Vmin_,Wmin_

#def compute_gamma(V,W,P):
#    S1prime = S1_deriv(V,P)
#    S2prime = S2_deriv(W,P)
#    return P.F*P.F*np.power(S1prime+S2prime,2)-4.0*P.M*P.M*S1prime*S2prime

def eigenvalues(V,W,P):
    S1prime = S1_deriv(V,P)
    S2prime = S2_deriv(W,P)
    beta = P.beta/P.p
    a = P.alpha*(-1.0+P.R1*P.F*S1prime)
    b = P.alpha*P.R1*P.M*S2prime
    #a = P.alpha*(-1.0+P.F*S1prime)
    #b = P.alpha*P.M*S2prime
    c = beta*P.R2*P.M*S1prime
    #c = P.beta*P.M*S1prime
    d = P.beta*(1.0+P.R2*P.F*S2prime)
    gamma = np.power(a+d,2)-4.0*b*c
    lambda_r = a-d
    #print("a=%f b=%f c=%f d=%f gamma=%f lambda_r=%f   S1=%F S2=%f"%(a,b,c,d,gamma,lambda_r,S1prime,S2prime))
    #gamma = P.F*P.F*np.power(S1prime+S2prime,2)-4.0*P.M*P.M*S1prime*S2prime
    #lambda_r = 0.5*(P.F*(S1prime-S2prime)-2.0)
    return lambda_r,gamma
   
    

