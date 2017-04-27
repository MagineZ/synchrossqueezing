import numpy as np
import math

def sqCWTbase(t, x, freqlow, freqhigh, alpha, opts, Smooth, Hemi):
    """
    %
    % Synchrosqueezing transform ver 0.5 (2015-03-09)
    % You can find more information in 
    %	http://sites.google.com/site/hautiengwu/
    %
    % Example: [~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
    %	time: 	time of the signal
    %	xm: 	the signal to be analyzed
    %	[lowfreq, highfreq]: the frequency range in the output time-frequency representation. For the sake of computational efficiency.
    %	alpha:	the frequency resolution in the output time-frequency representation
    %	opts:	parameters for the CWT analysis. See below
    %	tfr/tfrtic:	the CWT and its scale tic
    %	tfrsq/tfrsqtic: the SST-CWT and its frequency tic
    %
    % by Hau-tieng Wu v0.1 2011-06-20 (hauwu@math.princeton.edu)
    %		  v0.2 2011-09-10
    %		  v0.3 2012-03-03
    %		  v0.4 2012-12-12
    %		  v0.5 2015-03-09
    %% you can play with these 4 parameters, but the results might not be
    %% that different if the change is not crazy
    """
    Ex = np.mean(np.power(abs(x),2))
    Gamma = 1.0e-8*Ex  
    """% originally it was 1e-6*Ex"""    

    """%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    dt = t[1] - t[0]

    
    """%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """%% Continuous wavelet transform""" 
    tfr, tfrtic = CWT(t, x, opts)
    Dtfr = (-1j/2/math.pi/dt)*np.concatenate((tfr[:,1:] - tfr[:,0:-1], np.array([tfr[:,-1]-tfr[:,-2]]).T),axis = 1) 
    Dtfr[(abs(tfr) < Gamma)] = float('nan')
    omega = Dtfr/tfr

    """%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """%% Synchro-squeezing transform"""

    tfrsq, tfrsqtic = SQ(tfr, omega, freqlow, freqhigh, alpha, Smooth, Hemi)
    tfr = tfr
    tfrsq = tfrsq
    return tfr, tfrtic, tfrsq, tfrsqtic    




"""%====================================================================="""
"""%% function for CWT-based SST"""
def SQ(tfd, omega, tfrsqticlow, tfrsqtichigh, alpha, Smooth, Hemi):
    
    nvoice = 32
    scale = 2
    omega = abs(omega)
    nscale,n = tfd.shape
    nalpha = np.floor((tfrsqtichigh - tfrsqticlow)/alpha)
    tfrsq = np.zeros((int(nalpha),n),dtype=complex)
    tfrsqtic = np.arange(1,nalpha+1)*alpha + tfrsqticlow
    ntfrsqtic = len(tfrsqtic)
    
    Mid = round(ntfrsqtic/2)-1
    Delta = 20*np.power(tfrsqtic[1]-tfrsqtic[0],2) 
    weight = np.exp(-np.power(tfrsqtic[Mid-10:Mid+11]-tfrsqtic[Mid],2)/Delta)
    weight = weight / sum(weight) 
    weightIDX = np.arange(Mid-10,Mid+11) - Mid

    for b in range(1,n+1):
        for kscale in range(1,nscale+1):
            qscale = scale*np.power(2,(kscale/nvoice))
            
            if (np.isfinite(omega[kscale-1,b-1]) and (omega[kscale-1,b-1]>0)):
                k = int(np.floor( ( omega[kscale-1,b-1] - tfrsqticlow )/ alpha )+1)

                if (np.isfinite(k) and (k > 0) and (k < ntfrsqtic-1)):
                    ha = tfrsqtic[k]-tfrsqtic[k-1]


                    if Smooth:
                        IDXb = np.where(np.logical_and((k+weightIDX < ntfrsqtic-1) , (k+weightIDX > 0)))
                        IDXa = k+weightIDX[IDXb]
                        IDXa = IDXa.astype(int)
                        
                        if Hemi:
                            if (tfd[kscale-1,b-1].real > 0):
                                tfrsq[IDXa-1,b-1] = tfrsq[IDXa-1,b-1] + weight[IDXb]*np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                            else:
                                tfrsq[IDXa-1,b-1] = tfrsq[IDXa-1,b-1] - weight[IDXb]*np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                    	
                        else:
                            tfrsq[IDXa-1,b-1] = tfrsq[IDXa-1,b-1] + np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                	
                    else:
                        if Hemi:
                            if (tfd[kscale-1,b-1].real>0):
                                tfrsq[k-1,b-1] = tfrsq[k-1,b-1] + np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                            else:
                                tfrsq[k-1,b-1] = tfrsq[k-1,b-1] - np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                                     
                        else:
                            tfrsq[k-1,b-1] = tfrsq[k-1,b-1] + np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                                 
    return tfrsq, tfrsqtic