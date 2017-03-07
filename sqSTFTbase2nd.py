import numpy as np
import math

def sqSTFTbase2nd(x, lowFreq, highFreq, alpha, tDS, h, Dh, DDh, online):
    """
    % Synchrosqueezing modifed from tfrrsp.m, by Hau-tieng Wu, 2013 
    %
    %	computes the STFT and its SST
    %
    %   Example:
    %
    %	Hz = 32 ; t=[1/Hz:1/Hz:32]' ;
    %	x=cos(2*pi*(4*t+cos(t/2))) ;
    %	[h, Dh] = hermf(71, 1, 6) ;
    %		%% get the TF representation of the signal x with the frequency range
    %		%% [0,1, 0.4]*Hz with the frequency resolution 0.001*Hz
    %	[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0.1, 0.4, 0.001, 1, h', Dh');
    %	imageRTF(t, tfrsqtic*Hz, abs(tfrsq)) ;
    %
    %		%% the first version can be recovered by
    %   [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0, 0.5, 0.5/length(x), 1, h', Dh');
    %
    %
    %======================================
    %	X     : analysed signal.
    %	[lowFreq, highFreq] : frequency range \subset [0 0.5]. For the sake of computational efficiency.
    %	alpha : the resolution in the frequency axis
    %	tDS   : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
    %	H     : frequency smoothing window, H(0) being forced to 1
    %   DH    : differentiation of H	
    %	TFR   : STFT
    %	TFRSQ : synchrosqueezed STFT 
    %
    %	F. Auger, May-July 1994, July 1995.
    %	Copyright (c) 1996 by CNRS (France).
    %
    %	------------------- CONFIDENTIAL PROGRAM -------------------- 
    %	This program can not be used without the authorization of its
    %	author(s). For any comment or bug report, please send e-mail to 
    %	f.auger@ieee.org 
    """
    x=np.array([x])
    xrow,xcol = x.shape
    t = np.arange(1,xcol+1)
    tLen = len(np.arange(1,xcol+1,tDS))

    """ for tfr """
    N = len(np.arange(-0.5+alpha,0.5,alpha))+1
    """ for tfrsq """
    Lidx = int(np.ceil( (N/2)*(lowFreq/0.5) )+1) 
    Hidx = int(np.floor( (N/2)*(highFreq/0.5)))
    fLen = Hidx - Lidx + 1
    """===================================================================="""
    """check input signals"""
    if (xrow!=1):
        raise ValueError('X must have only one row')
    elif (highFreq > 0.5):
        raise ValueError('TopFreq must be a value in [0, 0.5]')
    elif (tDS < 1) or (tDS%1): 
        raise ValueError('tDS must be an integer value >= 1')

    hrow,hcol = h.shape
    Lh = int((hcol-1)/2) 
    if ((hrow!=1)or((hcol%2)==0)):
        raise ValueError('H must be a smoothing window with odd length')

    ht = np.arange(-Lh,Lh+1)
    """===================================================================="""
    """run STFT and reassignment rule"""
    if online: 
        tfr = np.zeros((100, int(N/2)),dtype=complex)     
        """ for h"""
        tfrsq = np.zeros((100, fLen),dtype=complex)
        tfrsq2nd = np.zeros((100, fLen),dtype=complex)
    else: 
        tfr = np.zeros((tLen, int(N/2)),dtype=complex) 
        """for h"""
        tfrsq = np.zeros((tLen, fLen),dtype=complex)
        tfrsq2nd = np.zeros((tLen, fLen),dtype=complex) 

    tfrtic = np.linspace(0, 0.5, N/2)
    tfrsqtic = np.linspace(lowFreq, highFreq, fLen)
    
    """Ex = mean(abs(x(min(t):max(t))).^2);"""
    Ex = np.mean(np.square(np.absolute(x[0,np.amin(t):np.amax(t)+1])))
    Threshold = math.pow(10,-8)*Ex  
    """% originally it was 1e-6*Ex"""

    for tidx in range(1,tLen+1):
        """ti is the current time"""
        ti = t[(tidx-1)*tDS]
        """tau is the relevant index associated with ti"""
        tau = (np.arange(-np.amin([round(N/2)-1,Lh,ti-1]),np.amin([round(N/2)-1,Lh,xcol-ti])+1)).astype(int)
        """indices is the absolute index in the "evaluation window" """
        indices= ((N+tau)%N+1).astype(int)
        norm_h=np.linalg.norm(h[:,(Lh+tau).astype(int)])
        tf0 = np.zeros((1, N),dtype=complex)
        tf1 = np.zeros((1, N),dtype=complex)
        tf2 = np.zeros((1, N),dtype=complex)
        tfx0 = np.zeros((1, N),dtype=complex)
        tfx1 = np.zeros((1, N),dtype=complex)
        tf0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])/norm_h
        tf1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])/norm_h
        tf2[:,indices-1] = x[:,ti+tau-1]*np.conjugate(DDh[:,Lh+tau])/norm_h
        tfx0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])*ht[Lh+tau]/norm_h
        tfx1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])*ht[Lh+tau]/norm_h
        tf0 = np.fft.fft(tf0);tf0 = tf0[:,0:int(N/2)]
        tf1 = np.fft.fft(tf1);tf1 = tf1[:,0:int(N/2)]
        tf2 = np.fft.fft(tf2);tf2 = tf2[:,0:int(N/2)]
        tfx0 = np.fft.fft(tfx0);tfx0 = tfx0[:,0:int(N/2)]
        tfx1 = np.fft.fft(tfx1);tfx1 = tfx1[:,0:int(N/2)]	
        """% get the first order omega"""
        omega = np.round(N*np.imag(tf1/tf0)/(2.0*math.pi))
        """% get the 2nd order omega"""
        omega2nd = np.round(N * np.imag(tf1/tf0 - (tf0*tf2-tf1*tf1)/(tfx1*tf0-tfx0*tf1)*tfx0/tf0)/(2.0*math.pi))
        
        
        sst = np.zeros((1,fLen),dtype=complex)
        sst2nd = np.zeros((1,fLen),dtype=complex)
        
        for jcol in range(1,int(N/2)+1):
            if abs(tfr[0,jcol-1]) > Threshold:
                
                jcolhat = int(jcol - omega[:,jcol-1])
                jcolhat2nd = int(jcol - omega2nd[:,jcol-1])
                
                if (jcolhat < Hidx+1) and (jcolhat >= Lidx):
                    sst[0,jcolhat-Lidx] = sst[0,jcolhat-Lidx]+tf0[0,jcol-1]
                if (jcolhat2nd < Hidx+1) and (jcolhat2nd >= Lidx):
                    sst2nd[0,jcolhat2nd-Lidx] = sst2nd[0,jcolhat2nd-Lidx]+tf0[0,jcol-1]
                
                
        if online:
            tfr[0:99,:] = tfr[1:100,:]
            tfrsq[0:99,:] = tfrsq[1:100,:]
            tfrsq2nd[0:99,:] = tfrsq2nd[1:100,:]
            tfr[99,:] = tf0[0,0:int(N/2)]
            tfrsq[99,:] = sst
            tfrsq2nd[99,:] = sst2nd
        
        else:
            tfr[tidx-1,:] = tf0[0,0:int(N/2)]
            tfrsq[tidx-1,:] = sst
            tfrsq2nd[tidx-1,:] = sst2nd

    return np.transpose(tfr), tfrtic, np.transpose(tfrsq), np.transpose(tfrsq2nd), tfrsqtic

