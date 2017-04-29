import numpy as np
def rsSTFTbase(x, lowFreq, highFreq, alpha, tDS, h, Dh, Causality):
    """
    %
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

    """% for tfr"""
    N = len(np.arange(-0.5+alpha,0.5,alpha))+1

    """% for tfrsq"""
    Lidx = round( (N/2)*(lowFreq/0.5)) + 1 
    Hidx = round( (N/2)*(highFreq/0.5)) 
    fLen = Hidx - Lidx + 1 



    """%===================================================================="""
    """%% check input signals"""
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

    Th = h*np.arange(-Lh,Lh+1)
    Dt = 1
    
    "=========================================================================="
    "run STFT and reassignment rule"
    tfr = np.zeros((tLen, int(N/2)),dtype=complex)
    tfrtic = np.linspace(0, 0.5, int(N/2)) 
    tfrrs = np.zeros((tLen, fLen),dtype=complex)
    tfrrstic = np.linspace(lowFreq, highFreq, fLen)
    
    Ex = np.mean(np.square(np.absolute(x[0,np.amin(t):np.amax(t)+1])))
    Threshold = 10e-8*Ex
    """% originally it was 1e-6*Ex"""


    for tidx in range(1,tLen+1):
        ti = t[(tidx-1)*tDS]
        tau = (np.arange(-np.amin([round(N/2)-1,Lh,ti-1]),np.amin([round(N/2)-1,Lh,xcol-ti])+1)).astype(int)
        norm_h=np.linalg.norm(h[:,(Lh+tau).astype(int)])
        indices= ((N+tau)%N+1).astype(int)
        tf0 = np.zeros((1, N),dtype=complex)
        tf1 = np.zeros((1, N),dtype=complex)
        tf2 = np.zeros((1, N),dtype=complex)
        tf0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])/norm_h
        tf1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])/norm_h
        tf2[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Th[:,Lh+tau])/norm_h
        tf0 = np.fft.fft(tf0);tf0 = tf0[:,0:int(N/2)]
        tf1 = np.fft.fft(tf1);tf1 = tf1[:,0:int(N/2)]
        tf2 = np.fft.fft(tf2);tf2 = tf2[:,0:int(N/2)]

        """% get the first order omega"""
        omega = np.zeros(tf1.size)
        omegaT = np.zeros(tf1.size)
        avoid_warn_x,avoid_warn_y = np.where(tf0!=0)
        omega[avoid_warn_y] = np.round(np.imag(N*tf1[avoid_warn_x,avoid_warn_y]/tf0[avoid_warn_x,avoid_warn_y]/(2.0*np.pi)))
        omegaT[avoid_warn_y] = np.round(np.real(tf2[avoid_warn_x,avoid_warn_y]/tf0[avoid_warn_x,avoid_warn_y]/Dt))

        """%rs = np.zeros((1,fLen) ;"""

        for jcol in range(1,int(N/2)+1):
            if abs(tfr[0,jcol-1]) > Threshold:
                """%% to keep the causality, use max(tidxhat,tidx-Lh) instead of max(tidxhat,1)"""
                tidxhat = int(tidx + omegaT[jcol-1])
                if Causality:
                    tidxhat = min(max(max(tidxhat, tidx-Lh), 1), tLen)
                else:
                    tidxhat = min(max(tidxhat, 1), tLen)

                jcolhat = int(jcol - omega[jcol-1]) 
                """jcolhat = rem(rem(jcolhat-1%N)+N,N)+1;"""
                if(jcolhat <= Hidx) and (jcolhat >= Lidx):
                    tfrrs[tidxhat-1,jcolhat-Lidx] = tfrrs[tidxhat-1,jcolhat-Lidx] + tf0[0,jcol-1]


        tfr[tidx-1,:] = tf0[0:int(N/2)] 


    return np.transpose(tfr), tfrtic, np.transpose(tfrrs), tfrrstic