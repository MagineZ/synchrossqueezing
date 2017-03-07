# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 16:35:23 2017

@author: Pei-Chun Su
"""
"""
x=xm
lowFreq=0
highFreq=0.5
alpha=FrequencyAxisResolution
tDS=1
Smooth=0
Hemi=0
WinLen=WindowLength
dim= NoWindowsInConceFT
supp=WindowBandwidth
h,Dh,_ = hermf(WinLen, dim, supp) 
h=h[0,:]
Dh=Dh[0,:] 
rv = np.random.randn(1, dim) + np.sqrt(-1+0j)*np.random.randn(1, dim)
rv = rv/np.linalg.norm(rv)
rh = np.dot(rv , h)  
rDh = np.dot(rv , Dh)
            
"""
import numpy as np
import math
def sqSTFTbase(x, lowFreq, highFreq, alpha, tDS, h, Dh, Smooth, Hemi):
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
    Lidx = round( (N/2)*(lowFreq/0.5) )+1 
    Hidx = round( (N/2)*(highFreq/0.5))
    fLen = Hidx - Lidx + 1

    "=========================================================================="
    "check input signals"
   
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

    "=========================================================================="
    "run STFT and reassignment rule"
    tfr = np.zeros((tLen, int(N/2)),dtype=complex)
    tfrtic = np.linspace(0, 0.5, int(N/2)) 
    tfrsq = np.zeros((tLen, fLen),dtype=complex)
    tfrsqtic = np.linspace(lowFreq, highFreq, fLen)


    """Ex = mean(abs(x(min(t):max(t))).^2);"""
    Ex = np.mean(np.square(np.absolute(x)))
    Threshold = math.pow(10,-8)*Ex  
    """% originally it was 1e-6*Ex"""

    Mid = round(len(tfrsqtic)/2)
    Delta = 20*np.square(tfrsqtic[1]-tfrsqtic[0]) 
    weight = np.exp(-np.square(tfrsqtic[Mid-11:Mid+10]-tfrsqtic[Mid-1])/Delta)
    weight = weight/np.sum(weight)
    weightIDX = np.arange(Mid-10,Mid+11,1) - Mid

    for tidx in range(1,tLen+1):
        ti = t[(tidx-1)*tDS] 
        tau = (np.arange(-np.amin([round(N/2)-1,Lh,ti-1]),np.amin([round(N/2)-1,Lh,xcol-ti])+1)).astype(int)
        norm_h=np.linalg.norm(h[:,(Lh+tau).astype(int)])
        """%norm_h = h(Lh+1)""" 
        indices= ((N+tau)%N+1).astype(int)
        tf0 = np.zeros((1, N),dtype=complex)
        tf1 = np.zeros((1, N),dtype=complex) 
        tf0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])/norm_h
        tf1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])/norm_h
        tf0 = np.fft.fft(tf0);tf0 = tf0[:,0:int(N/2)];tf1 = np.fft.fft(tf1);tf1 = tf1[:,0:int(N/2)]
        """% get the first order omega"""
        omega = np.zeros(tf1.shape)
        _,avoid_warn = np.nonzero(tf0)
        omega[:,avoid_warn]= np.round(np.imag(N*np.divide(tf1[:,avoid_warn],tf0[:,avoid_warn])/(2.0*math.pi)))
        sst = np.zeros((1,fLen),dtype=complex) 

        for jcol in range(1,int(N/2)+1):
            if abs(tfr[0,jcol-1]) > Threshold:
                jcolhat = int(jcol - omega[:,jcol-1])
                """%jcolhat = rem(rem(jcolhat-1,N)+N,N)+1;"""
                
                """if(jcolhat < Hidx + 1) & (jcolhat >= Lidx):
                sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) + tf0(jcol) ;
	    	      end"""     
                if (jcolhat <= Hidx) & (jcolhat >= Lidx):   
                    """%IDXa = unique(min(Hidx, max(Lidx, jcolhat-Lidx+1+weightIDX))) ; """
                    if Smooth:
                        IDXb = np.where((jcolhat-Lidx+1+weightIDX <= Hidx) & (jcolhat-Lidx+1+weightIDX >= Lidx))
                        IDXa = (jcolhat-Lidx+1+weightIDX[IDXb]).astype(int)
                        
                        if Hemi:
                            if np.real(tf0[0,jcol-1])>0:
                                sst[0,IDXa-1]=sst[0,IDXa-1]+tf0[:,jcol-1]*weight[IDXb]
                                """%sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) + tf0(jcol) ;
                                %sst(jcolhat,icol) = rtfr(jcolhat,icol) + tfr(jcol,icol) ;"""
                            else:
                                sst[0,IDXa-1]=sst[0,IDXa-1]-tf0[:,jcol-1]*weight[IDXb]
                                """%sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) - tf0(jcol) ;
                                %sst(jcolhat,icol) = rtfr(jcolhat,icol) - tfr(jcol,icol) ;"""
                        else:
                            sst[0,IDXa-1]=sst[0,IDXa-1]+tf0[:,jcol-1]*weight[IDXb]
                	
                    else:
                        if Hemi:
                            if np.real(tf0[0,jcol-1])>0:
                                sst[0,jcolhat-Lidx]=sst[0,jcolhat-Lidx]+tf0[0,jcol-1]
                                """%sst(jcolhat,icol) = rtfr(jcolhat,icol) + tfr(jcol,icol)"""
                            else:
                                sst[0,jcolhat-Lidx]=sst[0,jcolhat-Lidx]-tf0[0,jcol-1]
                                """%sst(jcolhat,icol) = rtfr(jcolhat,icol) - tfr(jcol,icol)"""
                        else:
                            sst[0,jcolhat-Lidx]=sst[0,jcolhat-Lidx]+tf0[0,jcol-1] 
		   			
        tfr[tidx-1, :] = tf0[0,0:int(N/2)]
        tfrsq[tidx-1,:]=sst
    return np.transpose(tfr), tfrtic, np.transpose(tfrsq), tfrsqtic     

	



