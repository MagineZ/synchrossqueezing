"""
import all needed packages
"""
import numpy as np
import math as math
"""
functions
"""

def hermf(N,M,tm):
    """
    hermf.m
    computes a set of orthonormal Hermite functions 
    (for use with tfrrsp_h.m) 
    
    input  : - N : number of points (must be odd) 
             - M : maximum order 
             - tm : half time support (>= 6 recommended) 
             
    output : - h : Hermite functions (MxN) 
             - Dh : H' (MxN) 
             - tt : time vector (1xN) 
    """
    
    dt = 2*tm/(N-1) ; 
    tt = np.linspace(-tm,tm,N) ; 
    g = np.exp(-np.square(tt)/2) ; 
    P = [np.ones(N)] 
    P = np.append(P, [2*tt], axis=0)

    for k in range(2,M+1):
        P = np.append(P,[2*np.multiply(tt,P[k-1]) - 2*(k-1)*P[k-2]],axis=0) ; 
     
    Htemp = [np.multiply(P[0],g)/math.sqrt(math.sqrt(math.pi)*math.gamma(1))*math.sqrt(dt)] ; 
    for k in range(1,M+1): 
        Htemp = np.append(Htemp, [np.multiply(P[k],g)/math.sqrt(math.sqrt(math.pi)
                           *math.pow(2,k)*math.gamma(k+1))*math.sqrt(dt)],axis = 0); 
       

    h = Htemp[0:M,:];          
    Dh = [(np.multiply(tt,Htemp[0]) - math.sqrt(2)*Htemp[1])*dt] ;
    for k in range(1,M): 
        Dh = np.append(Dh,[(np.multiply(tt,Htemp[k]) - math.sqrt(2*(k+1))*Htemp[k+1])*dt] ,axis=0);
      
    
    return h,Dh,tt


def dwindow(h):
    """
    %DWINDOW Derive a window.
    %	DH=DWINDOW(H) derives a window H.
    %
    %	Example : 
    %	 plot(dwindow(tftb_window(210,'hanning')))
    %
    %	See also WINDOW.
    
    %	F. Auger, August 1994, July 1995.
    %	Copyright (c) 1996 by CNRS (France).
    %
    %  This program is free software; you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation; either version 2 of the License, or
    %  (at your option) any later version.
    %
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    %
    %  You should have received a copy of the GNU General Public License
    %  along with this program; if not, write to the Free Software
    %  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
    """

    [hrow,hcol]=h.shape 
    if not (hrow==1):
        raise ValueError('h must have only one row')
    
    Lh = int((hcol-1)/2)
    step_height = (h[0,0]+h[0,-1])/2
    ramp = (h[0,-1]-h[0,0])/(hcol-1)
    h2 = np.append([0],h-step_height-ramp*np.arange(-Lh,Lh+1))
    h2 = np.append(h2,[0])
    Dh = (h2[2:hcol+2]-h2[0:hcol])/2+ramp
    Dh[0]+=step_height 
    Dh[-1]-=step_height
    return np.array(Dh)

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



def ConceFT_sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Second, Smooth, Hemi):
    """
    Usage: 
        [tfrsq, ConceFT, tfrsqtic] = sqSTFT(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp, MT)
        MT = 1: ordinary SST; MT > 1: ConceFT
        alpha: resolution in the frequency axis
        WinLen, dim, supp: for hermf.m
    Example:
        [tfrsq, ConceFT, tfrsqtic] = ConceFt_sqSTFT_C([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6, 10);
    """ 
    """       
    Multitapering
    generate the window for short time Fourier transform (STFT)
    """
    h,Dh,_ = hermf(WinLen, dim, supp)
                  

    print('Run ordinary STFT-SST (Smooth = ',Smooth,', Hemi = ',Hemi,')\n')
    tfr, tfrtic, tfrsq, tfrsqtic = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, np.conj(np.array([h[0,:]])), np.conj(np.array([Dh[0,:]])), Smooth, Hemi)
    ConceFT = tfrsq 


    if MT > 1:
        print('Complex sphere\n') 
        print('STFT-ConceFT total (Smooth = ',Smooth,', Hemi = ',Hemi,'): ',MT,'; now:     ')
        
        for ii in range(1, MT+1):
            print('\b\b\b\b') 
            print('%4d' % ii) 
            rv = np.random.randn(1, dim) + np.sqrt(-1+0j)*np.random.randn(1, dim)
            rv = rv/np.linalg.norm(rv)
            rh = np.dot(rv , h)  
            rDh = np.dot(rv , Dh)
            
            if not Second:
                _, _, tfrsqX, tfrsqtic = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, np.conj(rh), np.conj(rDh), Smooth, Hemi)
            else:
                _, _, _, tfrsqX, tfrsqtic = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, np.conj(rh), np.conj(rDh), dwindow(np.conj(rDh)), 0);
		
            ConceFT = ConceFT + tfrsqX 
                 
    

        ConceFT = ConceFT/(MT+1)
        print('\n') 

    return tfr, tfrtic, tfrsq, ConceFT, tfrsqtic

def imageSQ(ax, t, ytic, M, Qv, cm = 'Greys'):
    
    q = np.percentile(M,Qv)
    """ truncate the upper bound"""
    M[np.where(M>q)]=q
    """M = M/q"""
    """
    truncate the lower bound
    %m = quantile(Q, 0.002) ;
    %M(find(M<m)) = m ;
    """
    ax.pcolorfast(t, ytic, M, cmap = cm)