# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 00:16:51 2017

@author: NTU_Math
"""
import numpy as np
import scipy as sp
import wwhat    

class init:
    motherwavelet = 'Cinfc'
    CENTER = 1
    FWHM = 0.3
    
def CWT(t, x, opts = init):
    """%
    % Continuous wavelet transform ver 0.1
    % Modified from wavelab 85
    %
    % INPUT:
    % OUTPUT:
    % DEPENDENCY:
    %
    % by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
    % by Hau-tieng Wu 2012-12-23 (hauwu@math.princeton.edu)
    %
    """

    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """prepare for the input"""
    from math import sqrt, pi
    from numpy import power
    nvoice = 32;
    scale = 2;
    Oct = 1;
    
    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """start to do CWT"""
    
    n = len(x)
    """%% assume the original signal is on [0,L]."""
    """%% assume the signal is on [0,1]. Frequencies are rescaled to xi/L"""
    if (n%2==0):
        xi = np.concatenate((np.arange(0,n/2+1) , np.arange(-n/2+1,0)))
    else:
        xi = np.concatenate((np.arange(0,n/2) , np.arange(-n/2+1,-1)))
    xhat = np.fft.fft(x)
    
    noctave = np.floor(np.log2(n)) - Oct
    tfr = np.zeros((n,int(nvoice*noctave)),dtype=complex)
    kscale = 1
    tfrtic = np.zeros(int(nvoice*noctave))
    for jj in range(1 , int(nvoice*noctave)+1):
        tfrtic[jj-1] = scale*power(2,jj/nvoice)



    if (opts.motherwavelet=='morse-b' or opts.motherwavelet=='morse-c'):
        dim = opts.dim
        uFix = np.random.randn(dim,1)
        uFix = -sp.linalg.orth(uFix)


    for jo in range(1,int(noctave)+1): 
        """# of scales"""
        for jv in range(1,nvoice+1):
            qscale = scale * power(2,jv/nvoice)
            omega =  xi/qscale            
             
            if (opts.motherwavelet == 'morse'):
                windowq = wwhat(np.conj(omega), opts.beta, opts.gam, opts.k)
                windowq = np.conj(windowq.T)
                
            elif (opts.motherwavelet == 'morse-b'):
                u = np.random.randn(opts.dim,1)
                u = -sp.linalg.orth(u)
                W = np.zeros((len(omega), opts.dim))
                for ki in range(1,opts.dim+1):
                    W[:,ki-1] = wwhat(np.conj(omega), opts.beta, opts.gam, ki-1)
                
                windowq = np.dot(W ,u)
            
            elif (opts.motherwavelet == 'morse-c'):
                W = np.zeros(len(omega), opts.dim)
                for ki in range(1,opts.dim+1):
                    W[:,ki-1] = wwhat(np.conj(omega), opts.beta, opts.gam, ki-1)
            
                windowq = np.dot(W,uFix)
            
            elif (opts.motherwavelet == 'morse-a'):
                if(opts.k == 0):
                    windowq1 = wwhat(np.conj(omega),opts.beta,opts.gam,0)
                    windowq2 = wwhat(np.conj(omega),opts.beta,opts.gam,1)
                    windowq = ( windowq1 + windowq2 ) / sqrt(2)
                elif (opts.k == 1):
                    windowq1 = wwhat(np.conj(omega),opts.beta,opts.gam,0);
                    windowq2 = wwhat(np.conj(omega),opts.beta,opts.gam,1);
                    windowq = ( windowq1 - windowq2 ) / sqrt(2);
            
                windowq = np.conj(windowq.T)
            
            elif (opts.motherwavelet == 'Cinfc'):
                tmp0 = (omega-opts.CENTER)/opts.FWHM
                tmp1 = power(tmp0,2)-1
                windowq = np.exp(1/tmp1)
                windowq[omega >= (opts.CENTER+opts.FWHM)] = 0
                windowq[omega <= (opts.CENTER-opts.FWHM)] = 0
                windowq = np.array([windowq]).T

            elif (opts.motherwavelet == 'morlet'):
                windowq = 4*sqrt(pi)*np.exp(-4*power(omega-0.69*pi,2))-4.89098e-4*4*sqrt(pi)*np.exp(-4*power(omega,2))
                windowq = np.array([windowq]).T

            elif (opts.motherwavelet == 'gaussian'):
                psihat = lambda f: np.exp( -np.log(2)*power( 2*(f-opts.CENTER)/opts.FWHM,2) )
                windowq = psihat(omega)
                windowq = np.array([windowq]).T


            elif (opts.motherwavelet == 'meyer'):
                windowq = np.zeros(len(omega))
                int1 = np.logical_and(omega>=5/8*0.69*pi , omega<0.69*pi)
                int2 = np.logical_and(omega>=0.69*pi , omega<7/4*0.69*pi)
                meyeraux = lambda f: 35*power(f,4)-84*power(f,5)+70*power(f,6)-20*power(f,7)
                windowq[int1] = np.sin(pi/2*meyeraux((omega[int1]-5/8*0.69*pi)/(3/8*0.69*pi)))
                windowq[int2] = np.cos(pi/2*meyeraux((omega[int2]-0.69*pi)/(3/4*0.69*pi)))
                windowq = np.array([windowq]).T


            elif (opts.motherwavelet == 'BL3'):
                phihat = power(2*pi,-0.5)*power((np.sin(omega/4)/(omega/4)),4)
                phihat[0] = power(2*pi,-0.5)
                aux1 = 151/315 + 397/840*np.cos(omega/2) + 1/21*np.cos(omega) + 1/2520*np.cos(3*omega/2)
                phisharphat = phihat*power(aux1,-0.5)

                aux2 = 151/315 - 397/840*np.cos(omega/2) + 1/21*np.cos(omega) - 1/2520*np.cos(3*omega/2)
                aux3 = 151/315 + 397/840*np.cos(omega) + 1/21*np.cos(2*omega) + 1/2520*np.cos(3*omega)
                msharphat = power(np.sin(omega/4),4)*power(aux2,0.5)*power(aux3,-0.5)
                windowq = phisharphat*msharphat*np.exp(1j*omega/2)*(omega>=0)
                windowq = np.array([windowq]).T
                


            windowq = windowq/sqrt(qscale)
            what = windowq*np.array([xhat]).T
            w = np.fft.ifft(what.T)
            tfr[:,kscale-1] = w
            kscale+=1



        scale*=2




    """%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """%% calculate the constant for reconstruction
    %% TODO: calculate Rpsi for other mother wavelets"""
    xi = np.arange(0.05,10+1/10000,1/10000)
    
    if (opts.motherwavelet == 'gaussian'):   
        """%% Gaussian (not really wavelet)"""

        psihat = lambda f: np.exp( -np.log(2)*power( 2*(f-opts.CENTER)/opts.FWHM,2))
        windowq = psihat(xi)
        Rpsi = sum(windowq/xi)/10000


    elif (opts.motherwavelet == 'morlet'):
        windowq = 4*sqrt(pi)*np.exp(-4*power(xi-0.69*pi,2))-4.89098e-4*4*sqrt(pi)*np.exp(-4*power(xi,2))
        Rpsi = sum(windowq/xi)/10000

    elif (opts.motherwavelet == 'Cinfc'):
        tmp0 = (xi - opts.CENTER)/opts.FWHM
        tmp1 = power(tmp0,2)-1
        windowq = np.exp(1/tmp1)
        windowq[xi >= (opts.CENTER+opts.FWHM)] = 0
        windowq[xi <= (opts.CENTER-opts.FWHM)] = 0
        Rpsi = sum(windowq/xi)/10000
 
    else:
        """%% Normalization is not implemented for Other mother wavelets"""
        Rpsi = 1 
  
    tfr = tfr/Rpsi
    
    return tfr.T, tfrtic


