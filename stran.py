import numpy as np
import scipy as sp

def stran(h):
    """% Compute S-Transform without for loops

    %%% Coded by Kalyan S. Dash %%%
    %%% IIT Bhubaneswar, India %%%"""
    N = len(h) 
    """% h is a 1xN one-dimensional series"""

    nhaf=int(np.fix(N/2))
    odvn=1

    if (nhaf*2==N):
        odvn=0
    
    f = np.concatenate((np.arange(0,nhaf+1),np.arange(-nhaf+1-odvn,0)))/N
                      
    Hft=np.fft.fft(h)

    """%Compute all frequency domain Gaussians as one matrix"""

    invfk=np.conj(np.array([1/f[2:nhaf+2]])).T

    W=2*np.pi*np.matlib.repmat(f,nhaf,1)*np.matlib.repmat(invfk,1,N)
    
    G=np.exp((-np.power(W,2))/2) 
    
    """%Gaussian in freq domain"""
    
    """% End of frequency domain Gaussian computation"""

    """% Compute Toeplitz matrix with the shifted fft(h)"""

    HW=sp.linalg.toeplitz(np.conj(np.array([Hft[0:nhaf+1]])).T,Hft)

    """% Exclude the first row, corresponding to zero frequency"""

    HW=HW[1:nhaf+1,:]

    """% Compute Stockwell Transform"""

    ST=np.fft.ifft(HW*G,axis=1) 
    """%Compute voice"""
    
    """%Add the zero freq row"""
    st0=np.mean(h)*np.ones((1,N))
    
    ST=np.concatenate((st0,ST),axis = 0)
    
    return ST