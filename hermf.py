"""
Created on Sat Feb 11 16:35:23 2017

@author: Pei-Chun Su
"""


"""
import all needed packages
"""
import numpy as np
import math


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