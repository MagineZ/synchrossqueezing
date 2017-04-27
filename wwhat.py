import numpy as np

def wwhat(f,beta,gam,k):
    """% This function calculates morse complex wavelets (in frequency domain)
    % It is very similar to morsefreqs"""
    from scipy.special import gamma
    from math import sqrt, pi
    from numpy import power
    c=-1+((2*beta+1)/gam)
    con=sqrt(gam*power(2,((2.*beta+1)/gam))*gamma(k+1))
    con2=sqrt(pi*gamma(k+((2.*beta+1)/gam)))
    con=sqrt(pi/2)*(con/con2)
    arr=2*(power(2*pi*abs(f),gam))
    m=con*power(2*pi*abs(f),beta)*np.exp(-power(2*pi*abs(f),gam))*laggen(arr,k,c)
    m=sqrt(4*pi)*m
    
    for j in range(0,len(f)):
        if (f[j]<0):m[:,j]=0
    
    return m

def laggen(x,k,c):
    """
    % compute generalized Laguerre poly L_k^c(x)
    %
    %
    """
    from scipy.special import gammaln
    from numpy import power
    l=len(x)
    sn=-1;
    s=np.zeros((1,l))
    for m in range(0,k+1):
        sn = -sn
        ga = gammaln(k+c+1)-gammaln(k-m+1)-gammaln(c+m+1)-gammaln(m+1)
        ga = np.exp(ga)
        s=s+sn*ga*power(x,m)
    return s
