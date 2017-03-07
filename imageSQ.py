import numpy as np
import matplotlib.pyplot as plt


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
    
    
