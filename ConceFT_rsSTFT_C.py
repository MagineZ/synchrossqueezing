import numpy as np
def ConceFT_rsSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT):
    
    """
    %
    % Usage: 
    % 	tfrsq, ConceFT, tfrsqtic = sqSTFT(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp, MT)
    %
    % MT = 1: ordinary SST; MT > 1: ConceFT
    % alpha: resolution in the frequency axis
    % WinLen, dim, supp: for hermf.m
    %
    % Example:
    % 	tfrsq, ConceFT, tfrsqtic = sqSTFT([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6, 10);
    """
    

    """%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
    """generate the window for short time Fourier transform (STFT)"""
    
    h, Dh, _ = hermf(WinLen, dim, supp) 
    """%======================================="""
    print('Run ordinary RS\n')
    
    tfr, tfrtic, tfrrs, tfrrstic = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, np.conj(np.array([h[0,:]])),  np.conj(np.array([Dh[0,:]])), 1)

    """%======================================="""
    ConceFT = tfrrs

    if MT > 1:
        print('Complex sphere\n') 
        """%% Conceft"""
        ConceFT = np.zeros(tfrrs.shape) 
        print('RS-ConceFT total: ',MT,'; now:     ')
        
        for ii in range(1, MT+1):
            print('\b\b\b\b') 
            print('%4d' % ii) 
            rv = np.random.randn(1, dim) + np.sqrt(-1+0j)*np.random.randn(1, dim)
            rv = rv/np.linalg.norm(rv)
            rh = np.dot(rv , h)  
            rDh = np.dot(rv , Dh)
            _, _, tfrrs, tfrrstic = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, np.conj(rh), np.conj(rDh), 1)
            
            ConceFT = ConceFT + tfrrs
        
        ConceFT = ConceFT/MT 
        print('\n') 
    
    return tfr, tfrtic, tfrrs, ConceFT, tfrrstic
