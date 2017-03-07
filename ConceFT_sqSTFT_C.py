import numpy as np
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
