def ConceFT_CWT(t, x, lowfreq, highfreq, alpha, MT, opts, Smooth, Hemi):
     
    """%% ordinary SST"""
    print('Run ordinary CWT-SST (Smooth = ',Smooth,', Hemi = ',Hemi,')\n')
    tfr, tfrtic, tfrsq, tfrsqtic = sqCWTbase(t, x, lowfreq, highfreq, alpha, opts, Smooth, 0)

    """%===========================
    %% get the ConceFT"""	

    ConceFT = tfrsq


    if MT > 1:
        print('CWT-ConceFT total (Smooth = ',Smooth,', Hemi = ',Hemi,'): ',MT,'; now:     ')
        
        for ii in range(1,MT+1):
            print('\b\b\b\b') 
            print('%4d' % ii)
            _, _, tfrsqX, tfrsqtic = sqCWTbase(t, x, lowfreq, highfreq, alpha, opts, Smooth, Hemi)
            
            ConceFT = ConceFT + tfrsqX 
	
        print('\n')

        ConceFT = ConceFT/(MT+1)
        
    return tfr, tfrtic, tfrsq, ConceFT, tfrsqtic


