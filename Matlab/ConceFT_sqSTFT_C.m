function [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Second, Smooth, Hemi) ;

%
% Usage: 
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFT(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp, MT)
%
% MT = 1: ordinary SST; MT > 1: ConceFT
% alpha: resolution in the frequency axis
% WinLen, dim, supp: for hermf.m
%
% Example:
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFT([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6, 10);


N = length(x) ;

%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       	%% generate the window for short time Fourier transform (STFT)
[h, Dh, ~] = hermf(WinLen, dim, supp) ;

%=======================================

fprintf(['Run ordinary STFT-SST (Smooth = ',num2str(Smooth),', Hemi = ',num2str(Hemi),')\n']) ;
[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', Smooth, Hemi) ;

%=======================================

ConceFT = tfrsq ;


if MT > 1

	disp('Complex sphere') ;
		%% Conceft

	fprintf(['STFT-ConceFT total (Smooth = ',num2str(Smooth),', Hemi = ',num2str(Hemi),'): ',num2str(MT),'; now:     ']) ;

    for ii = 1: MT
		fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',ii) ; fprintf([tmp]) ;
		rv = randn(1, dim) + sqrt(-1)*randn(1, dim) ; rv = rv ./ norm(rv) ;
		rh = rv * h ; 
		rDh = rv * Dh ;

		if ~Second
			[~, ~, tfrsqX, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, rh', rDh', Smooth, Hemi);
		else
			[~, ~, ~, tfrsqX, tfrsqtic] = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, rh', rDh', dwindow(rDh'), 0);
		end

	 	ConceFT = ConceFT + tfrsqX ;
    end

    ConceFT = ConceFT ./ (MT+1) ;
	fprintf('\n') ;

end
