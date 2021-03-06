function [tfr, tfrtic, tfrrs, ConceFT, tfrrstic] = ConceFT_rsSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT) ;


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

fprintf(['Run ordinary RS\n']) ;
[tfr, tfrtic, tfrrs, tfrrstic] = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, h(1,:)', Dh(1,:)', 1);


%=======================================
ConceFT = [] ;

if MT > 1

	disp('Complex Sphere') ;
	%% Conceft
    ConceFT = zeros(size(tfrrs)) ;

	fprintf(['ConceFT-RS total: ',num2str(MT),'; now:     ']) ;
    for ii = 1: MT
		fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',ii) ; fprintf([tmp]) ;
		rv = randn(1, dim) + sqrt(-1)*randn(1, dim) ; rv = rv ./ norm(rv) ;
		rh = rv * h ; 
		rDh = rv * Dh ;

		[~, ~, tfrrs, tfrrstic] = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, rh', rDh', 1);

	 	ConceFT = ConceFT + tfrrs ;
    end

    ConceFT = ConceFT ./ MT ;
	fprintf('\n') ;

end

end
