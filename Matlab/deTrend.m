function [x, trend] = deTrend(t, z0, LowFreq) ;

if mod(size(t,1), 2)
	error('The signal should be of even length') ;
end

Hz = 1 ./ (t(2)-t(1)) ;
alpha = Hz/2/1000 ;

opts.motherwavelet = 'Cinfc' ;
opts.CENTER = 1 ;
opts.FWHM = 0.2 ;

[~, ~, tfrsq, ~] = sqCWTbase(t, z0, LowFreq, Hz/2, alpha, opts, 0, 0);

x = 2 * alpha * real( sum(tfrsq, 1) ) ;
trend = z0 - x' ;

x = x(:) ;
trend = trend(:) ;
