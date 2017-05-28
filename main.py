"""
%% in this simulation example, how to run ConceFT is demonstrated with the conceFT for the synchrosqueezed STFT. The same flow could be appied to conceFT for the synchrosqueezed CWT or others. Please see the code for details.

	%% the sampling rate for the simulated signal
"""
import numpy as np
import math
import statsmodels.api as sm
from functions import *    
Hz = 100 
""" the sampling time of the simulated signal"""
time = np.arange(0,16,1/Hz)
"""the number of sampling points of the simulated signal"""

N = len(time) 

"""fix the random seed for the reproducibility issue"""
np.random.seed(1)


"""the amplitude modulation of the simulated signal"""
"""simulate 2 oscillatory components with dynamics"""

am1 = sm.nonparametric.lowess(np.cumsum(np.random.randn(N))/ Hz,time,  frac=1/8)
am1 = 2 + am1[:,1]/ np.amax(np.absolute(am1[:,1])) 
am2 = sm.nonparametric.lowess(np.cumsum(np.random.randn(N))/ Hz,time,  frac=1/8)
am2 = 2 + am2[:,1]/ np.amax(np.absolute(am2[:,1])) 
am1[0:500] = 0 
am2[-601:] = 0 

"""time,am1,am2 are 1-d array here"""

"""the instantaneous frequency of the simulated signal"""
if1 = sm.nonparametric.lowess(np.cumsum(np.random.randn(N))/ Hz,time,  frac=1/4)
if1 = 10 + 6*if1[:,1]/ np.amax(np.absolute(if1[:,1])) 
if2 = sm.nonparametric.lowess(np.cumsum(np.random.randn(N))/ Hz,time,  frac=3/16)
if2 = math.pi + 3*if2[:,1]/ np.amax(np.absolute(if2[:,1])) 
phi1 = np.cumsum(if1) / Hz 
phi2 = np.cumsum(if2) / Hz 

"""the simulated signal."""
s1 = am1 * np.cos(2*math.pi*phi1) 
s2 = am2 * np.cos(2*math.pi*phi2)  
clean = s1 + s2 

if1[0:500] = np.NaN
if2[-601:] = np.NaN
am1[0:500] = np.NaN
am2[-601:] = np.NaN


"""add noise (Gaussian white noise)"""
sigma = 1 
"""sqrt( var(clean)*10.^( -snrdb /10 ) )"""
noise = np.random.standard_t(4, size=N) 
noise = sigma * noise  
np.var(noise)
snrdb = 20 * np.log10(np.std(clean)/np.std(noise)) 
print('snrdb = ',snrdb,'\n') 

"""simulated observed time series"""
xm = clean + noise 
"""setup parameters for the SST or ConceF"""

"""number of chosen orthonormal windows for ConceFT """
NoWindowsInConceFT = 2 
"""number of random linear combinations of chosen windows"""
NoConceFT = 20 
"""
the window length. Ideally, it should be chosen so that
roughly 7-10 oscillations (ignore the multiples) are 
included in the window.
"""
WindowLength = 377 
"""
this is the bandwith of the chosen window. See hermf.py
in the attached code for details.
"""
WindowBandwidth = 10 
SamplingRate = Hz 
"""
Setup the frequency range for the analysis
The allowed range is 0-0.5
This part might be tricky. This 0-0.5 limitation is 
setup under the assumption that the sampling rate is 1Hz
After the analysis, the sampling rate should be adjusted
so that the result is associated with the original sampling rate. 
In this example, the true range is [0, 0.5]*SamplingRate
"""
HighFrequencyLimit = 0.5 
LowFrequencyLimit = 0 
"""the frequency axis resolution in the final time-frequency representation"""
FrequencyAxisResolution = 0.001 

"""
	%% call the main code, which is the ConceFT based on 
	%% synchrosqueezed short time Fourier transform (STFT)
	%% Output:
	%% tfr: STFT result
	%% tfrtic: frequency axis tic for the STFT
	%% tfrsq: synchrosqueezed STFT (it is equivalent to running ConceFT only one time)
	%% ConceFT: ConceFT of synchrosqueezed STFT.
	%% tfrsqtic: frequency axis tic for the tfrsq and ConceFT
"""
tfr, tfrtic, tfrsq, ConceFT, tfrsqtic = ConceFT_CWT(time,xm, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, 1, opts, 0,0) 

"""
	%% plot the time frequency representation determined by
	%% ConceFT. .995 is the quantile truncation set to avoid 
	%% possible outliers in the final analysis result.
	%% see the code for details.
"""
import matplotlib.pyplot as plt
fig = plt.figure()
ax1 = plt.subplot(1, 2, 1)
imageSQ(ax1, time, tfrsqtic*SamplingRate, np.abs(tfrsq), 99.5,'Greys')
plt.title('SST')


ax2 = plt.subplot(1, 2, 2)
imageSQ(ax2, time, tfrsqtic*SamplingRate, np.abs(ConceFT), 99.5, 'Greys')
plt.title('ConceFT')
fig.savefig('plot.png')