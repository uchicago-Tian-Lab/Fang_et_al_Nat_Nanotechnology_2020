# -*- coding: utf-8 -*-
import numpy as np
import numpy.fft as fft
from skimage import io

MAX_FREQ = 1/2 #The maximum frequency we could observe is every other frame
CELL_ROI = 0 #The cell being examined (starting from index 0, i.e., pixel value - 1)
T_INT = 0.105 #The time interval between frames (in seconds)
SIG_NAME = "Signals"
SUB_INT = 100 #The number of time frames considered in each spectral analysis

def SpectralAnalyzer(signal,t_int):
    # ###!!!Only use this if the power spectrum is dominated by low-frequency noise!
    # It fits a line to the periodic signal and subtracts it, so long-term trends are suppressed
    # dummyx = np.arange(0,np.shape(signal)[1])
    # lfit = np.polyfit(dummyx,signal,1)
    # lval = np.polyval(lfit,dummyx)
    # signal = signal - lval
    sigfft = fft.fft(signal)
    n = np.shape(sigfft)[0]
    fullpow = np.abs(sigfft[1:n//2])**2
    freq = np.arange(1,n//2)/(n/2)*MAX_FREQ/t_int
    return fullpow, freq

#Initialize structures for the remainder of the analysis
D_SERIES = ['24','25','27','29','31','32','33','35','36','37','38','39','40','41']
d_names = []
d_freqs = []
d_pows = []
maxfreqs = []
maxpows = []
for d in D_SERIES:
    f_name = SIG_NAME+"_"+d+".tif"
    img = io.imread(f_name)
    n_frames = np.shape(img)[1]
    n_slices = n_frames//SUB_INT
    for j in range(n_slices):
        d_names.append(d+"_"+str(j))
        t_range = np.arange(j*100,(j+1)*100-1)
        sig = img[CELL_ROI,t_range] #Only looking at the signal from the interesting cell
        power,freq = SpectralAnalyzer(sig,T_INT)
        d_freqs.append(freq)
        d_pows.append(power)
        maxfreqs.append(freq[np.where(power == np.max(power))[0].item()])
        maxpows.append(np.max(power))