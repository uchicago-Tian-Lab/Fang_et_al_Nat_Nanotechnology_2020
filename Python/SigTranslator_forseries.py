# -*- coding: utf-8 -*-
import numpy as np #Numpy makes handling arrays much easier
import matplotlib.pyplot as plt #Matplotlib helps us to visualize and plot data
from skimage import io
from scipy.stats import pearsonr
from skimage.measure import label

#D_SERIES contains the numbers appended to the end of each video; SUB_INT is the time interval to analyze
D_SERIES = ['24','25','27','29','31','32','33','35','36','37','38','39','40','41']
SUB_INT = 100
mask = io.imread('Mask.tif')
labs,num_cells = label(mask, background=None, return_num=True, connectivity=1)
io.imsave('Labs.tif',labs.astype('float32')) 

d_corr = []
d_names = []
d_lbs = []
for d in D_SERIES:
    nm = "dFoF-50fps_"+d+".tif"
    lb = "Labs_"+str(d)+".tif"
    d_names.append(nm)
    d_lbs.append(lb)
    stack = io.imread(nm)   
    labs = io.imread(lb)
    num_cells = np.max(labs)
    signals = np.zeros([num_cells,np.shape(stack)[0]])
    for i in range(1,num_cells+1): #Our first loop!    
        cellmask = np.asarray(labs == i,stack.dtype)#Take only the region (cell) of the image indexed by i
        cellmask[cellmask == 0] = np.nan #To make the averaging work, turn all 0s into NaNs
        for j in range(np.shape(stack)[0]):
            myfuncell = cellmask*stack[j]#By default, multiplication in np is elementwise
            myfuncell[myfuncell == np.inf] = np.nan
            signals[i-1,j] = np.nanmean(myfuncell)
            
    #Show array containing summarized data
    g_name = "Signals_"+str(d)+".tif"
    plt.imshow(signals, cmap="gray")
    plt.axis('off')#This suppresses the pesky axis labels.
    io.imsave(g_name,signals.astype('float32'))
    plt.close()

    n_frames = np.shape(signals)[1]
    n_slices = n_frames//SUB_INT
    for i in range(n_slices):
        t_range = np.arange(i*100,(i+1)*100-1)
        sub_sigs = signals[:,t_range]
        corr,_ = pearsonr(sub_sigs[0],sub_sigs[1])
        d_corr.append(corr)
       
    # Use this for situations with several different clusters of oscillating cells    
    # sigmat = pd.DataFrame(signals.T)
    # corrmat = sigmat.corr()
    # c_name = "Corrmat_"+str(d)+".tif"
    # plt.imshow(corrmat)
    # io.imsave(c_name,corrmat.to_numpy().astype('float32'))   