#utility routines used by spec_combine.py for cleaning and processing ESI images

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.interpolate import interp1d
import spec_utils as u
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from PyAstronomy import pyasl as p
import argparse as ap
#----------------------------------------------------------------
def locate(array,value):
    idx = (np.abs(array-value)).argmin()
    if 0 < idx < len(array)-1:
        return idx
    else:
        return -1
#----------------------------------------------------------------
def maskim(im, thresh_hi=None, thresh_lo=None, col=None, row=None, replace=0., nn=1):
    if col is not None:
        col_arr = [col] if len(np.shape(col))==0 else col
        for col in col_arr: im.data[:,col] = np.median(np.median((im.data[:,col_arr[0]-nn:col_arr[0]],im.data[:,col_arr[-1]+1:col_arr[-1]+nn+1]),axis=0),axis=1) if replace=='nbr' else replace
    if row is not None:
        row_arr = [row] if len(np.shape(row))==0 else row
        for row in row_arr: im.data[row,:] = np.median(np.median((im.data[:,row_arr[0]-nn:row_arr[0]],im.data[:,row_arr[-1]+1:row_arr[-1]+nn+1]),axis=0),axis=1) if replace=='nbr' else replace
    if thresh_hi is not None or thresh_lo is not None:
        if thresh_hi is None: thresh_hi = np.min(im.data)-0.01
        if thresh_lo is None: thresh_lo = np.max(im.data)+0.01
        im.data[(im.data < thresh_lo) & (im.data > thresh_hi)] = replace #anything higher than thresh_hi and and anything lower than thresh_lo will be masked
    return im
#----------------------------------------------------------------
path = '/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_1x1/IRAF/'
filename = sys.argv[1] if len(sys.argv)>0 else 'std1_cos.fits'
if filename[-5:] != '.fits': filename += '.fits'
binx, biny = 1,1
thresh_hi = None
thresh_lo = None
col = np.arange(448,465+1,1)
row = None
replace_val = 'nbr'
nn = 3
hide = 1
suffix = '_clean'

fullname = path + filename
data,head = u.get_fitsspec(fullname)
im = ap.Namespace()
im.data = np.array(data)

if not hide: 
    fig = plt.figure(figsize=(8,8))
    plt.imshow(np.log10(im.data),cmap='gray',aspect='auto',vmin= 1.2, vmax=2.2)
    plt.title('Before')
    plt.show(block=False)

im = maskim(im, col=col, row=row, thresh_hi=thresh_hi, thresh_lo=thresh_lo, replace=replace_val, nn=nn)
u.write_fits(fullname[:-5]+suffix+'.fits', im.data, head)

if not hide: 
    fig = plt.figure(figsize=(8,8))
    plt.imshow(np.log10(im.data),cmap='gray',aspect='auto',vmin= 1.2, vmax=2.2)
    plt.title('After')
    plt.show(block=False)
