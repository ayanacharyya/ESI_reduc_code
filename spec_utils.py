#utility routines used by spec_combine.py for combining and cleaning flux calibrated, dispersion corrected, echelle spectra

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.interpolate import interp1d
import readspec  as r
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from PyAstronomy import pyasl as p
#----------------------------------------------------------------
def locate(array,value):
    idx = (np.abs(array-value)).argmin()
    if 0 < idx < len(array)-1:
        return idx
    else:
        return -1
#----------------------------------------------------------------
def cleanspec_old(sp, thresh_hi=7e-17, thresh_lo=-4e-17):
    print 'Cleaning spectrum..'
    sp.badmask.loc[sp['fnu'].gt(thresh_hi)] = True
    sp.badmask.loc[sp['fnu'].lt(thresh_lo)] = True
    return 0
#----------------------------------------------------------------
def domask(spec, jj, check=False):
    spec.badmask.loc[spec.index[jj]] = True
    if check: 
        print 'bad'
        plt.axvline(spec['restwave'].values[jj], c='k')
#----------------------------------------------------------------
def wt_mean(data, error):
    if len(data) != len(error):
        print 'Error during weighted mean calculation: data and error should have same length'
        return -1
    data = np.array(data)
    error = np.array(error)
    inv_var = 1./error**2
    wd = (data*inv_var)/np.sum(inv_var)
    we = (error*inv_var)/np.sum(inv_var)
    return wd[0], we[0] 
#----------------------------------------------------------------
def get_least_error(data,error):
    if len(data) != len(error):
        print 'Error during weighted mean calculation: data and error should have same length'
        return -1
    id = np.where(error == np.min(error))[0][0]
    return data[id], error[id]
#----------------------------------------------------------------        
def wc(x,o,xpow,opow,xmin,xmax,omin,omax,c):
    nx=(2*x-(xmax+xmin))/(xmax-xmin)
    no=(2*o-(omax+omin))/(omax-omin)
    zx=[1,nx]
    zo=[1,no]
    for i in range(1,xpow-1):
        zx.append(((2*i+1) * nx * zx[i] - i * zx[i-1]) / (i+1)) 
    for i in range(1,opow-1):
        zo.append(((2*i+1) * no * zo[i] - i * zo[i-1]) / (i+1)) 
    s = np.sum(np.multiply(c,np.outer(zo,zx)))
    return s
#---------------------------------------------------------------- 
def writespec_txt(spec, fout, z='', filename=''):
    head = '#1D spectrum for object file '+filename+'\n\
    #made using python spec_combine.py\n\
    #by Ayan, Oct 2016\n\
    redshift z = '+str(z)+'\n'
    np.savetxt(fout, [], header=head, comments='#')
    spec.to_csv(fout, sep='\t',mode ='a', index=None)
    print 'Written dataframe to file ', fout
#------------------------------------------------------------------
def get_dispsol(wc_file, ap_max, show=False, col='b', is_wave_air=False):
    lines = open(wc_file,'r').readlines()
    shift, offset, slope = 0, 0, 1 #default values in case absent in ec* file
    for i,l in enumerate(lines):
        if '\tslope\t' in l:
            slope = float(l.split()[1])
        if '\toffset\t' in l:
            offset = float(l.split()[1])
        if '\tshift\t' in l:
            shift = float(l.split()[1])
        if '\tcoefficients\t' in l:
            start = i
    xpow = int(float(lines[start+2].strip()))
    opow = int(float(lines[start+3].strip()))
    xmin, xmax = float(lines[start+5].strip()), float(lines[start+6].strip())
    omin, omax = float(lines[start+7].strip()), float(lines[start+8].strip())
    c=np.reshape([float(lines[start+9+i].strip()) for i in range(opow*xpow)],(opow,xpow))
    ypix = np.arange(int(xmin),int(xmax)+1)
    wlist=[]
    for ap in range(1,ap_max+1):
        w=[]
        o = ap*slope + offset
        for x in ypix:
            w.append((wc(x,o,xpow,opow,xmin,xmax,omin,omax,c)+shift)/o)
        if is_wave_air:
            w = p.airtovac2(w, mode='ciddor') #to convert to vacuum wavelengths if initially it is air wavelength
        if show:
            plt.plot(ypix,w, c=col)
            plt.axhline(w[-1], linestyle='dotted',c='k') #
        wlist.append(w)
    if show:
        plt.xlabel('y-pixel')
        plt.ylabel('Observed wavelength(A)')
        plt.title('Dispersion solution')
        plt.show(block=False)
    return np.array(wlist)
#---------------------------------------------------------------- 