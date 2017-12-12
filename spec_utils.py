#utility routines used by spec_combine.py for combining and cleaning flux calibrated, dispersion corrected, echelle spectra

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from PyAstronomy import pyasl as p
#----------------------------------------------------------------
def get_fitsspec(filename):
    d = fits.open(filename)
    #print d.info()
    h = d[-1].header
    s = d[-1].data
    return s,h
#-------------------function for saving array as fits file--------------
def write_fits(filename, data, head):
    hdu = fits.PrimaryHDU(data, header=head)
    hdulist = fits.HDUList([hdu])
    if filename[-5:] != '.fits':
        filename += '.fits'
    hdulist.writeto(filename, clobber=True)
    print 'Written file', filename
#----------------------------------------------------------------
def mask_img(filename, hi_thresh=None, lo_thresh=None):
    s,h = get_fitsspec(filename)
    if hi_thresh is None and lo_thresh is None:
        print 'mask_img ERROR: Need either hi_thresh or lo_thresh or both to mask image.'
        return
    if hi_thresh is not None:
        s[s > hi_thresh] = 0
    if lo_thresh is not None:
        s[s < lo_thresh] = 0
    write_fits(filename,s,h)
#----------------------------------------------------------------
def get_txtspec(filename, wrmin = None, wrmax = None, show=True, check = False, sig_thresh=5., clean=True):
    spec =  pd.read_table(filename, delim_whitespace=True, comment="#", header=0, dtype=np.float64)
    if wrmin is None:
        wrmin = spec['restwave'].values[0]
    if wrmax is None:
        wrmax = spec['restwave'].values[-1]
    spec = spec[spec['restwave'].between(wrmin,wrmax)]
    if clean and 'badmask' in spec:
        spec.badmask = spec.badmask.astype(bool)
        spec = spec[~spec.badmask]
    if show:
        plot_spec(spec, fout=filename+' clean= '+str(clean), check = check, sig_thresh=sig_thresh)
    return spec
#----------------------------------------------------------------
def plot_spec(spec, fout=None, check = False, sig_thresh=5.):
    fig = plt.figure(figsize=(18,4))
    fig.subplots_adjust(hspace=0.7, top=0.94, bottom=0.1, left=0.06, right=0.98)
    wave_axis = 'restwave'
    flux_axis = 'flam'
    plt.plot(spec[wave_axis],spec[flux_axis])
    plt.plot(spec[wave_axis],spec[flux_axis+'_u'])
    if check: 
        plt.axhline(np.median(spec[flux_axis].values)*sig_thresh, linestyle='--', c='r') #
        plt.axhline(-np.median(spec[flux_axis].values)*sig_thresh, linestyle='--', c='r') #
    #plt.ylim(np.min(spec[flux_axis]), np.max(spec[flux_axis]))
    plt.ylim(-4e-17,9e-17)
    #plt.ylim(-50,60000) #
    plt.xlim(spec[wave_axis].values[0],spec[wave_axis].values[-1])
    plt.ylabel(flux_axis) #'(ergs/cm^2/s/Hz)' or (ergs/cm^2/s/A)
    plt.xlabel(wave_axis+' (A)')
    plt.title(fout)
    plt.show(block=False)
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
def domask(spec, jj, check=False, case=''):
    spec.badmask.loc[jj] = True
    if len(jj)== 0: jj = [jj]
    if check: 
        print 'bad due to case '+str(case)
        for j in jj: plt.axvline(spec['restwave'].values[j], c='k')
'''
#----------NOT USED ANYMORE: REDUNDANT----------------------------
def wt_mean(data, error, sky=None):
    if len(data) != len(error):
        print 'Error during weighted mean calculation: data and error should have same length'
        return -1
    data = np.array(data)
    error = np.array(error)
    if sky is not None:
        sky = np.array(sky)
    inv_var = 1./error**2
    wd = np.inner(data,inv_var)/np.sum(inv_var)
    we = np.inner(error,inv_var)/np.sum(inv_var)
    if sky is not None:
        ws = np.inner(sky,inv_var)/np.sum(inv_var)
        return wd, we, ws
    else:
        return wd, we

#---------NOT USED ANYMORE BECAUSE: UNSCIENTIIFIC--------------------
def get_least_error(data,error, sky=None):
    if len(data) != len(error):
        print 'Error during weighted mean calculation: data and error should have same length'
        return -1
    id = np.where(error == np.min(error))[0][0]
    if sky is None:
        return data[id], error[id]
    else:
        return data[id], error[id], sky[id]
'''
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
    head = '1D spectrum for object file '+filename+'\n\
    made using python spec_combine.py\n\
    by Ayan, Oct 2016\n\
    redshift z = '+str(z)+'\n\
    Columns are:\n\
    flam: in units of ergs/s/A/cm^2\n\
    flam_u: uncertainties in above\n\
    obswave: observed-frame vacuum wavelength in Angstrom\n\
    restwave: obswave converted to restframe (vacuum)\n\
    badmask: (if present) "True" indicates that the pixel is bad and hence masked. Note: generation of this mask is automated.\n\
    NOTE: For any known bad column region flam is set to 0. and flam_u to 99.\n\
    '
    np.savetxt(fout, [], header=head, comments='#')
    spec.to_csv(fout, sep='\t',mode ='a', index=None)
    print 'Written dataframe to file ', fout
#------------------------------------------------------------------
def get_dispsol(wc_file, ap_max, show=False, col='b', is_wave_air=False, label=None):
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
            if ap > 1: label = None
            plt.plot(ypix,w, c=col, label=label)
            plt.axhline(w[-1], linestyle='dotted',c='k') #
        wlist.append(w)
    if show:
        plt.xlabel('y-pixel')
        plt.ylabel('Observed wavelength(A)')
        plt.title('Dispersion solution')
        plt.show(block=False)
    return np.array(wlist)
#---------------------------------------------------------------- 