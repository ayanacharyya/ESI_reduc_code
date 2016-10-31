#utility routine for reading in 1d spectra created using IRAF
#and returning the header and data

from astropy.io import fits
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#----------------------------------------------------------------
def get_fitsspec(filename):
    d = fits.open(filename)
    print d.info()
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
    plt.plot(spec[wave_axis],spec['fnu'])
    plt.plot(spec[wave_axis],spec['fnu_u'])
    if check: 
        plt.axhline(np.median(spec['fnu'].values)*sig_thresh, linestyle='--', c='r') #
        plt.axhline(-np.median(spec['fnu'].values)*sig_thresh, linestyle='--', c='r') #
    #plt.ylim(np.min(spec['fnu']), np.max(spec['fnu']))
    plt.ylim(-4e-17,9e-17)
    #plt.ylim(-50,60000) #
    plt.xlim(spec[wave_axis].values[0],spec[wave_axis].values[-1])
    plt.ylabel('fnu(ergs/cm^2/s/Hz)')
    plt.xlabel(wave_axis+' (A)')
    plt.title(fout)
    plt.show(block=False)
#---------------------------------------------------------------- 
if __name__ == 'main':
    filename = sys.argv[1] if len(sys.argv)>1 else '/Users/acharyya/Documents/esi_2016b/2016aug27_2x1/IRAF/obj1_sc_cos2.fits'
    s,h = get_fitsspec(filename)
    #sp = get_txtspec(filename)