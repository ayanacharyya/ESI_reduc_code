#utility routine for combining and cleaning flux calibrated, dispersion corrected, echelle spectra of several orders
#by appropriately treating the overlapping regions

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.interpolate import interp1d
import readspec  as r
import spec_utils as u
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import time
start_time = time.time()
#----------------------------------------------------------------
global p, c, check, pause, inspect
show_disp_sol = 0
show_ord = 0
show_full = 1
clean = 0
check = 0
write = 1
inspect = 0
pause = 0.00001
do_wt_mean = 0
#----------------------------------------------------------------
ap_max = 8
sig_thresh = 5.
thresh = 7.
npix = 2
z = 1.420 #1.3294 # ##redshift is 1.3294 for S1723 and 1.420 for S2340
wrmin, wrmax = 4700./(1+z), 10200./(1+z) #default wavelength range 4700-10200 A in observed frame
#wrmin, wrmax = 4190, 4200#2405, 2415
path = '/Users/acharyya/Documents/esi_2016b/2016aug27_1x1/IRAF/'
wc_file = path + 'database/ec' + 'arcs.tr.ec'#'Arcs_sum.tr.ec' #
filename = path + 'obj1_cos_1dspeca_dispcor_fluxcal.fits'
fout = 'junk'#'s1723_arc_b.txt'
fout = '/Users/acharyya/Documents/esi_2016b/2016aug27_1x1/IRAF/reduced/' + fout
#----------------------------------------------------------------
s,h = r.get_fitsspec(filename)
spec = pd.DataFrame(columns=['obswave','fnu','fnu_u']) #pandas dataframe to hold the spectrum
wlist = u.get_dispsol(wc_file, ap_max, show=show_disp_sol)
wave = np.arange(wlist[0][1],wlist[-1][-2],0.2)
#---------------------------------------------------------------- 
flux,error=[],[]
for w in wave:
    f,e=[],[]
    for ap in range(ap_max):
        id = u.locate(wlist[ap],w)
        if id > 0 :
            f.append(s[0,ap,id])
            e.append(s[3,ap,id])
    if do_wt_mean:
        f,e = u.wt_mean(f,e)
    else:
        f,e = u.get_least_error(f,e)
    flux.append("{0:.2e}".format(f))
    error.append("{0:.2e}".format(e))
#---------------------------------------------------------------- 
spec['obswave'] = pd.Series(["{0:.2f}".format(w) for w in wave])
spec['fnu'] = pd.Series(flux)
spec['fnu_u'] = pd.Series(error)
spec.obswave = spec.obswave.astype(np.float64)
spec.fnu = spec.fnu.astype(np.float64)
spec.fnu_u = spec.fnu_u.astype(np.float64)
spec['restwave'] = spec['obswave']/(1.+z)
spec = spec[spec['restwave'].between(wrmin,wrmax)]
#---------------------------------------------------------------- 
if show_ord:       
    fig = plt.figure(figsize=(18,4))
    fig.subplots_adjust(hspace=0.7, top=0.94, bottom=0.1, left=0.06, right=0.98)
    cmap = plt.get_cmap('rainbow')
    colors = [cmap(i) for i in np.linspace(0, 1, ap_max)]
    for ap in range(ap_max):
        plt.plot(wlist[ap],s[0,ap,:],c=colors[ap])
        plt.plot(wlist[ap],s[3,ap,:])
    plt.ylim(-4e-17,20e-17)
    #plt.ylim(-50,60000) #
    plt.ylabel('fnu(ergs/cm^2/s/Hz)')
    plt.xlabel('Observed wavelength (A)')
    plt.show(block=False)
#---------------------------------------------------------------- 
if clean:
    if check: 
        r.plot_spec(spec, 'masking', check=True, sig_thresh=sig_thresh)
        p= plt.axvline(spec['restwave'].values[0], c='r')
    spec['badmask']  = False
    u.cleanspec(spec)
#---------------------------------------------------------------- 
if write:
    u.writespec_txt(spec, fout, z=z, filename=filename)
#---------------------------------------------------------------- 
if write and show_full:
    spec = r.get_txtspec(fout, wrmin=wrmin, wrmax=wrmax, sig_thresh=sig_thresh, check=check, clean=False)
    if clean:
        spec = r.get_txtspec(fout, wrmin=wrmin, wrmax=wrmax, sig_thresh=sig_thresh, check=check)
#---------------------------------------------------------------- 
print('Done in %s minutes' % ((time.time() - start_time)/60))

