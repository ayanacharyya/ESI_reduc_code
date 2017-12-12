#utility routine for combining and cleaning flux calibrated, dispersion corrected, echelle spectra of several orders
#by appropriately treating the overlapping regions

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.interpolate import interp1d
import spec_utils as u
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import astropy.convolution as con
import time
start_time = time.time()
#----------------------------------------------------------------
def diagnose(k, n):
    global p, pause, inspect
    p.remove()
    p= plt.axvline(spec['restwave'].values[k+n], c='r')
    print 'checking restwave['+str(k+n)+'] =', spec['restwave'].values[k+n], 'out of', spec['restwave'].values[-1]
    if inspect:
        plt.pause(pause)
#----------------------------------------------------------------
def rec(spec, thresh, k, n, nlim):
    global p,c,check,protected
    c += 1
    if check: diagnose(k,n)
    if n == nlim:
        if np.abs(spec['flam'].values[k+n]) > thresh:
            return 0
        else:
            return 1            
    else:
        if np.abs(spec['flam'].values[k+n]) > thresh:
            if check: print 'still suspect'
            if rec(spec, thresh, k, n+1, nlim):
                u.domask(spec, k+n, check=check, case=2)
                return 1
            else:
                if check: print 'innocent'
                protected = 1
                return 0 
        else:
            return 1
#----------------------------------------------------------------
def cleanspec(spec, thresh=7., sig_thresh=5., npix=2):
    global p, c, check, protected
    print 'Cleaning spectrum carefully. This might take a few minutes..'
    fm = np.median(spec['flam'].values)
    em = np.median(spec['flam_u'].values)
    ind,j, protected = [0]*3
    u.domask(spec, spec[spec.flam_u > sig_thresh*em].index, check=check, case=1) #checking if error is too much then get rid of it
    '''
    while j in range(len(spec)-npix):
        c = 0
        if check: diagnose(j,0)
        if not protected:
            if np.abs(spec['flam'].values[j]) > sig_thresh*fm: #checking if too high flam spanning <npix pixels then get rid of it
                if check: print 'suspect'
                if rec(spec, sig_thresh*fm, j, 1, npix):
                    u.domask(spec, j, check=check, case=2)
        elif np.abs(spec['flam'].values[j]) <= sig_thresh*fm:
            protected = 0        
        j += c+1
    
    dummy = spec['flam'].values[ind]
    for j in range(0,len(spec)):
        if check: diagnose(j,0)
        if np.abs(spec['flam'].values[j]-dummy) > thresh*fm*(j - ind):#checking if slope of rise (sharpness) is beyond certain thresh, then get rid of it
            u.domask(spec, j, check=check, case=3)
        else:
            ind = j
            dummy = spec['flam'].values[ind]          
    '''
    return 0
#---------For boxcar smoothing-------------------------------------------------------
def smooth(array, boxcar=501):
    return con.convolve(np.array(array), np.ones((boxcar,))/boxcar, boundary='fill', fill_value=np.nan)
#----------------------------------------------------------------
global p, c, check, pause, inspect
show_disp_sol = 0
show_ord = 0
show_ord_all = 0
show_before_dispsol = 0
order_to_show = int(sys.argv[1]) if len(sys.argv)>1 else 2
show_full = 0
show_full_all = 1
clean = 0
mask_bad_col = 1
check = 0
write = 0
inspect = 0
pause = 0.00001
is_wave_air = 1
#----------------------------------------------------------------
ap_max = 10
sig_thresh = 5.
thresh = 7.
npix = 2
biny = 1
nypix = 4096/biny
#wrmin, wrmax = None,None #
#----------------------------------------------------------------
z_arr = np.sort([1.32952, 1.420]*8) #1.420 # ##redshift is 1.3294 for S1723 and 1.420 for S2340
z1 = 3.625 #for J1050
z2 = 3.487 #for J1458
z_arr = np.append(z_arr, z1)
z_arr = np.append(z_arr, z2)
z_arr = np.append(z_arr, z1)
z_arr = np.append(z_arr, z1)
z_arr = np.append(z_arr, z2)
path = ['/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_1x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_1x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug28_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug28_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug28_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug28_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug28_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug28_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug28_2x1/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016a/2016apr02/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016a/2016apr02/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016a/2016apr03/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016a/2016apr03/IRAF/',\
'/Users/acharyya/Work/astro/obs_data/esi_2016a/2016apr03/IRAF/']

fn = ['obj1_cos_clean_1dspeca_dispcor_fluxcal.fits',\
'obj1_cos_clean_1dspeca_dispcor_fluxcal.fits',\
'obj1_sc_cos_1dspecb_dispcor_fluxcal.fits',\
'obj2_cos_clean_1dspeca_dispcor_fluxcal.fits',\
'obj1b_cos_clean_1dspeca_dispcor_fluxcal.fits',\
'obj11_cos_clean_1dspeca_dispcor_fluxcal.fits',\
'obj11_cos_1dspecb_dispcor_fluxcal.fits',\
'obj12_cos_clean_1dspeca_dispcor_fluxcal.fits',\
'obj345_sc_cos_1dspeca_dispcor_fluxcal.fits',\
'obj345_sc_cos_1dspecb_dispcor_fluxcal.fits',\
'obj8_sc_cos_1dspeca_dispcor_fluxcal.fits',\
'obj8_sc_cos_1dspecb_dispcor_fluxcal.fits',\
'obj212223_cos_1dspeca_dispcor_fluxcal.fits',\
'obj212223_cos_1dspecb_dispcor_fluxcal.fits',\
'obj24_cos_1dspeca_dispcor_fluxcal.fits',\
'obj24_cos_1dspecb_dispcor_fluxcal.fits',\
'obj1_med_cos_1dspeca_dispcor_fluxcal.fits',\
'obj2_med_cos_1dspeca_dispcor_fluxcal.fits',\
'obj1_med_cos_1dspec2_b_dispcor_fluxcal.fits',\
'obj1_med_cos_1dspec3_b_dispcor_fluxcal.fits',\
'obj2_med_cos_1dspecb_dispcor_fluxcal.fits']

fo = ['s1723_arc_a.txt',\
's1723_center_a.txt',\
's1723_counter_a.txt',\
's1723_side_a.txt',\
's1723_arc_b.txt',\
's1723_center_b.txt',\
's1723_counter_b.txt',\
's1723_side_b.txt',\
's2340_a1_a.txt',\
's2340_a2_a.txt',\
's2340_a4_a.txt',\
's2340_a3_a.txt',\
's2340_a1_b.txt',\
's2340_a2_b.txt',\
's2340_a4_b.txt',\
's2340_a3_b.txt',\
'J1050_arc_a.txt',\
'J1458_arc_a.txt',\
'J1050_2_b.txt',\
'J1050_3_b.txt',\
'J1458_arc_b.txt']

mask_lim_arr = [[[4345, 4370], [4500, 4600]],\
[[4345, 4370], [4464, 4528]],\
[],\
[[4345, 4370], [4484, 4572]],\
[[4345, 4370], [4500, 4600]],\
[[4345, 4370], [4500, 4572]],\
[],\
[[4345, 4370], [4500, 4572]],\
]
outpath = '/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_2x1/IRAF/reduced/'
#fo = ['junk']*16 #
#----------------------------------------------------------------
print 'galaxy_name\t\tRA\t\tDEC\t\tDATE-OBS\t\tUT_START\t\tEXPTIME' #
index_arr = [0]#,1,3,4,5,7] #np.arange(8) #np.arange(len(fo)) #

for ii,k in enumerate(index_arr):
    #try:
    wc_file = path[k] + 'database/ec' + 'thar2.ec' #'Arcs_sum.tu.ec'# #thar2.ec is Ayan's new wavelength calibrations; the file itself is in air wavelength but it gets converted to vacuum in spec_utils.py
    filename = path[k] + fn[k]
    fout = outpath + fo[k]
    fout = fout[:-4]+'_esi.txt'
    z = z_arr[k]
    wrmin, wrmax = 3750./(1+z), 10200./(1+z) #default wavelength range 4700-10200 A in observed frame
    #print '\nStarting', filename
    #----------------------------------------------------------------
    s,h = u.get_fitsspec(filename)
    print fo[k]+'\t\t\t\t'+str(h['RA'])+'\t\t'+str(h['DEC'])+'\t\t'+str(h['DATE-OBS'])+'\t\t'+str(h['UT'])+'\t\t'+str(h['ELAPTIME']) #
    spec = pd.DataFrame(columns=['obswave','flam','flam_u']) #pandas dataframe to hold the spectrum
    wlist = u.get_dispsol(wc_file, ap_max, show=show_disp_sol, is_wave_air=is_wave_air) #wavelength solution being grabbed here
    wave = np.arange(wlist[0][1],wlist[-1][-2],0.2)
    #----------------------------------------------------------------
    fluxes, errors = [],[]
    for ap in range(ap_max):
        func1 = interp1d(wlist[ap,:],s[0,ap,:],fill_value=0.,bounds_error=False)
        fluxes.append(func1(wave))
        func2 = interp1d(wlist[ap,:],s[3,ap,:],fill_value=0.,bounds_error=False)
        errors.append(func2(wave))
    fluxes = np.array(fluxes)
    errors = np.array(errors)
    weights = 1./errors**2
    weights[weights == np.inf] = 0.001 #very small value sinstead of 0, so that weighting does not lead to division by 0, 
                                        #it will not affect results as, if all errors for certain wavelength value is zero, all fluxes will also be zero
    flux = np.average(fluxes,weights=weights,axis=0)
    error = np.average(errors,weights=weights,axis=0)
    #---------------------------------------------------------------- 
    spec['obswave'] = pd.Series(["{0:.2f}".format(w) for w in wave])
    spec['flam'] = pd.Series(flux)
    spec['flam_u'] = pd.Series(error)
    spec.obswave = spec.obswave.astype(np.float64)
    spec.flam = spec.flam.astype(np.float64)
    spec.flam_u = spec.flam_u.astype(np.float64)
    spec['restwave'] = spec['obswave']/(1.+z)
    if wrmin is None:
        wrmin = spec['restwave'].values[0]
    if wrmax is None:
        wrmax = spec['restwave'].values[-1]
    spec = spec[spec['restwave'].between(wrmin,wrmax)]
    #---------------------------------------------------------------- 
    if show_ord or show_ord_all:       
        if not show_ord_all or ii==0:
            fig = plt.figure(figsize=(18,4))
            fig.subplots_adjust(hspace=0.7, top=0.94, bottom=0.08, left=0.06, right=0.98)
        cmap = plt.get_cmap('rainbow')
        colors = [cmap(i) for i in np.linspace(0, 1, ap_max)]
        pix = np.arange(nypix)+1
        if show_before_dispsol:
            plt.plot(pix,smooth(s[0,order_to_show-1,:], boxcar=11),c=colors[order_to_show-1], linestyle = ['solid', 'dotted', 'dashed'][ii%3] if show_ord_all else 'solid', label='order '+str(order_to_show))
            plt.plot(pix,smooth(s[3,order_to_show-1,:], boxcar=11),c=colors[order_to_show-1], linewidth=0.3, linestyle = ['solid', 'dotted', 'dashed'][ii%3] if show_ord_all else 'solid')
            plt.xlabel('y-pixel')
            #plt.xlim(1950,2050) #
            plt.ylabel('Counts')
        else:
            for ap in range(ap_max):            
                plt.plot(wlist[ap],smooth(s[0,ap,:], boxcar=101),c=colors[ap], linestyle = ['solid', 'dotted', 'dashed'][ii%3] if show_ord_all else 'solid', label=fo[k] if ap==0 else None)
                plt.plot(wlist[ap],smooth(s[3,ap,:], boxcar=101),c=colors[ap], linewidth=0.3, linestyle = ['solid', 'dotted', 'dashed'][ii%3] if show_ord_all else 'solid')

            plt.xlabel('Observed wavelength (A)')
            plt.xlim(wrmin*(1+z), wrmax*(1+z))
            #plt.xlim(8000,9500) #
        
            if 'fluxcal' in fn[k]:
                plt.ylim(0,8e-17)
                plt.ylabel('flam(ergs/cm^2/s/A)')
            else: 
                if 'std' in filename: plt.ylim(0, 6e4) #for plotting non-flux-calibrated spectra
                else: plt.ylim(0, 500)
                plt.ylabel('Counts')
        plt.legend()
        plt.show(block=False)
    #----------------------------------------------------------------
    if mask_bad_col:
        for mask_lim in mask_lim_arr[k]:
            spec.ix[spec['obswave'].between(mask_lim[0], mask_lim[1]), 'flam'] = 0.0
            spec.ix[spec['obswave'].between(mask_lim[0], mask_lim[1]), 'flam_u'] = 99.
    #---------------------------------------------------------------- 
    if clean:
        if check: 
            u.plot_spec(spec, 'masking', check=True, sig_thresh=sig_thresh) #one plot to remain visible
            u.plot_spec(spec, 'masking', check=True, sig_thresh=sig_thresh) #aother to show the masks overlaid
            p= plt.axvline(spec['restwave'].values[0], c='r')
        spec['badmask']  = False
        cleanspec(spec)
    if write:
        u.writespec_txt(spec, fout, z=z, filename=filename)
    #---------------------------------------------------------------- 
    if show_full or show_full_all:
        #spec = u.get_txtspec(fout, wrmin=wrmin, wrmax=wrmax, sig_thresh=sig_thresh, check=check, clean=clean)
        if not show_full_all or ii==0:
            fig = plt.figure(figsize=(18,4))
            fig.subplots_adjust(hspace=0.7, top=0.94, bottom=0.08, left=0.06, right=0.98)
        p = plt.plot(spec['obswave'], smooth(spec['flam'], boxcar=51), label=fo[k])
        plt.plot(spec['obswave'], smooth(spec['flam_u'], boxcar=51), c=p[0].get_color(), linewidth=0.3)
        plt.xlabel('Observed wavelength (A)')
        plt.xlim(wrmin*(1+z), wrmax*(1+z))
        plt.xlim(3750,10e3) #

        if 'fluxcal' in fn[k]:
            plt.ylim(0,6e-17)
            plt.ylabel('flam(ergs/cm^2/s/A)')
        else: 
            if 'std' in filename: plt.ylim(0, 1e5) #for plotting non-flux-calibrated spectra
            else: plt.ylim(0, 500)
            plt.ylabel('Counts')
        plt.legend()
        plt.show(block=False)
    '''    
    except(e):
        print e
        print 'Failed at', fn[k]
        continue
    '''
#---------------------------------------------------------------- 
print('Done in %s minutes' % ((time.time() - start_time)/60))
