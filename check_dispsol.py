#utility routine for checking if dispersion solution obtained using IRAF for ESI is good enough and 
#comparing it with dispersion solution from ESI website

import numpy as np
import matplotlib.pyplot as plt
import spec_utils as u
#-------------------From Jane Rigby's esi_wavesol.py---------------------------------
#Dictionary of ESI wavelength coefficients, from http://www2.keck.hawaii.edu/inst/esi/QuickRef.html#Wavelengths
ed = {}
ed['15'] = np.array((4077.46, 0.154482, -1.140e-6, -3.106e-10))
ed['14'] = np.array((4366.24, 0.165128, -2.0205e-6, 5.71e-10))
ed['13'] = np.array((4699.50, 0.179043, -1.912e-6, -8.44e-11))
ed['12'] = np.array((5088.55, 0.194456, -2.140e-6, 4.00e-11))
ed['11'] = np.array((5549.09, 0.212052, -2.365e-6, -1.23e-10))
ed['10'] = np.array((6101.46, 0.233675, -2.593e-6, -1.105e-10))
ed['9']  = np.array((6776.99, 0.259847, -2.826e-6, -1.90e-10))
ed['8']  = np.array((7621.60, 0.29266, -3.203e-6, -2.77e-10))
ed['7']  = np.array((8707.59, 0.334496, -3.6815e-6, -2.58e-10))
ed['6']  = np.array((10156,   0.39, -4.25e-6, 0))

def approx_wave_esi(order, Yval) :
    ''' arguments are order (6 through 15) and Y pixel coordinate 
    returns approx wavelength in Angstroms'''
    coeff = ed[order]
    wave = coeff[0] + coeff[1]*(Yval-2048) + coeff[2]*(Yval-2048)**2 + coeff[3]*(Yval-2048)**3
    return(wave)
#------------------------------------------------------------------------------------
ap_max, z, restwave, ymark = 8, 1.3294, 3729.9, 2006
wc_file = '/Users/acharyya/Documents/esi_2016b/2016aug27_1x1/IRAF/database/ecarcs.tr.ec'
wlist = u.get_dispsol(wc_file, ap_max, show=True)

wc_file = '/Users/acharyya/Documents/esi_2016b/2016aug27_2x1/IRAF/database/ecArcs_sum.tr.ec'
wlist = u.get_dispsol(wc_file, ap_max, show=True, col='r')

ypix = np.arange(1,4096+1)
for o in range(6,13+1):
    w=[]
    for y in ypix:
        w.append(approx_wave_esi(str(o), y))
    plt.plot(ypix,w,c='g')
plt.axvline(ymark, linestyle='--', c='k')
plt.axhline(restwave*(1+z), linestyle='--', c='k')
plt.show(block=False)