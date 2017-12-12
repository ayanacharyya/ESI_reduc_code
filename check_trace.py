#utility routine for checking if trace obtained using IRAF is good enough

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from scipy.interpolate import interp1d
import spec_utils as u
import warnings
warnings.filterwarnings("ignore")
#----------------------------------------------------------------
fig = plt.figure(figsize=(8,8))
fig.subplots_adjust(hspace=0.7, top=0.94, bottom=0.05, left=0.1, right=0.95)
input = sys.argv[1] if len(sys.argv)>1 else 'std1_cos_clean'
filename = '/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_1x1/IRAF/'+input+'.fits'#'/Users/acharyya/Work/astro/obs_data/esi_2017a/2017feb24/IRAF/P0440_23_cos.fits'#
s,h = u.get_fitsspec(filename)
binx, biny = 1,1
#----------------------------------------------------------------        
count = 0
ap_file = '/Users/acharyya/Work/astro/obs_data/esi_2016b/2016aug27_1x1/IRAF/database/ap'+input#'/Users/acharyya/Work/astro/obs_data/esi_2017a/2017feb24/IRAF/database/aplast'#
lines = open(ap_file,'r').readlines()
o,left,right,ymin,ymax,coeff,off=[],[],[],[],[],[],[]
for i,l in enumerate(lines):
    if '\taperture\t' in l:
        off.append(float(lines[i+2].split()[1]))
        left.append(float(lines[i+3].split()[1]))
        right.append(float(lines[i+4].split()[1]))
        o.append(float(str(lines[i+19]).strip()))
        ymin.append(float(str(lines[i+20]).strip()))
        ymax.append(float(str(lines[i+21]).strip()))
        c=[]
        for j in range(int(o[-1])): c.append(float(str(lines[i+22+j]).strip()))
        coeff.append(c)
        count += 1
#------------------------------------------------------------- 
def ap(x,o,xmin,xmax,c):
    o = int(o)
    n=(2*x-(xmax+xmin))/(xmax-xmin)
    z=[1,n]
    for i in range(1,o-1):
        z.append(((2*i+1) * n * z[i] - i * z[i-1]) / (i+1)) 
    s = np.sum(np.multiply(c,z))
    return s
#---------------------------------------------------------------- 
#y = float(sys.argv[1]) if len(sys.argv)>1 else 1484
#print 'For y=',y,', x=',ap(y,o,xmin,xmax,c)
if 'std' in filename: vmin,vmax = 1.8, 2.2
elif 'Flat' in filename: vmin, vmax = 2.3, 3.
else: vmin, vmax = 1.5, 1.7

plt.imshow(np.log10(s),cmap='gray',aspect='auto',vmin=vmin, vmax=vmax)
for i in range(count):
    y = np.linspace(1,4096/biny,1000) #np.linspace(ymin[i],ymax[i],1000)
    aleft,aright=[],[]
    for j in range(len(y)):
        aleft.append(ap(y[j],o[i],ymin[i],ymax[i],coeff[i])+off[i]+left[i])
        aright.append(ap(y[j],o[i],ymin[i],ymax[i],coeff[i])+off[i]+right[i])
    p=plt.plot(aleft,y,lw=0.3)
    plt.plot(aright,y,c=p[0].get_color(),lw=0.3)
plt.ylim(0,4096/biny)
plt.xlim(0,2232/binx)
print 'Trace for', filename
plt.show(block=False)