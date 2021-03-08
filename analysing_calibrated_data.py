# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 08:06:17 2021

@author: nomou
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import A1_tools as a1
from scipy.optimize import curve_fit
from scipy.special import erf
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs import WCS
#%% Loading the image
img = np.loadtxt('modified_image.txt')
mu = 3.41838272e+03
noise = 1.09965805e+01
t2 = mu+3.5*noise
gain = 3.1
#%% Loading the calibrated data
fname = 'run1_calibrated_stats.txt'
new_sources = np.transpose(np.loadtxt(fname,unpack=True))

#%% Classify between stars and galaxies
flattening = np.transpose(new_sources)[4]
galaxies = []
stars    = []
for k in range(len(new_sources)):
    flattening = new_sources[k,4]
    i0,j0 = int(new_sources[k,0]), int(new_sources[k,1])
    I0 = img[i0,j0]
    if flattening >= 0.3:
        galaxies.append(new_sources[k])
    else:
        stars.append(new_sources[k])

galaxies = np.array(galaxies)
stars    = np.array(stars)
# Print number in each class
print('Number of objects classified as stars: '+str(len(stars))+'\nNumber of objects classified as galaxies: '+str(len(galaxies)))

#%% Galaxy umber count plot
m = np.transpose(galaxies)[5]
npix = np.transpose(galaxies)[2]
flux = np.transpose(galaxies)[3]

# Create the histogram
values, base = np.histogram(m,range=(11,19),bins=32)
values = values/0.04111 # area in square degree
# Create the cumulative histogram
cumul = np.cumsum(values)
err = np.sqrt(cumul)

# Calculating uncertainties
err_m = np.sqrt(0.02**2+2.36*npix*noise**2/flux**2+0.38/flux)

var_bins = []
for k in range(len(base)):
    var_bink = 0
    for i in range(len(m)):
        qi = 0.5*err_m[i] * (erf((base[k]-m[i])/(np.sqrt(2)*err_m[i]))+erf((m[i])/(np.sqrt(2)*err_m[i])))
        var_bink += qi*(1-qi)
    var_bins.append(var_bink)

err_bins = np.sqrt(np.array(var_bins))/0.04111

# Remove bins with less than 3 objects
wanted_indx = np.where(cumul>1)
cumul = cumul[wanted_indx]
base  = base[wanted_indx]
err_bins = err_bins[wanted_indx]


# Take the log10 of N(m)
y = np.log10(cumul)

# Attempt to fit the data
def f(m,a,c):
    return a*m + c

# Plot the number count    

fig, ax = plt.subplots(constrained_layout=True)

ax.errorbar(base, cumul,yerr=err_bins,capsize=2,ecolor='black',fmt='.',
            mec='black',mfc='none',ms=3,elinewidth=1,label='Data')
ax.set_yscale('log')

yticks = np.logspace(1,4,5)
ytick_labels = []
for k in range(len(yticks)):
    ytick_labels.append('%.1f'%np.log10(yticks[k]))
    
ax.set_yticks(yticks)
ax.set_yticklabels(ytick_labels)

ax.minorticks_off()
ax.set_title('R-band galaxy number count')

def f1(x):
    return x
def f2(x):
    return x

second_yaxis = ax.secondary_yaxis('right', functions=(f1, f2))
second_yaxis.set_ylabel('N [0.25 mag$^{-1}$ deg$^{-2}$]')

second_yticks = np.linspace(0,400,20)
second_ytick_labels = []
for k in range(len(second_yticks)):
    second_ytick_labels.append('%.1f'%second_yticks[k])
    
second_yaxis.set_yticks(second_yticks)
second_yaxis.set_yticklabels(second_ytick_labels)

base_1 = base[:4]
y_1 = y[:4]

popt,pcov = curve_fit(f,base_1,y_1)

x = np.linspace(base_1[0],base_1[-1]+0.2,1000)
ax.plot(x,10**f(x,*popt),'--',color='black',linewidth=0.7,label='a=1.20')

base_2 = base[5:-5]
y_2 = y[5:-5]

popt,pcov = curve_fit(f,base_2,y_2)

x = np.linspace(base_2[0]-1,base_2[-1]+1,1000)
ax.plot(x,10**f(x,*popt),'-.',color='black',linewidth=0.7,label='a=0.33')

print(popt)
ax.set_xlabel('R mag (AB)')
ax.set_ylabel('log$_{10}$(N)')
ax.legend()

#%% Creating catalog
radii = []
for k in range(len(galaxies)):
    i0, j0 = int(galaxies[k,0]),int(galaxies[k,1])
    radii.append(a1.find_radius(img,t2,(i0,j0))*0.258)
radii = np.array(radii)

# Convert to sky coords
hdulist = fits.open('mosaic.fits')
w = WCS(hdulist[0].header)
i_indx, j_indx = np.transpose(galaxies)[0],np.transpose(galaxies)[1]
sky_coords = pixel_to_skycoord(i_indx,j_indx,w)
RA, Dec = [],[]
for i in range(len(sky_coords)):
    RA.append(sky_coords[i].ra.deg)
    Dec.append(sky_coords[i].dec.deg)

m = np.transpose(galaxies)[5]
flattening = np.transpose(galaxies)[4]

np.savetxt('catalog.txt', np.transpose(np.array([RA,Dec,m,err_m,radii,flattening])),
           header='RA (deg),Dec (Deg),R-Mag (AB),Mag uncertainty (AB),Radius (as), Flattening',
           delimiter=',',fmt='%1.2f')

# Shorter version for report
np.savetxt('catalog_30.txt', np.transpose(np.array([RA,Dec,m,err_m,radii,flattening]))[:30],
           header='RA (deg),Dec (Deg),R-Mag (AB),Mag uncertainty (AB),Radius (as), Flattening',
           delimiter=',',fmt='%1.2f')

#%% Plotting the first 30 images
fig,axs =  plt.subplots(6,5)
for i in range(6):
    for j in range(5):
        R = radii[int(i*5+j)]
        i0,j0 = i_indx[int(i*5+j)],j_indx[int(i*5+j)]
        axs[i,j].pcolormesh(img[int(i0-3*R):int(i0+3*R),int(j0-3*R):int(j0+3*R)],
           cmap='Greys_r')
        axs[i,j].axis('off')
