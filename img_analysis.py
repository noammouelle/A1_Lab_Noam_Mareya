# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 18:36:21 2021

@author: nomou
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import A1_tools as a1
from scipy.optimize import curve_fit

#%% Loading the modified image and printing basic info

img = np.loadtxt("modified_image.txt")
#%%
N = len(img[0])*len(img)
n = img.reshape(N).tolist().count(0)

print('Number of pixels in the image: '+str(N)+'\nNumber of masked pixels: '+str(n)
       ,'\nFraction of masked pixels: '+str(n/N))

#%% Plot the histogram to estimate average background and standard deviation

# Flattening the image onto a single axis
img_flattened = img.reshape(len(img)*len(img[0]))

# Creating a histogram
lo = 3350; hi = 3600
img_flattened = img_flattened[(img_flattened>=lo) & (img_flattened<=hi)]
occurences, vals0, others = plt.hist(img_flattened,
                                   bins=len(np.unique(img_flattened)),
                                   histtype='step',color='black',
                                   label='Data')

# Removing last bin edge for plotting purposes
vals = vals0[:len(vals0)-1]

# Fitting a Gaussian to the distribution
popt, pcov = curve_fit(a1.gauss,vals,occurences,p0=[3e5,3420,13])
err = np.sqrt(np.diag(pcov))

# Printing the results of the fit
avg_bg = popt[1]
std_bg = popt[2]

print('AVG +/- STD: '+'%.3f'%avg_bg+' +/- '+'%.3f'%std_bg)
print('Uncertainty on AVG: '+'%.3f'%err[1]+'\nUncertainty on STD: '+'%.3f'%err[2])

# Plotting the corresponding curve
x =  np.linspace(min(vals),max(vals),1000)
plt.plot(x,a1.gauss(x,*popt),'--',color='black',label='Gaussian fit',linewidth=1)
plt.legend()

plt.title('Background pixel value distribution')
plt.xlabel('Value')
plt.ylabel('Occurences')
plt.xlim(3380,3460)

#%% Defining the thresholds

t1 = avg_bg + 4*std_bg # Detection threshold
t2 = avg_bg + 3.5  *std_bg # Expansion threshold

#%% Detecting objects

img2 = np.copy(img)

positions    = []
num_pixels   = []
photon_count = []
flattening   = []

r0 = a1.find_brightest(img2) # Current brightest pix
n = 0
while img2[r0] > t1:
    # Keep track of the loop
    print(n)
    n += 1
    # Find the indices of the pixels, radius of the source and peak count
    indices, is_at_boundary = a1.detect_object(img2,t2,r0)
    r = a1.find_radius(img2,t2,r0)
    
    if (r >= 3) and (img2[r0]<8000) and (is_at_boundary == False):
        # Store the number of pixels in the source
        num_pixels.append(len(indices))
        # Store the position of the source
        positions.append(r0)
        # Store the total photon counts from the source
        photon_count.append(a1.photon_count(img2,indices))
        # Store the flattening
        flattening.append(a1.find_flattening(img2,t2,r0))
        
    # Delete from the image
    a1.delete(img2,indices)
    # Update r0
    r0 = a1.find_brightest(img2)

print('Number of objects detected: '+str(len(positions)))

#%% Save the data
# Columns : i   j   Number of pixels    Photon Count    Flattening
positions_arr = np.transpose(np.array(positions))
i = positions_arr[0]
j = positions_arr[1]

num_pixels_arr     = np.array(num_pixels)
photon_count_arr   = np.array(photon_count)
flattening_arr     = np.array(flattening)

raw_stats     = np.transpose(np.array([i,j,num_pixels_arr,photon_count_arr,flattening_arr]))

np.savetxt("run2_raw_stats.txt",raw_stats)

#%% Calibration of all the sources
img_og = np.copy(img)
median_bg = []
for k in range(len(raw_stats)):
    # Keep track of the loop
    print(k)
    # Load the kth line of the data
    object_k_stats = raw_stats[k]
    
    img_temp = np.copy(img_og)
    r = (int(object_k_stats[0]),int(object_k_stats[1]))
    
    R_obj = a1.find_radius(img_temp,t2,r)
    R_disk = 3*R_obj
    
    pixels_obj, other = a1.detect_object(img_temp,t2,r)
    a1.delete(img_temp,pixels_obj)
    
    pixel_vals = a1.vals_disk(img_temp,R_disk,r)
    
    pixel_vals = pixel_vals[pixel_vals != 0]
    
    median_bg.append(np.median(pixel_vals))

#%% Modify the original data
calibrated_stats = np.c_[np.copy(raw_stats),np.zeros(len(raw_stats))]

for k in range(len(raw_stats)):
    # Obtain the photometry of each source
    calibrated_stats[k,3] = raw_stats[k,3] - raw_stats[k,2]*median_bg[k]
    
    # Add magnitude to the table
    mag = 2.530e1 - 2.5*np.log10(calibrated_stats[k,3])
    if (mag != np.nan) and (mag != np.inf):
        calibrated_stats[k,5] = mag
    else:
        calibrated_stats = np.delete(calibrated_stats,(k),axis=0)

# Print the number of sources which were discarded
print('Number of discarded sources: '+str(len(raw_stats)-len(calibrated_stats)))
# Save the calibrated data
np.savetxt('run2_calibrated_stats.txt',calibrated_stats)
