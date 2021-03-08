#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:32:29 2021

@author: mareya
"""

import aperture_tools as ap
import A1_tools as a1
import numpy as np
import matplotlib.pyplot as plt
import time


def variable_aperture2(data,min_brightness,min_radius,a,b,cal,e_cal):
    """
    
    Parameters
    ----------
    data : Dataset as array
    min_brightness : The minimum pixel value for which will be counted as a source
    a : Difference between inner annulus radius and aperture radius
    b : Difference between outer annulus radius and aperture radius
    cal : calibration stored as MAGZPT
    e_cal : error on calibration stored as MAGZRR

    Returns
    -------
    data : Dataset after masks
    count : Galaxy count
    

    """
    
    count = 0
    
    storage = []
    
    while np.amax(data) > min_brightness:
        
        # Save current max value 
        max_value = np.amax(data)
        
        # Zip into list of array indices for these max value(s)
        l = list(zip(*np.where(data == max_value)))
        
        for i in range(len(l)):
            
            r = ap.find_radius(data,l[i][0],l[i][1],min_brightness)
            
            if r < min_radius:
                
                # Too small to be a galaxy
                # Cover up and continue
                data = ap.data_after_mask(data,l[i][0],l[i][1],r)
                
                #print('\nNot big enough to be a galaxy')
                
            else:
            
                # Local background
                bkg,e_b = ap.annulus_background(data,l[i][0],l[i][1],r+a,r+b)
            
                # Total flux in aperture, corrected with local background
                c_flux,e_c = ap.corrected_flux(data,l[i][0],l[i][1],r,r+a,r+b)
            
                data = ap.data_after_mask(data,l[i][0],l[i][1],r)
                
                if c_flux > 0:
                    
                    # Galaxy found, add to count
                    count += 1
                    
                    # Calculate magnitude
                    m = cal - 2.5*np.log10(c_flux)
                    
                    # Calculate the error on the magnitude
                    e_m = np.sqrt((e_cal)**2 + ((-(5)/(2*np.log(10)*c_flux))*(e_c))**2)
                    
                    
                    #print('\nGalaxy at (%s,%s)! Count = %s. Mag = %s. Local bkg = %s.' % (l[i][1],l[i][0],c_flux,m,bkg))
                    storage.append(np.array([count,l[i][1]+1,l[i][0]+1,c_flux,e_c,bkg,e_b,m,e_m]))
                    
                #else: 
                    #print('\nWeirdness at (%s,%s). Count = %s. Local bkg = %s.' % (l[i][1],l[i][0],c_flux,bkg))
    
    
    np.savetxt('catalogue_3465_3_9.txt', storage, header="ID\tx\ty\tFlux\tFlux error\tBkg\tBkg error\tm\tm error",delimiter='\t',fmt='%0.4e')

    return data,count

# =============================================================================
# TEST ON PATCH WITH ONLY 1 GALAXY
# =============================================================================

# Load clean data for entire image
data = ap.load_clean_data('mosaic.fits')

# Define a small patch with only 1 galaxy
#small_patch = ap.patch(data,390,490,830,930,0)

small_patch = ap.patch(data,410,460,860,910,0)

#%%

cal = ap.load_keyword_value('mosaic.fits', 'MAGZPT')
cal_error = ap.load_keyword_value('mosaic.fits', 'MAGZRR')

start_time = time.time()
data_small,count = variable_aperture2(small_patch,3465,3,3,9,cal,cal_error)
print('Number of galaxies fcalound = %s' %count)
fin = time.time() - start_time
print("Total time = %s seconds ---" % fin)

#%%

plt.imshow(small_patch,interpolation='nearest',origin='lower')
plt.xlabel('x (pixels)')
plt.ylabel('y (pixels)')
cbar = plt.colorbar()
cbar.set_label('Pixel value')
#plt.savefig('star-1-galaxy-before.png',bbox_inches='tight')
plt.show()

#%%

plt.imshow(data_small,interpolation='nearest',origin='lower')
plt.xlabel('x (pixels)')
plt.ylabel('y (pixels)')
cbar = plt.colorbar()
cbar.set_label('Pixel value')
#plt.savefig('star-1-galaxy-after.png',bbox_inches='tight')
plt.show()

#%%

# =============================================================================
# TEST ON PATCH WITH MANY GALAXIES
# =============================================================================

# Load clean data for entire image
#data = ap.load_clean_data('mosaic.fits')

# Define a bigger patch with a few galaxies
small_patch1 = ap.patch(data,650,970,1610,2070,20)

#%%

start_time = time.time()
data_small1,count1 = variable_aperture3(small_patch1,3465,3,3,9,cal,cal_error)
print('Number of galaxies found = %s' %count1)
fin = time.time() - start_time
print("Total time = %s seconds ---" % fin)

#%%

plt.imshow(np.log(small_patch1),interpolation='nearest',origin='lower')
plt.xlabel('x (pixels)')
plt.ylabel('y (pixels)')
cbar = plt.colorbar()
cbar.set_label('ln(pixel value)')
#plt.savefig('star-14-before.png',bbox_inches='tight')

#%%

plt.imshow(np.log(data_small1),interpolation='nearest',origin='lower')
plt.xlabel('x (pixels)')
plt.ylabel('y (pixels)')
cbar = plt.colorbar()
cbar.set_label('ln(pixel value)')
#plt.savefig('star-14-after.png',bbox_inches='tight')


#%%
# =============================================================================
# RUN ON ENTIRE IMAGE
# =============================================================================

# Load clean data for entire image
data = ap.load_clean_data('mosaic.fits')


plt.imshow(data,interpolation='nearest',origin='lower')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()

#%%

mu = 3418.564628830125
sigma = 11.700641756493768
min_value = mu + 4*sigma

#%%
final_data,total_count = variable_aperture2(data,min_value,3,3,9,cal,cal_error)
print('Number of galaxies found = %s' %total_count)

#%%

plt.imshow(np.log(final_data),interpolation='nearest',origin='lower')
plt.title('Before (log)')
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
#plt.savefig('after-small-patch-log.png',bbox_inches='tight')


