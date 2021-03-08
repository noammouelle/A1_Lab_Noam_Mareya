#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 21:56:57 2021

@author: mareya
"""
#%%

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt

#%%

def load_data(filename):
    """
    Parameters
    ----------
    filename : FITS filename in string format

    Returns
    -------
    Data from FITS file

    """
    return fits.open(filename)[0].data


def load_clean_data(filename):
    """
    Parameters
    ----------
    filename : FITS filename in string format (i.e. "mosaic.fits")

    Returns
    -------
    data : Returns same size array with borders, stars and other artefacts 
    equal to 0
    
    """
    hdulist = fits.open(filename)
    data = hdulist[0].data
    hdulist.close()
    
    # Central bleed
    data[0:4611,1427:1447] = 0
    data[3136:3296,1384:1500] = 0
    data[2950:3475,1420:1454] = 0
    # Borders
    data[0:119,0:2570] = 0
    data[0:4611,0:117] = 0
    data[4514:4611,0:2570] = 0
    data[0:4611,2471:2570] = 0
    # Horizontal bleed, bottom
    data[217:228,1298:1538] = 0
    data[119:141,1289:1538] = 0
    data[231:249,1389:1477] = 0
    data[313:324,1026:1704] = 0
    data[313:358,1639:1650] = 0
    data[422:439,1101:1652] = 0
    data[422:455,1025:1046] = 0
    data[139:172,1406:1428] = 0
    data[139:151,1446:1456] = 0
    data[227:269,1416:1457] = 0
    data[323:334,1359:1503] = 0
    data[323:360,1445:1454] = 0
    data[438:452,1370:1428] = 0
    data[438:453,1446:1486] = 0
    data[333:361,1396:1428] = 0
    data[333:361,1442:1477] = 0
    data[359:390,1446:1456] = 0
    data[449:476,1409:1460] = 0
    # Stars
    data[4079:4119,547:573] = 0
    data[3292:3346,756:794] = 0
    data[3202:3419,771:780] = 0
    data[2748:2791,955:991] = 0
    data[2704:2835,967:979] = 0
    data[2221:2359,887:922] = 0
    data[4016:4050,1446:1467]=0
    data[3707:3803,2120:2149]=0
    data[3383:3445,2454:2573]=0
    data[2281:2339,2116:2148]=0
    data[1399:1455,2085:2093]=0
    data[1413:1437,2076:2102]=0
    # Middle star
    data[2862:3519,1113:1745] = 0
    # Star new
    data[4065:4142,518:601] = 0
    data[3206:3433,680:889] = 0
    data[2699:2848,898:1055] = 0
    data[2216:2370,831:981] = 0
    data[3704:3821,2084:2192] = 0
    data[2272:2355,2092:2178] = 0
    
    return data

# Return a value from a keyword in a FITS file
def load_keyword_value(filename,keyword):
    """
    Parameters
    ----------
    filename : FITS filename
    keyword : Keyword name whose value is desired, in string format

    Returns
    -------
    Value at keyword in FITS header
    
    """
    hdulist = fits.open(filename)
    header = hdulist[0].header
    hdulist.close()
    return header[keyword]

# Take a clean dataset and make a small patch with borders removed

def patch(data,ylow,yhigh,xlow,xhigh,size):
    """
    Parameters
    ----------
    data : Data as ndarray
    ylow : Lower boundary for patch on y-axis
    yhigh : Upper boundary for patch on y-axis
    xlow : Lower boundary for patch on x-axis
    xhigh : Upper boundary for patch on x-axis
    size : Width of border to be masked to 0

    Returns
    -------
    patch : Smaller rectangular patch of data to given dimensions

    """
    
    y_length = yhigh - ylow
    x_length = xhigh - xlow
    
    patch = data.copy()
    
    patch = patch[ylow:yhigh,xlow:xhigh]
    patch[0:size,0:x_length] = 0
    patch[0:y_length,0:size] = 0
    patch[y_length-size:y_length,0:x_length] = 0
    patch[0:y_length,x_length-size:x_length] = 0
    
    return patch

          
#%%

# FUNCTION THAT WILL RETURN AN ARRAY WITH THE SAME SIZE AS INPUT DATASET WITH 
# 1s IN A FILLED CIRCLE OF SPECIFIC LOCATION AND RADIUS AND 0s EVERYWHERE ELSE

def circle(data,y,x,r):
    """
    Parameters
    ----------
    data : Dataset with pixel values as an ndarray
    y : y-coordinate (i) of centre of circle
    x : x-coordinate (j) of centre of circle
    r : Desired radius of circle

    Returns
    -------
    circle_data : Logical 1/0 array of the same size as data with 1s filling
    in the circle of radius r and 0s everywhere else

    """
    # t is the tuple of max value (y,x)
    # y = l[0][0] = t[0]
    # x = l[0][1] = t[1]
    # where l is the list of max value array points
    # m = l[0] is a tuple of one array value
    # data[t] returns that max value

    square = []
    
    for i in range(int(-r),int(r+1)):
        for j in range(int(-r),int(r+1)):
            square.append((y+i,x+j)) # makes a list of tuples
    
    circle_data = data.copy()
    
    for i in range(len(square)):
        R = (y-square[i][0])**2 + (x-square[i][1])**2
        if R <= r**2:
            circle_data[square[i]] = 1
        else:
            circle_data[square[i]] = 0
    
    circle_data[circle_data != 1] = 0
    
    return circle_data

#%%

# FUNCTION THAT WILL RETURN TOTAL COUNT INSIDE A CIRCULAR APERTURE

def total_flux(data,y,x,r):
    """
    Parameters
    ----------
    data : Dataset with pixel values as an ndarray
    y : y-coordinate (i) of centre of circle
    x : x-coordinate (j) of centre of circle
    r : Desired radius of circle

    Returns
    -------
    total_flux : Sum of pixel values in circle of radius r centred at (x,y)
    aperture_area : Area of circular aperture

    """
    """
    square = []
    
    for i in range(int(-r),int(r+1)):
        for j in range(int(-r),int(r+1)):
            square.append((y+i,x+j)) # makes a list of tuples
    
    circle_data = data.copy()
    
    for i in range(len(square)):
        R = (y-square[i][0])**2 + (x-square[i][1])**2
        if R <= r**2:
            circle_data[square[i]] = 1
        else:
            circle_data[square[i]] = 0
    
    circle_data[circle_data != 1] = 0
    
    """
    
    circle_data = circle(data,y,x,r)
    
    data = np.where(circle_data == 0, 0, data)
    
    total_flux = np.sum(data)
    
    aperture_area = np.sum(circle_data)
    
    return total_flux, aperture_area


# FUNCTION THAT WILL TAKE A DATASET AND A PIXEL POINT, AND RETURN THE DATASET
# WITH A MASK OVER THE CIRCULAR APERTURE OF A GIVEN RADIUS (I.E. REPLACE PIXEL
# VALUES IN CIRCLE WITH 0)

def data_after_mask(data,y,x,r):
    """
    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    """
    # t is the tuple of max value (y,x)
    # y = l[0][0] = t[0]
    # x = l[0][1] = t[1]
    # where l is the list of max value array points
    # m = l[0] is a tuple of one array value
    # data[t] returns that max value
    
    square = []
    
    for i in range(int(-r),int(r+1)):
        for j in range(int(-r),int(r+1)):
            square.append((y+i,x+j)) # makes a list of tuples
    
    circle_data = data.copy()
    
    for i in range(len(square)):
        R = (y-square[i][0])**2 + (x-square[i][1])**2
        if R <= r**2:
            circle_data[square[i]] = 1
        else:
            circle_data[square[i]] = 0
    
    # Return array with only 1s where the circle is and 0s everywhere else
    circle_data[circle_data != 1] = 0
    
    data = np.where(circle_data==1,0,data)
        
    return data


def annulus_background(data,y,x,r1,r2):
    
    # Inner circle with radius r1
    circle_data_small = circle(data,y,x,r1)
    
    # Outer circle with radius r2
    circle_data_big = circle(data,y,x,r2)
    
    # Create a copy of big circle array
    annulus_data = circle_data_big.copy()
    
    # Replace 0s where small circle would be
    annulus_data = np.where(circle_data_small == 1, 0, circle_data_big)
    
    # Number of pixels in annulus
    annulus_area = np.sum(annulus_data)
    
    # Put zeros everywhere apart from annulus in data
    data = np.where(annulus_data == 0, 0, data)
    
    # Flatten into 1d array
    data1d = data.flatten()
    
    # Remove all zero values for just the annulus pixel values
    data1d = data[data != 0]
    
    # Filter data by clipping any values outside of 2 sigma
    data1d_filtered = sigma_clip(data1d, sigma=2)
    
    # Compress the sigma-clipped data set to remove masked elements
    data1d_filtered = data1d_filtered.compressed()
    
    # Find median of this
    med = np.median(data1d_filtered)
    
    # Calculate error on median by working out standard deviation of original 
    # unfiltered data and then dividing by annulus area
    error_b = np.std(data1d) / np.sqrt(annulus_area)
    
    return med, error_b


# Not using a sigma-clip

def annulus_background2(data,y,x,r1,r2):
    
    circle_data_small = circle(data,y,x,r1)
    
    circle_data_big = circle(data,y,x,r2)
    
    annulus_data = circle_data_big.copy()
    
    annulus_data = np.where(circle_data_small == 1, 0, circle_data_big)
    
    #annulus_area = np.sum(annulus_data)
    
    data = np.where(annulus_data == 0, 0, data)
    
    data1d = data.flatten()
    
    data1d = data[data != 0]
    
    med = np.median(data1d)
    
    return med


# Calculate corrected flux

def corrected_flux(data,y,x,r,r1,r2):
    t_flux, aperture_area = total_flux(data,y,x,r)
    
    local_background, error_b = annulus_background(data,y,x,r1,r2)
    
    c_flux = t_flux - (aperture_area * local_background)
    
    error_c = np.sqrt(aperture_area + aperture_area**2) * error_b
    
    return c_flux, error_c
    
    
    
#%%

# =============================================================================
# VARIABLE APERTURE
# =============================================================================

# Function to return data points in an annulus with all other values in 
# dataset = 0

def annulus(data,y,x,r1,r2):
    
    circle_data_small = circle(data,y,x,r1)
    
    circle_data_big = circle(data,y,x,r2)
    
    annulus_data = circle_data_big.copy()
    
    annulus_data = np.where(circle_data_small == 1, 0, circle_data_big)
    
    #annulus_area = np.sum(annulus_data)
    
    data = np.where(annulus_data == 0, 0, data)
    
    return data

#%%

def find_radius(data,y,x,min_brightness):
    
    # Positive x direction (right)
    r_pixel = (y,x)
    
    while data[r_pixel] > min_brightness:
        r_pixel = r_pixel[:1] + (r_pixel[1]+1,)
    
    """
    # Negative x direction (left)
    l_pixel = (y,x)
    
    while data[l_pixel] > min_brightness:
        l_pixel = l_pixel[:1] + (l_pixel[1]-1,)
        
    # Positive y direction (up)
    u_pixel = (y,x)
    
    while data[u_pixel] > min_brightness:
        u_pixel = (u_pixel[0]+1,) + u_pixel[1:]
        
    # Negative y direction (down)
    d_pixel = (y,x)   
        
    while data[d_pixel] > min_brightness:
        d_pixel = (d_pixel[0]-1,) + d_pixel[1:]
    
    
    furthest = max(r_pixel,l_pixel,u_pixel,d_pixel)
    """
    
    radius = np.sqrt((y-r_pixel[0])**2 + (x-r_pixel[1])**2)
    
    # Check that by increasing radius by 1, there will not be any values above
    # the min_brightness
    
    for i in range(100):
        annulus_values = np.where(annulus(data,y,x,radius,radius+1)==0,0,data)
        if np.amax(annulus_values) > min_brightness:
            radius +=1
        else:
            return radius
            
    return radius
    


