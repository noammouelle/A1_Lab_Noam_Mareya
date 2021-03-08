# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:08:30 2021

@author: nomou

Contains useful functions
"""
import numpy as np
from astropy.io import fits

def gauss(x,A,mu,std):
    return A*np.exp(-(x-mu)**2/(2*std**2))

def logn(n,a):
    return np.log(a)/np.log(n)

def load_img(filename):
    return fits.open(filename)[0].data

def alive(img,threshold,r):
    '''
    In: i,j
    Out: lives = True if above threshold, zero = True is zero
    '''
    lives = False
    zero = False
    if img[r] >= threshold:
        lives = True
    elif img[r] == 0:
        zero = True
    
    return lives, zero

def find_brightest(img):
    '''
    Finds the indices of the brightest object in the image.
    '''
    result = np.where(img==np.amax(img))
    imax = result[0][0]
    jmax = result[1][0]

    return (imax, jmax)

def check_neighbours(img,threshold,r):
    '''
    Returns an array containing tuple indices of all the live neighbours of a
    cell
    '''
    # Returns true if a zero was found among the neighbours
    found_zero = False
    
    imax = len(img)-1
    jmax = len(img[0])-1
    neighbours = []
    i0 = r[0] ; j0 = r[1]
    # Check each neighbour in turn
    i = i0 + 1 ; j = j0
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    i = i0 - 1 ; j = j0
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    i = i0     ; j = j0 + 1
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    i = i0     ; j = j0 - 1
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    i = i0 + 1 ; j = j0 + 1
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    i = i0 - 1 ; j = j0 + 1
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    i = i0 + 1 ; j = j0 - 1
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    i = i0 - 1 ; j = j0 - 1
    if i<=imax and j<=jmax:
        lives, zero = alive(img,threshold,(i,j))
        if lives == True:
            neighbours.append((i,j))
        elif zero == True:
            found_zero = True
        
    return neighbours, found_zero
        

def detect_object(img,threshold,r0):
    '''
    Takes the brightest pixel and exmpand a region around it until the boundary
    of the object is found.
    '''
    # Avoids boundary effects
    is_at_boundary = False
    # Will contain the tuple coordinates of all "live" cells
    live_cells = [r0]    
    # Will contain the coordinates of all the newly born cells so that
    # neighbours can be checked for
    new_cells = [r0]
    
    n = len(live_cells)
    n_new = 2
    while n_new > n:
        n = n_new
        new_cells_temp = []
        for k in range(len(new_cells)):
            r = new_cells[k]
            # Check neighbours for kth cell
            neighbours_cell_k, found_zero = check_neighbours(img,threshold,r)
            if found_zero == True: is_at_boundary = True
            # Update the list of newly born cells, and avoid duplicates
            new_cells_temp += neighbours_cell_k
            new_cells_temp = list(dict.fromkeys(new_cells_temp))
            
        new_cells = new_cells_temp
        live_cells += new_cells
        live_cells = list(dict.fromkeys(live_cells))
        n_new = len(live_cells)
    
    return live_cells, is_at_boundary

def delete(img,indices):
    '''
    Sets the pixels at site indicated by "indices" to zero.
    '''
    for k in range(len(indices)):
        img[indices[k]] = 0
    
def photon_count(img, indices):
    '''
    Returns the average value of the pixels indicated by indices
    '''
    N = len(indices)
    tot = 0
    for i in range(N):
        tot += img[indices[i]]
    return tot

def circle(R,crit):
    '''
    Outputs the indices of pixels forming a circle of radius R
    '''
    i_arr, j_arr = [], []
    for i in range(R+1):
        for j in range(R+1):
            if np.abs(i**2+j**2 - R**2) <= crit:
                i_arr.append(i)
                j_arr.append(j)
                
    i_arr, j_arr = np.array(i_arr),np.array(j_arr)
    
    out_i = np.append(i_arr,[-i_arr,i_arr,-i_arr])
    out_j = np.append(j_arr,[j_arr,-j_arr,-j_arr])
    
    return out_i, out_j
    
def local_median(img,R,r):
    '''
    Defines an annulus of radius R around position r and returns the median
    value of those pixels
    '''
    i0, j0 = r[0], r[1]
    indices_i, indices_j = circle(R,R)
    
    imax = len(img)-1
    jmax = len(img[0])-1
    
    i_temp, j_temp = i0+indices_i, j0+indices_j
    i,j = [],[]
    for ii in range(len(i_temp)):
        if (0<i_temp[ii]<=imax) and (0<j_temp[ii]<=jmax):
            i.append(i_temp[ii])
            j.append(j_temp[ii])
            
    i,j = np.array(i),np.array([j])
            
    vals = img[i, j]
    return np.median(vals)

def detect_object_2(img,R,r):
    '''
    Takes the brightest pixel and exmpand a region around it until the boundary
    of the object is found.
    Threshold is now defined using the annulus method.
    '''
    threshold = local_median(img,R,r) + 5*13
    live_cells = detect_object(img,threshold,r)
    
    return live_cells

def find_orientation(img,threshold,r):
    '''
    Returns the orientation of an object centered around (i0,j0)
    0 :  S - N
    1 :  W - E
    2 : SW - NE
    3 : SE - NW
    '''
    ceil = threshold
    i0 = r[0]
    j0 = r[1]
    
    imax, jmax = len(img)-1,len(img[0])-1
    
    orientation = 0
    
    pixel_string = []
    
    val0 = img[i0,j0]
    
    # Check in vertical direction
    temp_pixel_string = []
    
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        i += 1
        val = img[i,j]
    
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        i -= 1
        val = img[i,j]
    
    pixel_string = temp_pixel_string
    
    # Check in horizontal direction
    temp_pixel_string = []
    
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        j += 1
        val = img[i,j]
    
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        j -= 1
        val = img[i,j]
        
    if len(temp_pixel_string) > len(pixel_string):
        pixel_string = temp_pixel_string
        orientation = 1
    
    # Check in SW - NE
    temp_pixel_string = []
    
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        i += 1
        j += 1
        val = img[i,j]
        
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        i -= 1
        j -= 1
        val = img[i,j]
        
    if len(temp_pixel_string) > len(pixel_string):
        pixel_string = temp_pixel_string
        orientation = 2
        
    # Check in SE - NW
    temp_pixel_string = []
    
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        i += 1
        j -= 1
        val = img[i,j]
        
    i, j = i0, j0
    val = val0
    while (val > ceil) and (i<imax) and (j<jmax):
        temp_pixel_string.append((i,j))
        i -= 1
        j += 1
        val = img[i,j]
        
    if len(temp_pixel_string) > len(pixel_string):
        pixel_string = temp_pixel_string
        orientation = 3
    
    return orientation, pixel_string

def find_radius(img,threshold,r):
    '''
    Returns the min radius containing the source fully
    '''
    orientation, pixel_string = find_orientation(img,threshold,r)
    l = (len(pixel_string)-1)/2
    if (orientation == 2) or (orientation == 3):
        R = np.round(l*np.sqrt(2))
    else:
        R = l
    '''
    if R<6:
        R = 6
    '''
    return int(R)

def vals_disk(img,R0,r0):
    '''
    Returns an array with the values of the pixels within the disk of radius R
    '''
    R = int(np.round(R0))
    i0 = r0[0] ; j0 = r0[1]
    vals = []
    for i in range(-R,R+1):
        for j in range(-R,R+1):
            if i**2 + j**2 <= R**2:
                vals.append(img[i0+i,j0+j])
    return np.array(vals)

def find_flattening(img,threshold,r):
    '''
    Finds the flattening of the object: (a-b)/a
    '''
    orientation, pix = find_orientation(img,threshold,r)
    a = (len(pix)-1)/2
    
    i0, j0 = r[0],r[1]
    val0 = img[i0,j0]
    
    pixb = []
    i, j = i0, j0
    val = val0
    
    if orientation == 0:
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            j += 1
            val = img[i,j]
    
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            j -= 1
            val = img[i,j]
            
    if orientation == 1:
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            i += 1
            val = img[i,j]
    
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            i -= 1
            val = img[i,j]
            
    if orientation == 2:
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            i += 1
            j -= 1
            val = img[i,j]
        
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            i -= 1
            j += 1
            val = img[i,j]
    
    if orientation == 3:
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            i += 1
            j += 1
            val = img[i,j]
        
        i, j = i0, j0
        val = val0
        while val > threshold:
            pixb.append((i,j))
            i -= 1
            j -= 1
            val = img[i,j]
        
    l = (len(pixb)-1)/2
    if orientation == 3 or orientation == 4:
        b = np.round(np.sqrt(2)*l)
    else:
        b = l
    
    return (a-b)/a