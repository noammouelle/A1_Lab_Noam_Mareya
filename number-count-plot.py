#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 22:50:24 2021

@author: mareya


Use this to load the catalogue.txt data and plot N(m) against m

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#%%

# =============================================================================
# READ IN CATALOGUE DATA AND CALCULATE N(m)
# =============================================================================

# Read in data from catalogue text file
id_no,x,y,flux,flux_error,bkg,bkg_error,m,err_m = np.loadtxt('catalogue_full_error_3465_3_9.txt', delimiter="\t", unpack=True)

# Calculate N(m) for each value of m
N_m = []

for i in range(len(m)):
    mag = m[i]
    
    count = 0
    for j in range(len(m)):
        
        # Add 1 to the count if a magnitude is less (brighter) than mag
        if m[j] < mag:
            count += 1
            
    N_m.append(count)
    
#%%

# =============================================================================
# PLOT SCATTER TO VIEW
# =============================================================================

y = np.log10(N_m)

plt.plot(m,y,'x',mew=2,ms=10,color='red')
plt.xlabel('m')
plt.ylabel('log(N(m))')
#plt.title('Small patch number count plot',y=1.02)
plt.show()

#%%

# =============================================================================
# BEST FIT LINE
# =============================================================================

def line(x, a, b):
    return a * x + b

y = np.log10(N_m)

mask = np.isfinite(y)
y = y[mask]

m = m[mask]

popt, pcov = curve_fit(line, m, y)

print("a =", popt[0], "+/-", pcov[0,0]**0.5)
print("b =", popt[1], "+/-", pcov[1,1]**0.5)

x = np.linspace(m[1], max(m), 1000)


plt.plot(m,y,'o',ms=8,mew=2, mfc='none',color='red',label='Galaxies')
plt.plot(x, line(x, popt[0], popt[1]), lw=2, color='black',label='Best fit line')
plt.xlabel(r'$\mathrm{m}$')
plt.ylabel(r'$\mathrm{log(N(m))}$')
#plt.xlim(10,18)
plt.title(r'$\mathrm{a = %.3f \pm %.3f, b = %.3f \pm %.3f}$' % (popt[0],pcov[0,0]**0.5,popt[1],pcov[1,1]**0.5),y=1.02)
plt.legend()
#plt.savefig('small-patch-number-count-plot.png',bbox_inches='tight')
plt.show()

#%%

N_m = np.asarray(N_m)

m_bright = np.where(m < 11.5)

bright_values = list(zip(*np.where(m <= 11.5)))

m_bright = (m[bright_values]).flatten()

y_bright = (y[bright_values]).flatten()

popt1, pcov1 = curve_fit(line, m_bright, y_bright)

print("a =", popt1[0], "+/-", pcov1[0,0]**0.5)
print("b =", popt1[1], "+/-", pcov1[1,1]**0.5)

x = np.linspace(m[1], max(m), 1000)


plt.plot(m,y,'o',ms=8,mew=2, mfc='none',color='red',label='Galaxies')
plt.plot(x, line(x, popt1[0], popt1[1]), lw=2, color='black',label='Best fit for brightest')
plt.xlabel(r'$\mathrm{m}$')
plt.ylabel(r'$\mathrm{log(N(m))}$')
#plt.xlim(10,18)
plt.title(r'$\mathrm{a = %.3f \pm %.3f, b = %.3f \pm %.3f}$' % (popt1[0],pcov1[0,0]**0.5,popt1[1],pcov1[1,1]**0.5),y=1.02)
plt.legend()
#plt.savefig('small-patch-number-count-plot.png',bbox_inches='tight')
plt.show()

#%%

# =============================================================================
# 
# =============================================================================

plt.hist(m,cumul=True,bins=20)
plt.yscale('Log')

#%%

# =============================================================================
# SECOND PART
# =============================================================================

# =============================================================================
# WITH ERRORS
# =============================================================================

# Read in data from catalogue text file
id_no,x,y,flux,flux_error,bkg,bkg_error,m,m_error = np.loadtxt('catalogue_patch_14_galaxies_error.txt', delimiter="\t", unpack=True)


# Calculate N(m) for each value of m
N_m = []

for i in range(len(m)):
    mag = m[i]
    
    count = 0
    for j in range(len(m)):
        
        # Add 1 to the count if a magnitude is less (brighter) than mag
        if m[j] < mag:
            count += 1
            
    N_m.append(count)



#%%

# =============================================================================
# PLOT SCATTER TO VIEW
# =============================================================================

y = np.log10(N_m)

yerr = 0.6 * m_error

xerr = m_error

plt.errorbar(m,y,yerr,xerr,'x',mew=2,ms=10,color='red')
plt.xlabel('m')
plt.ylabel('log(N(m))')
#plt.xlim(9,20)
#plt.ylim(0,4)
#plt.title('Small patch number count plot',y=1.02)
plt.show()

#%%

# =============================================================================
# BEST FIT LINE
# =============================================================================

def line(x, a, b):
    return a * x + b

y = np.log10(N_m)

mask = np.isfinite(y)
y = y[mask]

m = m[mask]

popt, pcov = curve_fit(line, m, y)

print("a =", popt[0], "+/-", pcov[0,0]**0.5)
print("b =", popt[1], "+/-", pcov[1,1]**0.5)

x = np.linspace(m[1], max(m), 1000)


plt.plot(m,y,'o',ms=8,mew=2, mfc='none',color='red',label='Galaxies')
plt.plot(x, line(x, popt[0], popt[1]), lw=2, color='black',label='Best fit line')
plt.xlabel(r'$\mathrm{m}$')
plt.ylabel(r'$\mathrm{log(N(m))}$')
#plt.xlim(10,18)
plt.title(r'$\mathrm{a = %.3f \pm %.3f, b = %.3f \pm %.3f}$' % (popt[0],pcov[0,0]**0.5,popt[1],pcov[1,1]**0.5),y=1.02)
plt.legend()
#plt.savefig('small-patch-number-count-plot.png',bbox_inches='tight')
plt.show()

#%%

N_m = np.asarray(N_m)

m_bright = np.where(m < 11.5)

bright_values = list(zip(*np.where(m <= 11.5)))

m_bright = (m[bright_values]).flatten()

y_bright = (y[bright_values]).flatten()

popt1, pcov1 = curve_fit(line, m_bright, y_bright)

print("a =", popt1[0], "+/-", pcov1[0,0]**0.5)
print("b =", popt1[1], "+/-", pcov1[1,1]**0.5)

x = np.linspace(m[1], max(m), 1000)


plt.plot(m,y,'o',ms=8,mew=2, mfc='none',color='red',label='Galaxies')
plt.plot(x, line(x, popt1[0], popt1[1]), lw=2, color='black',label='Best fit for brightest')
plt.xlabel(r'$\mathrm{m}$')
plt.ylabel(r'$\mathrm{log(N(m))}$')
#plt.xlim(10,18)
plt.title(r'$\mathrm{a = %.3f \pm %.3f, b = %.3f \pm %.3f}$' % (popt1[0],pcov1[0,0]**0.5,popt1[1],pcov1[1,1]**0.5),y=1.02)
plt.legend()
#plt.savefig('small-patch-number-count-plot.png',bbox_inches='tight')
plt.show()


#%%

# =============================================================================
# SOMETHING FROM NOAM'S ANALYSIS
# =============================================================================

# Create the histogram
values, base = np.histogram(m,bins=14)
# Create the cumulative histogram
cumul = np.cumsum(values)
err = np.sqrt(cumul)
# Take the log10 of N(m)
y_cumul = np.log10(cumul)
err = np.log10(err)
# Plot the number count
plt.scatter(base[:-1], y_cumul, color='none', edgecolor='blue')
#plt.errorbar(base[:-1],y_cumul,yerr=err,capsize=2,ecolor='blue',fmt=' ')

#%%
# Attempt to fit the data
def f(m,a,c):
    return a*m + c

base = base[13:-9]
y_cumul = y_cumul[13:-9]

popt,pcov = curve_fit(f,base[:-1],y_cumul)

x = np.linspace(base[0],base[-1],1000)

plt.plot(x,f(x,*popt),color='red',lw=2,label='Line fit')
plt.plot(m,y_cumul[:-1])

print(popt)
plt.xlabel('magnitude')
plt.ylabel('log$_{10}$(N(m))')