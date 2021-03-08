# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 13:39:47 2021

@author: nomou
"""

import numpy as np
import matplotlib.pyplot as plt
import A1_tools as a1

#%%
img = a1.load_img('mosaic.fits')

nrows = len(img)
ncols = len(img[0])
N = nrows*ncols

#%%
fig_og=plt.figure()
ax_og=fig_og.add_subplot()
og = ax_og.pcolormesh(img)
cb_og = fig_og.colorbar(og)
#%%
img2 = np.copy(img)

#%%

# Central bleed
img2[0:4611,1427:1447] = 0
img2[3136:3296,1384:1500] = 0
img2[2950:3475,1420:1454] = 0
# Borders
img2[0:119,0:2570] = 0
img2[0:4611,0:117] = 0
img2[4514:4611,0:2570] = 0
img2[0:4611,2471:2570] = 0
# Horizontal bleed, bottom
img2[217:228,1298:1538] = 0
img2[119:141,1289:1538] = 0
img2[231:249,1389:1477] = 0
img2[313:324,1026:1704] = 0
img2[313:358,1639:1650] = 0
img2[422:439,1101:1652] = 0
img2[422:455,1025:1046] = 0
img2[139:172,1406:1428] = 0
img2[139:151,1446:1456] = 0
img2[227:269,1416:1457] = 0
img2[323:334,1359:1503] = 0
img2[323:360,1445:1454] = 0
img2[438:452,1370:1428] = 0
img2[438:453,1446:1486] = 0
img2[333:361,1396:1428] = 0
img2[333:361,1442:1477] = 0
img2[359:390,1446:1456] = 0
img2[449:476,1409:1460] = 0
# Stars
img2[4079:4119,547:573] = 0
img2[3292:3346,756:794] = 0
img2[3202:3419,771:780] = 0
img2[2748:2791,955:991] = 0
img2[2704:2835,967:979] = 0
img2[2221:2359,887:922] = 0
img2[4016:4050,1446:1467]=0
img2[3707:3803,2120:2149]=0
img2[3383:3445,2454:2573]=0
img2[2281:2339,2116:2148]=0
img2[1399:1455,2085:2093]=0
img2[1413:1437,2076:2102]=0
# Middle star
img2[2862:3519,1113:1745] = 0
# Star new
img2[4065:4142,518:601] = 0
img2[3206:3433,680:889] = 0
img2[2699:2848,898:1055] = 0
img2[2216:2370,831:981] = 0
img2[3704:3821,2084:2192] = 0
img2[2272:2355,2092:2178] = 0
# Corners, added on 16/02/2021
img2[4400:4518,115:2475] = 0
img2[112:4518,115:300]   = 0
img2[112:4518,2364:2475] = 0
img2[112:205,109:2475]   = 0

#%%
fig2 = plt.figure()
ax2 = fig2.add_subplot()
plt2 = ax2.pcolormesh(img2)
cb2 = fig2.colorbar(plt2)
#%% Save the array amiright
np.savetxt("modified_image.txt",img2)