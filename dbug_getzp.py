#!/usr/bin/env python

'''

'''
from astropy.io import ascii
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

# read in positions and mags of Landolt standards
stand = ascii.read('mystandards.csv', data_start=1, delimiter=',')
sfile = stand['col1']
landolt_rmag = stand['col5']
ximage = np.array(stand['col6'],'f')
yimage = np.array(stand['col7'],'f')
zp = stand['col8']


# read in SE cats
secats = ['nshdztrc7131t0019o00.cat','nshdztrc7131t0020o00.cat','nshdztrc7131t0073o00.cat','nshdztrc7131t0074o00.cat']
# read in one catalogs

cat = fits.getdata(secats[0],2)
    
# create array with nrow = nstandards, and format like SE cat
newarray = np.zeros(len(stand),dtype = cat.dtype)

# loop through SE catalogs

for c in secats:
    cat = fits.getdata(c,2)
    # for appropriate standards
    prefix = c.split('.')[0]
    for i,s in enumerate(sfile):
        if s.find(prefix) > -1:
            # find location of landolt standards in SE cats
            distance = np.sqrt((ximage[i] - cat['X_IMAGE'])**2 + (yimage[i] - cat['Y_IMAGE'])**2)
            imatch = np.where(distance == min(distance))
            # populate new array with matched lines in SE cat
            newarray[i] = cat[imatch]
                              


# plot measured mag vs known mag
testmag = newarray['MAG_PETRO']
plt.figure()
#plt.plot(landolt_rmag, zp+newarray['MAG_APER'][:,4]-landolt_rmag,'bo')
#plt.plot(landolt_rmag, zp+newarray['MAG_PETRO']-landolt_rmag,'bo')
residual = zp + testmag - landolt_rmag
plt.plot(landolt_rmag[0:5], residual[0:5],'bo',label='PG0918+029')
plt.plot(landolt_rmag[5:13], residual[5:13],'go',label='RU_149A')
plt.plot(landolt_rmag[13:16], residual[13:16],'co',label='PG1528')
plt.plot(landolt_rmag[16:-1], residual[16:-1],'mo',label='PG1633')
plt.legend()
plt.xlabel('Landolt R mag',fontsize=14)
plt.ylabel('Measured mag',fontsize=14)
