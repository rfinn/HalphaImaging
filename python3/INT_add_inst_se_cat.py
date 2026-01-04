#!/usr/bin/env python

'''
adding INSTRMNT keyword to header to all WFC*PA.cat files for INT data so scamp can fit the geometric distortion for each chip separately.

for some reason, sextractor won't add this into the header

USAGE:

python ~/github/Virgo/programs/INT_add_inst_se_cat.py 

'''


from astropy.io import fits
from astropy import wcs
import glob
import os

matchstrings = ['WFC*1PA.cat','WFC*2PA.cat','WFC*3PA.cat','WFC*4PA.cat']
instruments = ['INTWFC01','INTWFC02','INTWFC03','INTWFC04']
print(f"\nWorking on directory {os.getcwd()}")
for i in range(len(matchstrings)):
    files = glob.glob(matchstrings[i])
    print('chip ',i+1,' updating ',len(files),' cats')
    j=0
    for f in files:
        #print(j)
        # read in image and header for the image that needs to be updated
        
        hdu = fits.open(f)
        #hdu.verify()

        # trying again after implementing two commands that might fix the issue with CD1_1
        hdu[2].header.set('INSTRMNT',instruments[i],'ccd number')

        # add the filter as well
        t = f.split('.')
        hdu[2].header.set('FILTER',t[1], 'photmetric filter')
        #if 'WFC.r' in f:
        #    hdu[2].header.set('FILTER','r', filter,'photmetric filter')
        #elif 'WFC.Ha6657' in f:
        #    hdu[2].header.set('FILTER','Ha6657', filter,'photmetric filter')
        #elif 'WFC.Halpha' in f:
        #    hdu[2].header.set('FILTER','H', filter,'photmetric filter')
        hdu.writeto(f,overwrite=True)#,output_verify='ignore')
        hdu.close()
        j+= 1
