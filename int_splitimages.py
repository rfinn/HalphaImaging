#!/usr/bin/env python

'''
Read in INT multiextension WFC images and split them into 4 images

'''

from astropy.io import fits
import glob
import os


files = glob.glob('r*.fit')
if os.path.exists('ORIGINALS'):
    print('directory ORIGINALS already exists')
else:
    os.mkdir('ORIGINALS')
for f in files:
    filestring ,t = f.split('.')
    a = fits.open(f)
    fits.writeto(filestring+'_1.fits',a[1].data,header=a[0].header+a[1].header)
    fits.writeto(filestring+'_2.fits',a[2].data,header=a[0].header+a[2].header)
    fits.writeto(filestring+'_3.fits',a[3].data,header=a[0].header+a[3].header)
    fits.writeto(filestring+'_4.fits',a[4].data,header=a[0].header+a[4].header)
    a.close()
    os.rename(f,'ORIGINALS/'+f)
    
