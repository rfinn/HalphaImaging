#!/usr/bin/env python

'''
GOAL: 
- align r and H-alpha images that are created by theli

USAGE:

%run ~/github/HalphaImaging/int_align_images.py --image1 coadd-Ha.fits --image2 coadd-r.fits --weight2 coadd-r.weight.fits


'''


import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import CCDData
import glob
import time
# an astropy module to reproject images
from reproject import reproject_interp
from reproject import reproject_exact
from ccdproc import wcs_project

import argparse

start_time = time.perf_counter()
parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--image1', dest = 'image1', default = 'test-ha.fits', help = 'Image to serve as reference')
parser.add_argument('--image2', dest = 'image2', default = 'test-r.fits', help = 'Image to align to reference')
parser.add_argument('--weight2', dest = 'weight2', default = None, help = 'Weight map for image2')

args = parser.parse_args()

hdu1 = fits.open(args.image1)
hdu2 = fits.open(args.image2)
#hdu1 = CCDData.read(args.image1,unit='adu/s')
#hdu2 = CCDData.read(args.image2,unit='adu/s')

if args.weight2 is not None:
    hdu2w = fits.open(args.weight2)
    #hdu2w = CCDData.read(args.weight2,unit='adu')

print('\t shifting image')
im2new, im2footprint = reproject_interp(hdu2[0], hdu1[0].header,shape_out=hdu1[0].data.shape)
#wcsout = WCS(hdu1)
#im2new = wcs_project(hdu2,target_wcs=wcsout,target_shape=hdu1.data.shape)
#im2new.write(args.image2.split('.fits')[0]+'-shifted.fits', overwrite=True)
#hdu3 = fits.open(args.image2.split('.fits')[0]+'-shifted.fits')
newheader = hdu2[0].header
# update wcs to image 1
wcskeys = ['NAXIS1','NAXIS2','CRVAL1','CRVAL2','CRPIX1','CD1_1','CD1_2','CRPIX2','CD2_1','CD2_2']
for k in wcskeys:
    newheader.set(k, value=hdu1[0].header[k])

newheader.set('HAIMAGE',args.image1)
fits.writeto(args.image2.split('.fits')[0]+'-shifted.fits', im2new, header=newheader, overwrite=True)
hdu2.close()
#hdu3.close()
end_time = time.perf_counter()
print('\t elapsed time = ',end_time - start_time)

# fix header keyword
if args.weight2 is not None:
    print('\t shifting weight image')
    im2wnew, im2wfootprint = reproject_interp(hdu2w[0], hdu1[0].header)
    newheader = hdu2w[0].header
    # update wcs to image 1
    wcskeys = ['NAXIS1','NAXIS2','CRVAL1','CRVAL2','CRPIX1','CD1_1','CD1_2','CRPIX2','CD2_1','CD2_2']
    for k in wcskeys:
        newheader.set(k, value=hdu1[0].header[k])
    
    #im2wnew = wcs_project(hdu2w[0],WCS(hdu1[0].header),target_shape=hdu1[0].data.shape)
    fits.writeto(args.weight2.split('.fits')[0]+'-shifted.fits', im2wnew,newheader, overwrite=True)
    hdu2w.close()
hdu1.close()
end_time = time.perf_counter()
print('\t total time = ',end_time - start_time)
