#!/usr/bin/env python

'''
GOAL: 
- align r and H-alpha images that are created by theli

USAGE:

%run ~/github/HalphaImaging/int_align_images.py --image1 coadd-Ha.fits --image2 coadd-r.fits --weight2 coadd-r.weight.fits

can run on multiple in parallel 

parallel --eta --ungroup python ~/github/HalphaImaging/python3/MOS_align_images.py --image1 NGC5846_0{1}Ha.fits --image2 NGC5846_0{2}R.fits ::: {1..7} ::: {1..7}

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
#not using this one
#from ccdproc import wcs_project

# this should find the best wcs, and by default, N is up, E to left - yay!
from reproject.mosaicking import find_optimal_celestial_wcs

import argparse



start_time = time.perf_counter()
parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--image1', dest = 'image1', default = 'test-ha.fits', help = 'Image to serve as reference')
# mosaic images don't have weight images, so just saving the footprint as the weight image
#parser.add_argument('--weight1', dest = 'weight1', default = None, help = 'Weight map for image1')
#parser.add_argument('--image2', dest = 'image2', default = 'test-r.fits', help = 'Image to align to reference')
#parser.add_argument('--weight2', dest = 'weight2', default = None, help = 'Weight map for image2')

args = parser.parse_args()

image2 = args.image1.replace('Ha.fits','R.fits')
hdu1 = fits.open(args.image1)
hdu2 = fits.open(image2)
#hdu1 = CCDData.read(args.image1,unit='adu/s')
#hdu2 = CCDData.read(args.image2,unit='adu/s')

#if args.weight1 is not None:
#    hdu1w = fits.open(args.weight1)
    #hdu2w = CCDData.read(args.weight2,unit='adu')
#if args.weight2 is not None:
#    hdu2w = fits.open(args.weight2)
    #hdu2w = CCDData.read(args.weight2,unit='adu')

# Find optimal shape of output image
wcs_out, shape_out = find_optimal_celestial_wcs([hdu1,hdu2])
header_out = wcs_out.to_header()


# construct output filenames
ra = header_out['CRVAL1']
dec = header_out['CRVAL2']

dateobs = hdu1[0].header['DATE-OBS'].split('T')[0].replace('-','')
object = args.image1[0:10]

rimage_outname = f'VF-{ra:.3f}+{dec:.3f}-MOS-{dateobs}-{object}-R.fits'
haimage_outname = f'VF-{ra:.3f}+{dec:.3f}-MOS-{dateobs}-{object}-Ha4.fits'

rweight_outname = f'VF-{ra:.3f}+{dec:.3f}-MOS-{dateobs}-{object}-R.weight.fits'
haweight_outname = f'VF-{ra:.3f}+{dec:.3f}-MOS-{dateobs}-{object}-Ha4.weight.fits'




##
# SHIFT Halpha IMAGE
##
print(f'shifting image {args.image1} -> {haimage_outname}' )
im2new, im2footprint = reproject_interp(hdu1[0], wcs_out, shape_out=shape_out)

# we want to keep most of the info in the original header, and just update the WCS
newheader = hdu1[0].header
# update wcs to image 1
wcskeys = ['NAXIS1','NAXIS2','CRVAL1','CRVAL2','CRPIX1','CD1_1','CD1_2','CRPIX2','CD2_1','CD2_2']
reprojectkeys = ['NAXIS1','NAXIS2','CRVAL1','CRVAL2','CRPIX1','PC1_1','PC1_2','CRPIX2','PC2_1','PC2_2']
for i,k in enumerate(wcskeys):
    if k == 'NAXIS1':
        newheader.set(k, value=shape_out[1])
    elif k == 'NAXIS2':
        newheader.set(k, value=shape_out[0])
    else:
        newheader.set(k, value=header_out[reprojectkeys[i]])

newheader.set('FILTER','ha4')

fits.writeto(haimage_outname, im2new, header=newheader, overwrite=True)
hdu1.close()

fits.writeto(haweight_outname, im2footprint,newheader, overwrite=True)
hdu1w.close()


##
# SHIFT R IMAGE
##
print()
print(f'shifting image {image2} -> {rimage_outname}' )

im2new, im2footprint = reproject_interp(hdu2[0], wcs_out, shape_out=shape_out)

# we want to keep most of the info in the original header, and just update the WCS
newheader = hdu2[0].header
# update wcs to image 1
wcskeys = ['NAXIS1','NAXIS2','CRVAL1','CRVAL2','CRPIX1','CD1_1','CD1_2','CRPIX2','CD2_1','CD2_2']
for k in wcskeys:
    if k == 'NAXIS1':
        newheader.set(k, value=shape_out[1])
    elif k == 'NAXIS2':
        newheader.set(k, value=shape_out[0])
    else:
        newheader.set(k, value=header_out[k])

newheader.set('HAIMAGE',args.image1)
newheader.set('FILTER','R')

fits.writeto(rimage_outname, im2new, header=newheader, overwrite=True)
hdu2.close()

fits.writeto(rweight_outname, im2wnew,newheader, overwrite=True)
hdu2w.close()




end_time = time.perf_counter()
print('\t elapsed time = ',end_time - start_time)

