#!/usr/bin/env python


'''

GOAL: create coadd images 

OVERVIEW:
* designed to complete final stages of coaddition for INT WFC images
* assumes processing is complete through astrometric correction.  
  - code will look for .head file

REFERENCES:
https://photutils.readthedocs.io/en/stable/background.html

https://reproject.readthedocs.io/en/stable/mosaicking.html

ccdproc wcs projection
https://ccdproc.readthedocs.io/en/latest/image_combination.html
'''

#from photutils import make_source_mask
import os
homedir = os.getenv("HOME")
import sys
print(homedir)

sys.path.append('/home/rfinn/github/HalphaImaging/python3/')
import imutils
import ccdproc as ccdp

from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp

from astropy.nddata import CCDData

##########################################################
### SUBTRACT MEDIAN FROM IMAGE
##########################################################

def subtract_median(ic):
    print('subtracting median from images')
    for hdu, fname in ic.hdus(return_fname=True):    
        # background subtraction
        hdu.data = imutils.subtract_median_sky(hdu.data)
        hdu.writeto('m'+fname,overwrite=True)
def reproject_images(ic):
    ##########################################################
    ### FIND BEST WCS FOR ALL IMAGES
    ##########################################################

    # find wcs for all images
    print('finding optimal output wcs')
    wcs_out, shape_out = find_optimal_celestial_wcs([h for h in ic.hdus()])


    ##########################################################
    ### REPROJECTING IMAGES
    ##########################################################

    print('Starting loop to reproject images...')

    for hdu, imagename in ic.hdus(return_fname=True):    
        print('processing ',imagename)

        print('subtracting median...')
        hdu.data = imutils.subtract_median_sky(hdu.data)

        # change to reproject_exact once code works
        print('reprojecting...')        
        reprodata, footprint = reproject_interp(hdu, wcs_out,shape_out=shape_out)
        exposuremask = ~footprint.astype('bool')
        newimg = CCDData(reprodata, wcs=wcs_out, unit='adu',mask=exposuremask)
        print('saving ','rm'+imagename)
        newimg.write('rm'+imagename,overwrite=True)
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ='Run gui for analyzing Halpha images')

    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC')
    args = parser.parse_args()

    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']

    ic = ccdp.ImageFileCollection(os.getcwd(), keywords=keys, glob_include=args.filestring+'*.fits',glob_exclude='*coadd*.fits')
    reproject_images(ic)
