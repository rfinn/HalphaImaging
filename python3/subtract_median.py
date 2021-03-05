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
import glob
from astropy.io import fits



##########################################################
### SUBTRACT MEDIAN FROM IMAGE
##########################################################

def subtract_median(files,overwrite=False):
    print('subtracting median from images')
    for fname in files:
        if not overwrite:
            print("{} -> m{}".format(fname,fname))
        else:
            print("{} -> {}".format(fname,fname))
        if os.path.exists("m"+fname) and not overwrite:
            print("m"+fname,' already exists.  moving to next file')
            continue
        
        hdu = fits.open(fname)

        # background subtraction
        hdu[0].data,median = imutils.subtract_median_sky(hdu[0].data)
        hdu[0].header.set('MEDSUB',value=median,comment='median subtraction')
        if overwrite:
            hdu.writeto(fname,overwrite=True)
        else:
            hdu.writeto("m"+fname,overwrite=True)
        hdu.close()
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ='Subtract the median from images, after masking out objects and growing mask.')

    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC')
    parser.add_argument('--overwrite', action = 'store_true', default = False, help = 'overwrite file?  the default is false, so that a new file with m prefix is created.')    
    args = parser.parse_args()

    #if args.hdi:
    #    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']
    #else:
    #    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']

    files = glob.glob(args.filestring+'*.fits')
    files.sort()
    print(files)
    subtract_median(files,overwrite=args.overwrite)
