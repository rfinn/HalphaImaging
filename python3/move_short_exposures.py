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


##########################################################
### GET RID OF SHORT EXPOSURES
##########################################################
def move_short_exposures(ic):
    # get rid of 15 sec r-band exposure, b/c can't figure out how to do the weighting!!!
    if not os.path.exists('shortexposure'):
        os.mkdir('shortexposure')
    for i,f in enumerate(ic.files):
        if ic.summary['exptime'][i] < 20.:
            os.rename(f,'shortexposure/'+f)
            froot = f.split('.fits')[0]
            os.system('mv '+froot+'* shortexposure/.')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ='Run gui for analyzing Halpha images')

    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC')
    args = parser.parse_args()

    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']

    ic = ccdp.ImageFileCollection(os.getcwd(), keywords=keys, glob_include=args.filestring+'*.fits',glob_exclude='*coadd*.fits')
    move_short_exposures(ic)
