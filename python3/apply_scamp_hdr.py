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
### UPDATE WCS HEADER IF .HEAD FILE EXISTS
##########################################################


def update_header(ic):

    # update image headers if .head file exists
    print('updating headers if .head file exists')
    for hdu, fname in ic.hdus(return_fname=True):
        filebasename = fname.split('.fits')[0]
        altheader = filebasename+'.head'
        # read in scamp header if it exists
        if os.path.exists(altheader):
            # build header from .head filt
            wcsinfo = imutils.convert_headfile_header(altheader)
            for key in wcsinfo:
                # skip some keys that were causing problems (not sure why or what they do)
                if key.startswith('PV'):
                    continue
                hdu.header.set(key,value=wcsinfo[key])
                # scamp put some weird numbers in here from phot calib process
                # resetting to more reasonable values
            hdu.header['SATURATE']=40000
            hdu.header['FLXSCALE']=1
            hdu.writeto(fname,overwrite=True)
        else:
            continue



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ='Run gui for analyzing Halpha images')

    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC')
    args = parser.parse_args()

    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']

    ic = ccdp.ImageFileCollection(os.getcwd(), keywords=keys, glob_include=args.filestring+'*.fits',glob_exclude='*coadd*.fits')
    update_header(ic)
