#!/usr/bin/env python

'''
GOAL


PROCEDURE


REQUIRED MODULES
ccdproc


'''

import glob
import os
from astropy.io import fits
import argparse

from astropy.nddata import CCDData
import ccdproc
from astropy.modeling import models

parser = argparse.ArgumentParser(description ='Subtract overscan bias and trim images.')
parser.add_argument('--filestring', dest = 'filestring', default = 'c', help = 'string to use to get input files (default = "c", which grabs all of the files "c*o00.fits")')
parser.add_argument('--irafbiassec', dest = 'irafbiassec', default= '[4100:4150, 1:4150]', help = 'biassec in iraf notation.  default is [4100:4150, 1:4150], which applies to HDI camera')
parser.add_argument('--iraftrimsec', dest = 'iraftrimsec', default= '[1:4095, 1:4109]', help = 'biassec in iraf notation.  default is [4100:4150, 1:4150], which applies to HDI camera')
#parser.add_argument('--gain', dest = 'gain', default= 1.3, help = 'gain in e-/ADU.  default is 1.3, which applies to HDI camera')
#parser.add_argument('--rdnoise', dest = 'rdnoise', default= 7.3, help = 'gain in e-/ADU.  default is 1.3, which applies to HDI camera')


args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)

poly_model = models.Polynomial1D(1)

for f in files:
    # read in image
    # was having trouble getting image into the format that ccdproc wants
    with fits.open(f) as hdu1:
        print 'working on ',f
        # convert data to CCDData format and save header
        ccd = CCDData(hdu1[1].data, unit='adu')
        header = hdu1[0].header

        # subtract overscan
        o_subtracted = ccdproc.subtract_overscan(ccd, fits_section = args.irafbiassec, model=poly_model)
        header['HISTORY'] = 'overscan subtracted '+args.irafbiassec

        # trim image
        trimmed = ccdproc.trim_image(o_subtracted, fits_section = args.iraftrimsec)
        header['HISTORY'] = 'trimmed '+args.iraftrimsec
        
        fits.writeto('tr'+f,trimmed,header, overwrite='True')
        header.append(card=('CCDSEC',args.iraftrimsec,'CCD SECTION'))
#        # remove cosmic rays
#        crmask, crimage = ccdproc.cosmicray_lacosmic(trimmed, gain = float(args.gain), readnoise = float(args.rdnoise))
#        header['HISTORY'] = 'Cosmic rays rejected using ccdproc.cosmicray_lacosmic '
        
        # save output as t*orgininal_name

#        fits.writeto('ztr'+f,crimage,header)
        

        # close images if necessary - happens automatically

        hdu1.close()
