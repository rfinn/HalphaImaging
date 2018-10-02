#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import glob

# an astropy module to reproject images
from reproject import reproject_interp


import argparse

parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--image1', dest = 'image1', default = 'test-ha.fits', help = 'Image to serve as reference')
parser.add_argument('--image2', dest = 'image2', default = 'test-r.fits', help = 'Image to align to reference')

args = parser.parse_args()

hdu1 = fits.open(args.image1)[0]
hdu2 = fits.open(args.image2)[0]

im2new, im2footprint = reproject_interp(hdu2, hdu1.header)

fits.writeto(args.image2.split('.fits')[0]+'-shifted.fits', im2new, hdu1.header, overwrite=True)

#hdu1.close()
#hdu2.close()
