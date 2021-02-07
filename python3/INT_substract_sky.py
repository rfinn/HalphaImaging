#!/usr/bin/env python

'''
Read in INT image (chip 4 only) and subtract median or mode

'''

from astropy.io import fits
import numpy as np
import argparse

parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--image', dest = 'image', default = 'test.fits', help = 'Image to subtract sky from')

args = parser.parse_args()

a = fits.getdata(args.image)
h = fits.getheader(args.image)
skysub = a - np.median(a)
outfile = args.image.split('.fits')[0]+"-sky.fits"
fits.writeto(outfile,skysub,header=h)
