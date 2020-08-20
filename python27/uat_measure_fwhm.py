#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import argparse
import glob
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description ='This program will measure the FWHM for a bunch of images using the SExtractor catalog.')
parser.add_argument('--filestring', dest = 'filestring', default = 'hcftr*o00.cat', help = 'string to use to get input files (default = "hcftr*o00.cat")')

args = parser.parse_args()

input_cats = glob.glob(args.filestring)
input_cats.sort()
#cat = fits.getdata(input_cats[2])
print 'Image         mean (median)'
for c in input_cats:
    # read in sextractor catalog
    # the catalog is in extension 2
    cat = fits.getdata(c,2)
    # select galaxies that are likely to be stars
    # limit magnitude range so that stars are not saturated and not too faint
    flag = (cat.MAG_AUTO > 10) &  (cat.FLAGS == 0) 
    #flag = np.ones(len(cat.MAG_AUTO),'bool')
    print '%20s: %5.1f (%5.1f)'%(c,np.mean(cat.FWHM_IMAGE[flag]),np.median(cat.FWHM_IMAGE[flag]))


'''
imstat hcftr7893t0061o00.fits nclip=3

daofind hcftr7893t0062o00.fits

psfmeasure hcftr7893t0062o00.fits imagecur =  hcftr7893t0062o00.fits.coo.2

'''



