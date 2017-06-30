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

cat = fits.getdata(input_cats[2])
'''
for c in input_cats:
    # read in sextractor catalog
    # the catalog is in extension 2
    cat = fits.getdata(c,2)
    plt.figure()
    
''' 
