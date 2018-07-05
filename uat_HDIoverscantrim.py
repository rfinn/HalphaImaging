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


parser = argparse.ArgumentParser(description ='Run sextractor, scamp, and swarp to determine WCS solution and make mosaics')
parser.add_argument('--filestring', dest = 'filestring', default = 'c', help = 'string to use to get input files (default = "c", which grabs all of the files "c*o00.fits")')

args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)



