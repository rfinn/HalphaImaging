#!/usr/bin/env python

'''
PURPOSE:

The goal of the program is to call uat_mask.py for all files that
match the input file string.

USAGE:

from within ipython type:

   %run ~/github/HalphaImaging/uat_mask.py --image 'A1367-140231-R.fits'

you just need to run this on R-band images.


PROCEDURE:


REQUIRED MODULES:
   os
   astropy
   numpy
   argsparse
   matplotlib
   scipy


'''

import os
import subprocess
import argparse
import glob




parser = argparse.ArgumentParser(description ='Run uat_mask.py on all images that match input string')
parser.add_argument('--string', dest = 'string', default = None, help = 'image string to match (e.g. A1367)')
parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/', help = 'Directory that contains uat_mask.py')
args = parser.parse_args()

search_string = args.string+'*-R.fits'
print search_string
input_images = glob.glob(search_string)

for image in input_images:
    print 'masking image ',image
    os.system('~/github/HalphaImaging/uat_mask.py --image '+image)
