#!/usr/bin/env python

'''
PURPOSE:

The goal of the program is to call uat_mask.py for all
R-band images with a prefix that matches the input string.  

USAGE:

~/github/HalphaImaging/uat_mask_all.py --string 'A1367'

This will find all images that match A1367*-R.fits and
run uat_mask.py on each image.

We use the same mask for the R and Halpha images, so you
just need to run this on R-band images.


PROCEDURE:

- construct match string from input prefix
- get list of files that match string
- call uat_mask.py for each image


REQUIRED MODULES:

os
argsparse
glob

'''

import os
import argparse
import glob
import sys

parser = argparse.ArgumentParser(description ='Run uat_mask.py on all images that match input string')
parser.add_argument('--cluster', dest = 'string', default = None, help = 'image string to match (e.g. A1367)')
parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/', help = 'Directory that contains uat_mask.py')
args = parser.parse_args()

search_string = args.string+'*-R.fits'
print search_string
input_images = glob.glob(search_string)

for image in input_images:
    print 'masking image ',image
    try:
        os.system(args.d+'uat_mask.py --image '+image)
    except SystemExit:
        sys.exit()
        
