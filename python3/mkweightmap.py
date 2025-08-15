#!/usr/bin/env python
"""
GOAL:
- in Becky's data, the chip gaps have a value of zero, but there is no associated weight map
- the goal of this program is to create a weight map from a mosaic mosaiced image

PROCEDURE:
- get image name from the argument line
- read in image
- create a flag with values > 0 set to true
- save the flag as a weight image

USAGE:

python ~/github/HalphaImaging/python3/uat_MOSmkweightmap.py

"""
from astropy.io import fits
import sys
import argparse
import numpy as np
import os


parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--filename', dest = 'filename', default = 'iwf000.fits', help = 'filename of input file')
parser.add_argument('--badval', dest = 'badval', type=int, default = 0, help = 'set value to use when identifying bad pixels.  default is zero.')

args = parser.parse_args()


# get image name from the argument line
filename = args.filename

# read in image
hdu = fits.open(args.filename)

# create a flag with values == args.badval set to zero
weight = hdu[0].data == int(args.badval)

# now the bad values are set to zero and good values are equal to 1
weight = np.array(~weight,dtype=np.int8)

# check if there is a bpm for mosaic
mfile = args.filename.replace('.fits','_bpm.fits')
if os.path.exists(mfile):
    print(f"found a mosaic pbm file to add to mask: {mfile}")
    # open mask
    mhdu = fits.open(mfile)

    # create array with good values equal to 1
    mosaic_weight = mhdu[0].data == 0
    
    # convert to integer array
    mosaic_weight = np.array(mosaic_weight, dtype=np.int8)

    # add mosaic mask to weight image
    weight = weight * mosaic_weight

    # close the mask file
    mhdu.close()


# save the flag as a weight image
outfile = filename.replace('.fits','.weight.fits')

# writing the output
fits.writeto(outfile, weight, header=hdu[0].header, overwrite=True)

'''
Commenting this out because the right way to do this is to convolve the mask 
in the same way the image is convolved.


if outfile.startswith('if'):
    if os.path.exists('g'+filename):
        # copy a version
        fits.writeto('g'+outfile, weight, header=hdu[0].header, overwrite=True)
'''


