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
from astropy.io import fits


parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--filename', dest = 'filename', default = 'iwf000.fits', help = 'filename of input file')
parser.add_argument('--badval', dest = 'badval', type=int, default = 0, help = 'set value to use when identifying bad pixels.  default is zero.')

args = parser.parse_args()


# get image name from the argument line
filename = args

# read in image
hdu = fits.open(args.filename)

# create a flag with values == args.badval set to zero
weight = hdu[0].data == int(badval)
weight = ~weight

# save the flag as a weight image
outfile = filename.replace('.fits','.weight.fits')

fits.writeto(outfile, weight, header=hdu[0].header, overwrite=True)



