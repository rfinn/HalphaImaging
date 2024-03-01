#!/usr/bin/env python

'''
Some INT images have missing header info

the most important fields are CRVAL1 and CRVAL2

going to write this so you can enter a ref image and then update the header of the
broken image using the coordinates in the ref image.


USAGE:

python ~/github/Virgo/programs/INT_fixheader.py --ref r1442821.fit --image r1442823.fit --fixall

'''


from astropy.io import fits
from astropy import wcs
import argparse
import os
import sys

parser = argparse.ArgumentParser(description ='Create a crazy big catalog from HL, AGC, NSA')
parser.add_argument('--image',dest = 'image', default=None,help='image with bad header')
parser.add_argument('--ref',dest = 'ref', default=None,help='reference image, to use RA and DEC from')
#parser.add_argument('--fields',dest = 'fields', nargs='+',default=['CRVAL1','CRVAL2','AIRMASS'],help='list of header fields to update. default is CRVAL1 and CRVAL2')

# even for MEF files, the main header is 0, and this is what comes in 
#parser.add_argument('--wcs',dest = 'wcs', action='store_true', default=False, help='set this if you need to update basic wcs fields (ra,dec,equinox,crval1,crval2,cd1_1,cd2_2)')
parser.add_argument('--fixall',dest = 'fixall', action='store_true', default=False, help='set this to copy the entire header')
#parser.add_argument('--postsplit',dest = 'postsplit', action='store_true', default=False, help='set this if you are fixing headers after the mef fits file was split into 4 separate files.')

args = parser.parse_args()


# read in image and header for the image that needs to be updated
hdu = fits.open(args.image)
hdu.verify()
## adding these two statements to try to get around FITS card error with CD1_1
#w = wcs.WCS(hdu[0].header)
#h.verify()

# get header from reference image
hdu_ref = fits.open(args.ref)
href = hdu_ref[0].header

if args.fixall:
    # find fields in reference header that are not in bad header
    # for INT data, the telescope telemetry fields are missing (rather than empty)
    # this may not work for other datasets...
    goodh = set(href)
    badh = set(hdu[0].header)
    fields = list(goodh.difference(badh))
    print('fields to update: ',fields)
else:
    fields = args.fields


# testing
sys.exit()
for f in fields:
    #print(f)
    if (f.find('DUMMY') > -1) | (f.find('HISTORY') > -1):
        continue
    try: 
        newval = href[f]
    except:
        print('error with ',f)
        continue

    hdu[0].header.set(f,newval)

    ## the following commands are an alternate way to update the header using sethead
    ## however, when I ran this on the virgo vms, it creates a new image *.fit,1
    ## and neither image seemed to have the updated header values
    ## not sure what's going on with that, so I chose to try to fix the issue with CD1_1
    ## because astropy has more documentation than wcstools
    '''
    try:
        os.system('sethead '+args.image+' '+f+'='+str(newval))
    except KeyError:
        os.system('sethead -k '+args.image+' '+f+'='+str(newval))
    '''
hdu.writeto(args.image,overwrite=True,output_verify='ignore')
