#!/usr/bin/env python


'''
GOAL:
- add conversion from ADU/s to flux density F_nu (erg/s/cm^2/Hz)
- add conversion from ADU/s to flux nu F_nu (erg/s/cm^2)

APPROACH:
- get ZP, filter, redshift from image header
- allow user to override header info in case the image header doesn't have the necessary info
- look up central wavelength of dictionary
- convert from ADU/s to flux density
- convert from flux density to flux
- convert from flux to luminosity (ergs/s)
- convert to SFR if halpha filter
- add values to image header

INPUTS:
- image (with header containing filter, ZP, redshift)

OPTIONAL INPUTS:
- filter
- ZP
- redshift

OUTPUT:
- image with updated header containing:
   - Fnu - conversion to flux density
   - nuFnu - conversion to flux
   - SFR - conversion to SFR (if Halpha filter)

'''

import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy.io import fits
from astropy import units as u
from astropy import constants as c
import argparse
import sys

center_wavelength = {'R': 6513.54, 'r':6292.28, 'ha':6574.74, 'ha4': 6620.52, 'ha8':6654.19, 'ha12':6698.53, 'ha16':6730.72}

fzeroAB = 3631.*u.Jy # Jy

parser = argparse.ArgumentParser(description = 'Calculate conversion from ADU/s to flux density, flux and SFR (if Halpha filter).  The program will read info from the image header (header keywords = FILTER, PHOTZP, REDSHIFT, ZDIST).  The user has the option to provide the required information (filter, ZP, redshift) at the command line.')
parser.add_argument('--image',dest = 'image', default = None, help = 'input image')
parser.add_argument('--filter',dest = 'filter', default = None, help = 'Image filter.  Options are: R, r, ha, ha4, ha8, ha12, ha16')
parser.add_argument('--ZP', dest = 'ZP', default=None,help = 'flux zeropoint in AB system')
parser.add_argument('--redshift', dest = 'redshift', default = None, help = 'galaxy redshift - required to add luminosity (and SFR if H-alpha filter) conversion to image header.')
parser.add_argument('--redshift_keyword', dest = 'redshift_keyword', default = 'REDSHIFT', help = 'header keyword for redshift.  Usually this is either REDSHIFT or ZDIST.')

args = parser.parse_args()

## get ZP, filter, redshift from image header
## allow user to override header info in case the image header doesn't have the necessary info
## look up central wavelength of dictionary



# read in image
imdat, header = fits.getdata(args.image, header=True)

# check to see if user provided filter
if args.filter:
    imfilter = args.filter
else:
    try:
        imfilter = header['FILTER']
    except KeyError:
        print('could not read filter from header!')
        print('please provide filter at the command line.')
        sys.exit(1)


# check to see if user provided ZP
if args.ZP:
    ZP = args.ZP
else:
    try:
        ZP = header['PHOTZP']
    except KeyError:
        print('could not read PHOTZP from header!')
        print('please provide ZP at the command line.')
        sys.exit(1)

# check to see if user provided redshift
if args.redshift:
    redshift = float(args.redshift)
else:
    try:
        redshift = float(header[args.redshift_keyword])
    except KeyError:
        print('could not read redshift from header!')
        print('please provide redshift at the command line.')
        sys.exit(1)

        
dl = cosmo.luminosity_distance(redshift)
clambda = center_wavelength[imfilter]*u.Angstrom
nu = c.c/clambda


## convert from ADU/s to flux density
'''
mAB = ZP - 2.5log10(flux)
mAB = -2.5log10(fnu/3631)
ZP = -2.5log10(fnu/3621) + 2.5log10(flux) = 2.5log(3631*flux/fnu)
3631*flux/fnu = 10**(ZP/2.5)
fnu = 3631*flux*10**(-ZP/2.5)
''' 
flux_density = fzeroAB*10.**(-1*ZP/2.5)

## convert from flux density to flux
flux = flux_density*nu 

## convert to luminosity
lumin = flux*4*np.pi*dl**2
## add values to image header
header.set('Fnu',flux_density.cgs.value,comment='conversion from image counts to erg/s/cm^2/Hz')
header.set('nuFnu',flux.cgs.value,comment='conversion from image counts to erg/s/cm^2')
header.set('log10(nuLnu)',np.log10(lumin.cgs.value),comment='conversion from image counts to erg/s (nuLnu); WMAP9')

## convert to SFR if halpha filter
if (imfilter.find('ha') > -1) | (imfilter.find('197') > -1):
    # kennicutt & evans
    # L in units of erg/s (not what I have!!!)
    # Cx = 41.27
    # where
    # log M* (Msun/yr) = log Lx - log Cx
    logsfr = np.log10(lumin.cgs.value) - 41.27

    header.set('log(SFR)',logsfr,comment='conversion from image counts to log(SFR (Msun/yr))')


## write out image with updated header
fits.writeto(args.image, imdat, header, overwrite=True)
