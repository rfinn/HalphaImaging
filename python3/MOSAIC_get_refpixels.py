#!/usr/bin/env python
"""
Generate the list of x and y ref pixels for theli config file
"""
from astropy.io import fits


sx = ''

sy = ''

infile = 'data08607.fits'

t = fits.open(infile)

for i in range(1,len(t)):
    sx = sx+f"[{i}]={t[i].header['CRPIX1']:.0f} "
    sy = sy+f"[{i}]={t[i].header['CRPIX2']:.0f} "

print(sx)

print(sy)
