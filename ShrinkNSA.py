#!/usr/bin/python
# coding: utf-8

# Creates a smaller version of the NSA catalog containing only the RA, DEC, and Z of objects with a DEC > -5


from astropy.io import fits
import numpy as np


nsadat =fits.getdata('nsa_v0_1_2.fits')


flag = [nsadat.DEC >-5]

d = nsadat.DEC[flag]
r = nsadat.RA[flag]
z = nsadat.Z[flag]

colr = fits.Column(name='RA',format=nsadat.formats[2],array=r)
cold = fits.Column(name='DEC',format=nsadat.formats[3],array=d)
colz = fits.Column(name='Z',format=nsadat.formats[11],array=z)
newcol = fits.ColDefs([colr,cold,colz])
hdu = fits.BinTableHDU.from_columns(newcol)
hdu.writeto('nsa_uat.fits',clobber=True)

