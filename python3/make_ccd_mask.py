#!/usr/bin/env python

'''
https://photutils.readthedocs.io/en/stable/background.html



'''

from astropy.io import fits
import ccdproc

def make_ccd_mask(flat1,flat2=None):
    ''' make bad pixel mask from flat image.  use ratio of flats if flat2 is given '''
    if flat2 is None:
        flatimage = fits.getdata(flat1)
    else:
        flat1 = fits.getdata(flat1)
        flat2 = fits.getdata(flat2)
        flatimage = flat1/flat2

    # use ccdmask to identify bad pixels

    t = ccdmask(flatimage)
