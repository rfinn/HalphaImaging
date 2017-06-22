#!/usr/bin/env python

import glob
from astropy.io import fits
import numpy as np 
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.modeling.models import Gaussian2D

#from run_sextractor import *
#The goal of this program is to have a python-based convolution routine


def get_fwhm():
    for i in range(nfiles):
        hdulist = fits.open(files[i])
        data = hdulist[2].data
    
        fwhm = np.mean(data.FWHM_IMAGE)
        fwhm_std = np.std(data.FWHM_IMAGE)
    
        image_fwhm[i] = fwhm
        image_fwhm_std = fwhm_std
        
def convolve_images():    
    for i in range(nfiles):
        hdulist = fits.open(files[i])
        data = hdulist[2].data
        coords = [data.X_IMAGE, data.Y_IMAGE]
        convolved_image = 'g'+ files[i]
        outfile = convolve(coords, kernel)
        fits.writeto(convolved_image,np.array(outfile))
    print np.shape(outfile)
#get all files for a  given object
prefix=raw_input('Give image prefix (EX: ifwcs_data08???.fits)')
#prefix = 'A1367_R.cat' ###I forgot how raw works so this was a test
files = glob.glob(prefix)
nfiles = len(files)
image_fwhm = np.zeros(nfiles,'f')
image_fwhm_std = np.zeros(nfiles,'f')   
#convolved_image = np.zeros(nfiles,'f')
get_fwhm()

#Need worst FWHM to convolve to
fwhm_max=np.max(image_fwhm)
print 'the largest FWHM = ',fwhm_max

#(sigma_out)^2 = (sigma_in)^2 + (sigma_filter)^2
sigma_filter = np.sqrt((fwhm_max/2.35)**2 - (1./2.35)**2)
kernel = Gaussian2DKernel(sigma_filter)

convolve_images()