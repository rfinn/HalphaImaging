#!/usr/bin/env python

import glob
import numpy as np
#from uat_astr_mosaic --s import *
from pyraf.iraf import gauss

# run sextractor

# read in output

# select unsaturated stars using : class_star > 0.9 and (10 < m < 13)

# measure mean and std of FWHM

# repeat for all images

def get_fwhm():
    for i in range(nfiles):
        im=image(files[i])
        image_fwhm[i] = im.fwhm
        image_fwhm_std[i] = im.fwhm_std

def convolve_images():    
    for i in range(nfiles):
        convolved_image = 'g'+files[i]
        gauss(input=files[i],output=convolved_image,sigma=sigma_filter[i])

# get all files for a  given object
prefix=raw_input('Give image prefix (EX: ifwcs_data08???.fits)')
files=glob.glob(prefix)
nfiles=len(files)
image_fwhm=np.zeros(nfiles,'f')
image_fwhm_std=np.zeros(nfiles,'f')


get_fwhm()


for i in range(nfiles): print i, files[i],image_fwhm[i],image_fwhm_std[i]

    
# get worst fwhm
fwhm_max=max(image_fwhm)
print 'the largest FWHM = ',fwhm_max



# convolve all images to worst seeing
# use pyraf.iraf.gauss (sigma = FWHM/2.35)
# (sigma_out)^2 = (sigma_in)^2 + (sigma_filter)^2
#
# (sigma_filter) = sqrt[(fwhm_out/2.35)^2 - (sigma_in/2.35)^2] 
sigma_filter = sqrt((fwhm_max/2.35)**2 - (image_fwhm/2.35)**2)
#convolve_images()
