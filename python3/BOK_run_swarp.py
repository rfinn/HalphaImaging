#!/usr/bin/env python

'''

ORGANIZING DATA
* working in /mnt/qnap_home/rfinn/Halpha/Bok

* moving short exposure time images to subdirectory junk



'''

import ccdproc as ccdp
import os
from astropy.io import fits



def combine_masks(weight_image,dq_image):
    '''
    combine the weight and data quality image
    combined weight = weight_image + 1000*data_quality_image
    then use the following flag in swarp
    WEIGHT_THRESH 1000

    INPUT: 
    * weight_image : this is the weight image
    * dq_image : data quality image, with nonzero indicating bad pixels

    OUTPUT:
    * combined_weight.fits : this is the combined weight image
    
    '''
    weight_hdu = fits.open(weight_image)
    dq_hdu = fits.open(dq_image)
    
    # loop over the 4 images, extensions 1-4
    for i in range(1,5):
        print('combining image number ',i)
        weight_hdu[i].data = weight_hdu[i].data + 1000*dq_hdu[i].data
    weight_hdu.writeto('combined_weight.fits',overwrite=True)

def run_swarp(image_list, weight_list):
    pass


def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i+1
    except UnboundLocalError:
        return 0


homedir = os.getenv("HOME")
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
runscamp=False
runswarp=True


rawdir = os.getcwd()
# get list of images (images only)


# sort images by location and filter
# alternatively could use object name,
# but not all are correct, so need to fix names


# create file list with r-band images
# create a file list with r-band weight images


# create file list with halpha images
# create a file list with halpha weight images


# run swarp to mosaic r-band

# run swarp to mosaic halpha

# alight r and halpha imaging

#keys = ['OBJECT', 'RA', 'DEC', 'FILTER', 'EXPTIME']
#ic = ccdp.ImageFileCollection(rawdir, keywords=keys, glob_include='*.fit*', glob_exclude='*test*.fit')

