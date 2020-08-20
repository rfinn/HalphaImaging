#!/usr/bin/env python

'''
GOAL:
  The goal of this program is to convolve the R and Halpha mosaics so that they have the same FWHM.
  This should be done BEFORE subtracting the continuum.

PROCEDURE:
  - before running, run sextractor.  For example:

       uat_astr_mosaic.py --s --filestring 'pointing-16*coadd.fits'
  - read in the sextractor catalogs and calculate the FWHM of each image
  - 


EXAMPLE:



INPUT/OUPUT:

REQUIRED MODULES:

EXTRA NOTES:

WRITTEN BY:
Rose Finn

EDITED BY:

'''


#!/usr/bin/env python

import glob
from astropy.io import fits
import numpy as np 
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.modeling.models import Gaussian2D
import argparse
from argparse import RawDescriptionHelpFormatter
from matplotlib import pyplot as plt


def get_fwhm(input_images): #measure FWHM of SE catalogs
    nfiles = len(input_images)
    image_fwhm = np.zeros(nfiles,'f')
    image_fwhm_std = np.zeros(nfiles,'f')   
    for i in range(len(input_images)):
        t = input_images[i].split('.fits')
        se_cat = t[0]+'.cat'
        data = fits.getdata(se_cat,2)
        # select unsaturated stars using : class_star > 0.9 and (10 < m < 13)
        flag = (data.CLASS_STAR > 0.9) & (data.MAG_AUTO > 10.) & (data.MAG_AUTO < 13.)
        image_fwhm[i] = np.mean(data.FWHM_IMAGE[flag])
        image_fwhm_std[i] = np.std(data.FWHM_IMAGE[flag])
    return image_fwhm, image_fwhm_std
    
#All of the above is new  stuff


#from run_sextractor import *
#The goal of this program is to have a python-based convolution routine
parser = argparse.ArgumentParser(description ='This code will convolve image cutouts that have bad focus so that we can get a more precise continuum subtraction.')
parser.add_argument('--prefix', dest = 'prefix', default = 'pointing-1',  help = 'Input the string of images to be convolved (before continuum subtraction, after cutouts). Enter prefix pointing (e.g. pointing-1)')
parser.add_argument('--test',dest = 'test', default = False, action='store_true')
args = parser.parse_args()



#search_prefix = args.string+'*-Ha.fits'
search_prefix = args.prefix+'*coadd.fits'
print search_prefix
input_images = glob.glob(search_prefix)

#convolved_image = np.zeros(nfiles,'f')
image_fwhm, image_fwhm_std = get_fwhm(input_images)

#Need worst FWHM to convolve to
fwhm_max=np.max(image_fwhm)
print 'the largest FWHM = ',fwhm_max

if args.test:
    #try to double the FWHM of one image
    sigma_filter = np.sqrt((2*image_fwhm[0]/2.35)**2 - (image_fwhm[0]/2.35)**2) # filter to double FWHM
    print sigma_filter
    kernal = Gaussian2DKernel(sigma_filter)
    input_image = fits.getdata(input_images[0])
    plt.imshow(input_image)
    outfile = convolve(input_image,kernel)
    plt.imshow(outfile)
    fits.writeto('test.fits',np.array(outfile))
else:


    # convolve all images to worst seeing
    # use pyraf.iraf.gauss (sigma = FWHM/2.35)
    # (sigma_out)^2 = (sigma_in)^2 + (sigma_filter)^2
    #
    # (sigma_filter) = np.sqrt[(fwhm_out/2.35)^2 - (sigma_in/2.35)^2] 
    sigma_filter = np.sqrt((fwhm_max/2.35)**2 - (image_fwhm/2.35)**2)
    #convolve_images()
    for i in range(len(input_images)):
        #if image_fwhm[i] == fwhm_max:
        #    continue
        imdata = fits.getdata(input_images[i])
        convolved_image = 'g'+ input_images[i]
        kernel = Gaussian2DKernel(sigma_filter[i])
        outfile = convolve(imdata, kernel)
        fits.writeto(convolved_image,np.array(outfile))
        print np.shape(outfile)




#convolve_images(kernal=kernal)

'''
#BELOW IS KELLYS STUFF
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
'''




'''
original sketch by RF


# run sextractor

# read in output



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




'''
