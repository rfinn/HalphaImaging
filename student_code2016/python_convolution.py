#!/usr/bin/env python

import glob
from astropy.io import fits
import numpy as np 
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.modeling.models import Gaussian2D
import argparse
from argparse import RawDescriptionHelpFormatter

#from run_sextractor import *
#The goal of this program is to have a python-based convolution routine
parser = argparse.ArgumentParser(description ='This code will convolve image cutouts that have bad focus so that we can get a more precise continuum subtraction.')
parser.add_argument('--prefix', dest = 'prefix', default = False, action = 'store_true', help = 'Input the string of images to be convolved (before continuum subtraction, after cutouts). Enter prefix pointing (e.g. pointing-1)')
parser.add_argument('--filter', dest = 'filter' default = Ha, help = 'Input filter (e.g. Ha or R) for input string of images to convolve for each input pointing prefix')
args = parser.parse_args()

#search_prefix = args.string+'*-Ha.fits'
search_prefix = args.prefix+'*-'+args.filter+'.fits'
print search_prefix
input_images = glob.glob(search_prefix)



def get_fwhm():
    for image in input_images:
        hdulist = fits.open(input_images[image])
        data = hdulist[2].data
    
        fwhm = np.mean(data.FWHM_IMAGE)
        fwhm_std = np.std(data.FWHM_IMAGE)
    
        image_fwhm[image] = fwhm
        image_fwhm_std = fwhm_std
        
    
#All of the above is new  stuff

def convolve_images():    
    for image in input_images:
        hdulist = fits.open(input_images[image])
        data = hdulist[2].data
        coords = [data.X_IMAGE, data.Y_IMAGE]
        convolved_image = 'g'+ input_images[image]
        outfile = convolve(coords, kernel)
        fits.writeto(convolved_image,np.array(outfile))
    print np.shape(outfile)

files = glob.glob(args.prefix)
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
