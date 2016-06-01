#!/usr/bin/env python

####################################################################
# PURPOSE
# This program interactively creates or edits an existing stellar 
# mask image using the iraf/pyraf imedit routine. It can be used
# to mask stars near a galaxy in a stellar mask created by maskstars.py
# or similar program or it can be used as a stand-alone routine. An
# Ha mask is made during the same process, as long as at least one area
# is masked. An Ha mask is expected for the photometry programs, so
# make sure an Ha mask is created by masking at least one region.
# 
# Example run command:
#  python handpick.py a9655 
#
# Note that ds9 must be running BEFORE you run handpick.
#
# INPUT
# Requires R and continuum-subtracted Ha images with SEEING provided
# as a header parameter. The R image name is provided (without .fits)
# in the command line as shown. An Halpha image named <galaxyName>_Ha, 
# e.g., a9655_Ha in above example is required to exist in the same
# directory. If editing an existing R mask, the mask image
# must have name <galaxyname>mask.fits, e.g., a9655mask.fits. If a 
# previous "handmask" exists of format # <galaxyname>handmask.fits exists, 
# it will be used in place of the "mask" image, so that this task may be 
# used iteratively.
#
# DESCRIPTION
# The process of using the program itself is fairly straightforward, with 
# commands displayed on the screen at pertinent points. The iraf task imedit
# is used to designate masked pixels. The default size of the stellar mask
# is 3*SEEING, where SEEING must be provided in the R image header. This
# is typically the minimum radius of a stellar mask, and should be increased
# where necessary, e.g. for saturated stars and background galaxies. Standard
# imedit tools are used to change the radius, e.g., +/- to increment by +/-1,
# :radius value to set to a particular value.
#
# The stdimage must be set in the login.cl file. Use "set stdimage=imt55" for 
# MOSAIC data, "set stdimage=imt4096" for HDI data.
#
# The program contains default z1 and z2 for use in displaying. You may change
# these as needed before running the program. If the display is not adequate,
# experiment with values using the display command, then edit the program to
# the optimal values.
#
# OUTPUT
#  Two versions of the mask for each image
#    R Image: [input_Image_name]+'handmask.fits'
#             [input_Image_name]+'handmasked.fits'
#   Ha Image: [input_Image_name]+'_Hahandmask.fits'
#             [input_Image_name]+'_Hahandmasked.fits'
#   The "handmask" files are the masks, with good pixels set to 0 and bad pixels
#   set to 1, which is the expectation for .pl files used by IRAF ellipse.
#   The "handmasked" files are created so that users may examine these separately
#   within ds9. They show the background image with masked regions superposed. 
#   These are useful for testing whether stellar outskirts are sufficiently masked

#Author: Ryan Muther, Union College, June-July 2013
#This program was written as a replacement for the IDL program of the same 
#name, originally written by Cullen Blake 1999.
# Updated May 8, 2014 by Rebecca Koopmann:
#   - added documentation

from pyraf import iraf
from astropy.io import fits
import numpy as np
import sys

iraf.set(stdimage='imt2048')
# Load the bad pixel mask and reverse its values so the mask exists in a sensible manner.
# Currently, the zeroes represent good pixels and the ones represent masked pixels.  
# This is the IRAF default, but less useful than the reverse when manipulating the masks
# and images. This reversal process allows for the creation of a masked image by simple 
# elementwise multiplication of the  mask and the image to be masked

def reverseMask(filename):
  mask = fits.getdata(filename+"mask.fits")
  
  manyZeros = np.zeros_like(mask)
  manyOnes = np.ones_like(mask)
  
  maskReversed = np.where((mask!=0),manyZeros,manyOnes)
  
  hdu=fits.PrimaryHDU(maskReversed)
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(filename+'MaskReversed.fits')

def maskImages(filename):
  #Read in the full width half maximum from the input file
  imageHeader = fits.open(filename+".fits")[0]
  fwhm = imageHeader.header['SEEING']

  ######CONFIG STUFF######
  z1_default=-12 #Change this value to alter the low end of the image display
  z2_default=1000 #Change this value to alter the high end of the image display
  ztrans_default='log'  #Change this to 'log' to use a logarithmic intensity plot in the display of the R image
  iraf.imedit.radius=3.0*fwhm  #Change this value to change imedit's default radius.
  ########################


  #Determine if a masked R image exists; create one if it does not.
  maskedExists=False
  try:
    with open(filename+'Masked.fits'):maskedExists=True
  except IOError:
    print 'The masked file does not exist.  Creating masked file...'
    mask = fits.getdata(filename+'MaskReversed.fits')
    base = fits.getdata(filename+'.fits')
    result = np.multiply(base,mask)
    
    hdu=fits.PrimaryHDU(result)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(filename+'Masked.fits')
    print 'Masked file written as testDataMasked.fits'
    
  
  #Display the image's Ha counterpart
  iraf.display(filename+'_Ha.fits',2)
  
  #Set the display parameters for the R image
  iraf.display.z1=z1_default
  iraf.display.z2=z2_default
  iraf.display.ztrans=ztrans_default
  iraf.imedit.command="display $image 1 erase=$erase fill=no order=0 >& dev$null"
  
  #Determine if a handmasked image exists
  handmaskedExists=False
  try:
    with open(filename+'_Handmasked.fits'):handmaskedExists=True
    print 'Handmasked image detected, loading it for further editing...'
  except IOError:
    print 'Creating new handmask image...'
  
  #Print some useful commands
  print 'Here are some pertinent commands used in imedit\ne  Fill circles   with constant value\n+/-  increase/decrease radius\nk  Fill circles with input data (spot healing, essentially)\nu  Undo the last change\nq  Finish editing the image and write the output\nQ Finish editing without writing output\ni  Reset editor to initial loaded image'
  
  #if a handmasked image does not exist, load the R image for editing otherwise load the handmask to edit it.
  iraf.display.zscale='no'
  iraf.display.zrange='no'
  if (not handmaskedExists):
    iraf.imedit(filename+'Masked.fits',filename+'_Handmasked.fits')
  else:
    iraf.imedit(filename+'_Handmasked.fits','')
  
  #Repeat the process for Ha.
  print 'Mask Ha.  Same commands apply'
  haHandmaskedExists=False
  iraf.display.zscale='yes'
  iraf.display.zrange='yes'
  try:
    with open(filename+'_HaHandmasked.fits'):haHandmaskedExists=True
    print 'Ha handmasked image detected, loading it for further editing...'
  except IOError:
    print 'Creating new Ha handmask image...'
    
  
  iraf.display(filename+'Masked.fits',2) 
  
  if (haHandmaskedExists):
    iraf.imedit(filename+'_HaHandmasked.fits','')
  else:
    iraf.imedit(filename+'_Ha.fits',filename+'_HaHandmasked.fits')

#This reads through a handmasked set of R and Ha images and extracts the #masks from it, writing it as (filename)handmask.fits and (filename
#)_HaHandmask.fits
def getMasks(filename):
  #Delete the old masks, if any exist.
  iraf.imdelete(filename+'_Handmask.fits')
  iraf.imdelete(filename+'_HaHandmask.fits')
  
  #Extract the mask from the masked R image, using a handmasked image if one exists otherwise using the masked image
  maskValue = 0.0
  try:
    maskedRImage = fits.getdata(filename+'_Handmasked.fits')
  except IOError:
    maskedRImage = fits.getdata(filename+'Masked.fits')
  manyZeros = np.zeros_like(maskedRImage)
  manyOnes = np.ones_like(maskedRImage)
  
  print 'Writing R mask'
  RmaskPixels = np.where((maskedRImage!=maskValue),manyZeros,manyOnes)
  
  hdu=fits.PrimaryHDU(RmaskPixels)
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(filename+'_Handmask.fits')
  
  #Do the same for the Ha image
  print 'Writing Ha mask'
  try:
    maskedHaImage = fits.getdata(filename+'_HaHandmasked.fits')
  except IOError:
    maskedHaImage = fits.getdata(filename+'_Ha.fits')
  print maskedHaImage.shape,manyZeros.shape
  HaMaskPixels = np.where((maskedHaImage!=maskValue),manyZeros,manyOnes)
  
  hdu=fits.PrimaryHDU(HaMaskPixels)
  hdulist = fits.HDUList([hdu])
  hdulist.writeto(filename+'_HaHandmask.fits')

######################
filename = sys.argv[1]
maskReversed=False
try:
  with open(filename+'MaskReversed.fits'):maskReversed=True
except IOError:
  print 'Reversing input mask'

if(not maskReversed):
  reverseMask(filename)

maskImages(filename)
getMasks(filename)
print('Removing superfluous files...')
iraf.imdelete(filename+'MaskReversed.fits')
print('Done!')
