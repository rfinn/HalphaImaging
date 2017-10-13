#This program takes an R image and a Ha image that have been aligned
# (for example by contSub.py) and an input continuum subtraction factor
# to create a continuum-subtracted Ha image
#
#Requires an R image, an Ha image, and the continuum subtraction
#factor to be passed as arguments to the program
#Images should be background-subtracted.
#
# This program may be used to fine-tune automatic continuum-subtraction
# as output by contSub.py
#
# Example run command:
# python contSubOnly.py a9655R a9655_Ha 0.54
#
# OUTPUT
#   Image [input_Image_name]+'cs.fits'
#
# Author: Rebecca Koopmann, Union College, February 2014
# based on contSub.py by Ryan Muther, Union College


import sys
import os
from pyraf import iraf
from astropy.io import fits
import numpy as np

#Attempts to delete a file, raising no error if the file does not exist
#Saves me from writing all the try/except blocks
def silentDelete(filename):
  try:
    os.remove(filename)
  except OSError:
    pass

#Scales and subtracts the R image from the Ha image, leaving the continuum
#subtracted Ha image behind.

def continuumReduce(imageName,imageHaName):

  avgScaleFactor = contsubfactor

  print("The scale factor is "+str(avgScaleFactor))
  
  #configure hedit
  iraf.hedit.add='yes'
  iraf.hedit.verify='no'
  
  iraf.imarith(imageName+'.fits',"*",avgScaleFactor,imageName+"_scaled.fits")
  iraf.imarith(imageHaName+'.fits',"-",imageName+"_scaled.fits",imageHaName+"cs.fits")
  iraf.hedit(imageHaName+'cs.fits','Rscale',avgScaleFactor) 
  
  #clear up superfluous images
  for f in (imageName+'_scaled.fits'):
    silentDelete(f)

  print('The continuum subtracted image is '+ imageHaName+'cs.fits')

#This bit just takes the arguments from the command line 
filename_R = sys.argv[1]
filename_Ha = sys.argv[2]
contsubfactor= sys.argv[3]

print('Subtracting continuum...')
continuumReduce(filename_R,filename_Ha)
