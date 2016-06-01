# This program creates stellar masks for R images, finding the stars in the
# image using daophot within iraf.  An  elliptical area to be excluded 
# should be specified as 'galcut.file' containing the x,y coordinates of the 
# center, the semimajor axis, the eccentricity and the position angle of 
# the ellipse in that order. Make the semimajor axis large enough to eliminate
# the area of the galaxy containing Halpha emission. (This is because daophot
# will identify bright HII regions in the R image as stars.)
#
# Requires an R image with SKYSIGMA, SEEING, and SKY
# Requires file named 'galcut.file' in same directory. 
# For example, output of 'more galcut.file' might show:
#  500 610 100 .4 34
# xcenter,ycenter, max semimajor axis, b/a, pa of ellipse
#
# Example run command:
# python maskstars.py 9655
#
# If too few or many stars are found in the image, adjust the
# findpars.threshold parameter.
#
# Output: [input_name]mask.fits :  mask of 0,1 values    (e.g., 9655mask.fits)
#         [input_name]Masked.fits : original image with mask values superposed  (e.g., 9655Masked.fits) 
#         [input_name].coo : coordinate file of all stars found
#          maskfile :    all stars outside aperture
#          maskfile.sat : all saturated stars outside aperture
#Author: Ryan Muther, Union College, July 2013 based on maskstars.cl        
#  originally written by Rebecca Koopmann
# Updated March 27, 2014 by R. Koopmann with changes:
#   - daophot.datamin corrected to -3*sigma
#   - cleaner end if no stars are found
#   - additional documentation

from __future__ import print_function
from astropy.io import fits
import numpy as np
from pyraf import iraf
import sys
import os
import math


#Attempts to delete a file, raising no error if the file does not exist
#Saves me from writing all the try/except blocks
def silentDelete(filename):
  try:
    os.remove(filename)
  except OSError:
    pass

def maskStars(imageName):
  #read in the fits data and header
  imageHeader = fits.open(imageName+".fits")[0]
  
  #Read in the skysigma, skyv, and seeing (fwhm) from the header and set up other relevant values
  sigma = imageHeader.header['SKYSIGMA']
  threshold = 3.0 * sigma
  
  skyv = imageHeader.header['SKY']
  
  fwhm = imageHeader.header['SEEING']
  aper = 3.0 * fwhm
  annul = 5 * fwhm
  
  #Read in the galcut file
  galcut = open("galcut.file")
  galcutString = file.read(galcut)
  galcut.close()
  
  #Parse the galcut file into its variables
  galList = galcutString.split(' ')
  xCenter = float(galList[0])
  yCenter = float(galList[1])
  a = float(galList[2]) #Semimajor Axis
  eccent = float(galList[3]) #Eccentricity
  posAngle = float(galList[4]) #Position angle
  
  posAngleRad = (posAngle+90)*3.14159/180 
  b = a * eccent #Semiminor Axis
  na = -1*a
  nb = -1*b
  
  #Create a version of the r image with the sky added back in for the daofind procedure
  iraf.imdelete("withsky")
  iraf.imarith(imageName,'+',skyv,"withsky")

  # Load daophot package
  iraf.digiphot()
  iraf.daophot()  

  print("A. Find stars using daofind")

#DAOFIND  searches  the  IRAF  images image for local density maxima,
#    with a  full-width  half-maxima  of  datapars.fwhmpsf,  and  a  peak
#    amplitude  greater  than  findpars.threshold  * datapars.sigma above
#    the local background, and writes a list of detected objects  in  the
#    file  output.  
# Here we set the fwhm and sigma equal to those listed in the header.
# We set the datamin to -3*sigma = -1*threshold defined above
# The findpars.threshold is set to 4.5*sigma by default. You may have
# to alter this value for images where too few/many stars are found.


  #Configure 'datapars', 'findpars', and 'daofind'
  iraf.datapars.fwhmpsf=fwhm
  iraf.datapars.sigma=sigma
  iraf.datapars.datamax='indef'
  iraf.datapars.datamin=(-1)*threshold
  iraf.datapars.ccdread='RDNOISE'
  iraf.datapars.gain="GAIN"
  iraf.datapars.exposure="EXPTIME"
  iraf.datapars.airmass="AIRMASS"
  iraf.datapars.filter="FILTER"
  iraf.datapars.obstime="TIME-OBS"

  iraf.findpars.threshold=4.5
  iraf.findpars.roundhi=3
  iraf.findpars.roundlo=-3

  iraf.daofind.verify='no'
  iraf.daofind.verbose='no'
  
  #create the variable for the coordinate file
  allStars = imageName+'.coo'
  
  #Remove a previous coordinate file if one exists
  silentDelete(imageName+'.coo')

  #Find the stars in the aligned R image
  iraf.daofind("withsky",allStars)
  print("The file containing the data of the stars in the aligned R-image is ",allStars)
  print(" ") 

  if os.path.getsize(allStars) > 2239 :

   print("B. Separate stars inside and outside galaxy")
  
   #Delete existing files
   for f in ("temp.file","maskfile","maskfile.sat"):
     silentDelete(f)
  
   #Extract the coordinates from the allStars file, redirecting stdout to do so
   sys.stdout=open("temp.file","w")
   iraf.pdump(allStars,"xcenter,ycenter",'yes')
   sys.stdout = sys.__stdout__
  


  #Read the file to get the xy coordinate pairs in a list
   coords = open("temp.file")
   coordList = list(coords)
   countout=0
  
  #Create strings of xy pairs to write to files
   outGalString = ""
   for pair in coordList:
     #for each entry in the list, parse it into xy value integer pairs
     values = pair.split("  ")
     x = float(values[0])
     y = float(values[1])
    
    #Determine if the star lies within the galaxy-centered ellipse
     inGalaxy = False
     deltaX = x - xCenter
     deltaY = y - yCenter
     cdx = abs(deltaX)
     cdy = abs(deltaY)
     if (cdx < a and cdy < a):
       xt=(deltaX*math.cos(posAngleRad))+(deltaY*math.sin(posAngleRad))
       yt=(deltaY*math.cos(posAngleRad))-(deltaX*math.sin(posAngleRad))
       if(xt <= a and xt >= na and yt <= b and yt >= nb):
         inGalaxy = True
      
     if (not inGalaxy):
       outGalString += str(x) + " " + str(y) + "\n"
       countout=countout+1
 
  #Write the coordinate list to a file
   mask = open("maskfile","w")
   mask.write(outGalString)
   mask.close()
   print('')
   print("  ",countout," stars found outside galaxy")
   print('')
   if countout!=0 :

    print("C. Do photometry to find saturated stars")
  
    for f in ("maskfile.mag","maskfile.satmag","maskfile.sat"):
     silentDelete(f)
  
  #Configure 'centerpars', 'fitskypars','photpars'
    iraf.datapars.datamax=150000
    iraf.datapars.datamin='indef'
    iraf.centerpars.calgorithm="none"

    iraf.fitskypars.salgorithm="mode"
    iraf.fitskypars.annulus=annul
    iraf.fitskypars.dannulus=10

    iraf.photpars.apertures=fwhm
    iraf.phot.verify='no'
    iraf.phot.verbose='no'
  
  #Call 'phot' to get the data of the stars 
    iraf.phot (imageName,"maskfile","maskfile.mag")

    boexpr="PIER==305"
    iraf.pselect("maskfile.mag","maskfile.satmag",boexpr)

    sys.stdout=open("maskfile.sat","w")
    iraf.pdump("maskfile.satmag","xcenter,ycenter","yes")
    sys.stdout = sys.__stdout__
  
    print("D. Make mask of stars outside galaxy")
  
  #Delete the rmask if one exists
    silentDelete("rmask.fits")
  
  #Configure 'imedit'
    iraf.imedit.cursor="maskfile"
    iraf.imedit.display='no'
    iraf.imedit.autodisplay='no'
    iraf.imedit.aperture="circular"
    iraf.imedit.value=-1000
    iraf.imedit.default="e"
  
  #Replace the pixel values near the star coordinates with value -1000 using imedit
    iraf.imedit(imageName,"redit",radius=aper)
 
    iraf.imcopy.verbose='no'
  #Do the same for the saturated stars, but with larger apertures
    satMaskSize = os.path.getsize("maskfile.sat")
    if satMaskSize < 1 : 
     print("   No saturated stars found")
    if (satMaskSize!=0):
     iraf.imedit.cursor="maskfile.sat"
     iraf.imedit.radius=2*aper
     iraf.imedit("redit","rmask")
    else:
     iraf.imcopy("redit","rmask")
  
  #Replace good pixels with value 0 in 'rmask'
    iraf.imreplace.lower=-999
    iraf.imreplace.upper='indef'
    iraf.imreplace("rmask",0)

  #Replace bad pixels with value 1 in 'rmask'
    iraf.imreplace.lower='indef'
    iraf.imreplace.upper=-1000
    iraf.imreplace("rmask",1)
  
   
  #Remove the mask file if one already exists
    silentDelete(imageName+"mask.fits")
    iraf.imcopy("rmask",imageName+"mask")

    print(" ")
    print("The mask image is ",imageName+"mask")
    print("The masked image is ",imageName+"Masked")
    print(" ")
 
   else :
    print("No stars found outside galaxy, no mask created.")

  if os.path.getsize(allStars) < 2314 :
   print("No stars found, no mask created.")



  #Delete intermediate images 
  for f in ("rmask.fits","redit.fits","maskimage.fits","maskfile.mag","maskfile.satmag","temp.file","withsky.fits","mask.fits"):
    silentDelete(f)
##################### Everything after this point reads the input from the
#command line and runs the above function on it.

filename = sys.argv[1]
maskStars(filename)
