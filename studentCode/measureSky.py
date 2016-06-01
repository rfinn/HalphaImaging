#   This program takes a pair of R and Ha images and performs the sky 
# subtraction process on them.  
# Requires a file in the same directory named 'coordslist" providing 
# x,y coordinates of 30-40  positions to measure sky. Example included 
# in programs may be used; edit if significant numbers fall near stars. 
# (Check by displaying image and then using 
# cl> tvmark <frame number> coordslist )
# The coordslist file was made by using iraf fitsky and 
# pdumping xcenter,ycenter.
# For example, in iraf or pyraf:
# fitsky 9655
# <spacebar> around image for 30-40 positions without stars, then q, w to save.
# pdump 9655.sky.1 xcenter,ycenter > coordslist
#
# Example run command:
# python measureSky.py a9655 a9655_Ha
#
#Author: Ryan Muther, July 2013, Union College

import sys
import os
from pyraf import iraf
import numpy as np

def silentDelete(filename):
  try:
    os.remove(filename)
  except OSError:
    pass

def skySub(imageName):
  #Set up iraf packages and variables for filenames
  skyIm = imageName+'.sky'
  outfile = 'mSky'+imageName+'.file'
  avgValLog = 'mSkyLog'+imageName+'.file'
  stDevLog = 'stDevLog'+imageName+'.file'
  outputIm = imageName+'.sky'
  iraf.apphot(_doprint=0)
  
  #Delete exisiting files 
  for f in (skyIm,outfile,'mSky.file',outputIm+'.fits',avgValLog):
    silentDelete(f)
  
  #Configure fitsky parameters
  iraf.fitsky.coords="coordslist"
  iraf.fitsky.interactive='no'
  iraf.fitsky.verify='no'
  iraf.datapars.datamax=150000
  iraf.datapars.datamin=-200
  iraf.fitskypars.salgorithm="mode"
  iraf.fitskypars.annulus=0
  iraf.fitskypars.dannulus=40
  iraf.fitskypars.shireject=2
  iraf.fitskypars.sloreject=2
  iraf.fitsky.verbose='no'
  
  #Call fitsky to compute the mean of the sky
  iraf.fitsky(imageName,output=outfile)
  
  #Use 'PDUMP" to extract only the mean sky values
  sys.stdout=open(avgValLog,"w")
  iraf.pdump(outfile,'mSky',"yes")
  sys.stdout=sys.__stdout__

  #Use pdump to get the standard deviation from the sky
  sys.stdout=open(stDevLog,"w")
  iraf.pdump(outfile,'stDev',"yes")
  sys.stdout=sys.__stdout__

  #Read in the mSky values file
  skyValues = open(avgValLog,'r')
  valueList = skyValues.readlines()
  
  #convert the value array to a numpy array
  valueListNp = np.asfarray(valueList)

  #compute average
  average = np.mean(valueListNp)
  
  #Print output 
  print(str(average)+" is the average")
  
  #Read in the stDev log, take the average
  stDevs = open(stDevLog,'r')
  stDevList = stDevs.readlines()
  stDevListNp = np.asfarray(stDevList)
  stDevAvg = np.mean(stDevListNp)
  
  #Configure imarith
  iraf.imarith.pixtype="real"
  iraf.imarith.calctype="real"
  
  #Use imarith to subtract the average sky value from the image
  iraf.imarith(imageName,"-",average,outputIm)
  
  #Add the average sky value and stDev to the image header
  iraf.hedit(outputIm,"SKY",average,add='yes',verify='no')
  iraf.hedit(outputIm,"SKYSIGMA",stDevAvg,add='yes',verify='no')
  
  #Print an output image
  print("Subtracted "+ str(average)+ " to create a new image called "+outputIm)
  
  #Delete acillary output files
  silentDelete(outfile)

#Takes commandline arguments and runs the above function on them
filename = sys.argv[1]
filename_Ha = sys.argv[2]

print('R...')
skySub(filename)
print('Ha...')
skySub(filename_Ha)
