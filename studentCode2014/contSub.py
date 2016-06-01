# This program takes sky-subtracted R and Ha images (Ha that has not 
# yet been continuum subtracted) and performs alignment and contiuum 
# subtraction to create a continuum-subtracted Ha image
#
# Requires an R image and an Ha image to be passed as arguments to the program
# Requires SEEING, SKYSIGMA header fields to exist in the image headers
# Images should be background-subtracted.
#
# The program uses daofind to find stars in the R image, selects the 25
# brightest, non-saturated stars and feeds the coordinates of the stars 
# into the imalign routine to align and trim the input images. Photometry
# is performed on the stars in each image and an average scale factor 
# computed as the average of the ratio of the Ha photometry to the R photometry.
# The standard deviation in the scale factor is computed and used to find
# an upper and lower scale factor. BE AWARE THAT THESE WILL SIMPLY PUT THE
# CONTINUUM FACTOR IN THE RIGHT BALLPARK. IT IS UP TO THE USER/ANALYZER TO
# FINE-TUNE THE CONTINUUM SUBTRACTION FOR EACH GALAXY BEFORE PHOTOMETRY.  
#
# Any existing continuum-subtracted images and files of default name will be
# deleted before running. If you need to save images from previous runs,
# rename from the default.
# 
# Example run command:
#    python contSub.py a9655.sky a9655_Ha.sky
# where a9655.sky and a9655_Ha.sky were output from the measureSky.py routine.
#
# FEATURE:
# Sometimes if images are different sizes, trimming in the imalign routine
# fails. If this happens, trim images to the same size. To do this, use
#       imhead yourimages
# to find the x and y sizes, then trim using imcopy on one or both images:
#       imcopy yourimage[1:smallest_xsize,1:smallest_y_size]
# where smallest_xsize and smallest_y_size are the minimum values output
# by imhead.
#
# OUTPUT
#   Images
#     R_final.fits: aligned R image
#     Ha_final.fits: aligned Ha image (not continuum subtracted)
#     Hacs_final.fits: continuum-subtracted Ha using average scale factor
#     Hacs_finalplus.fits: continuum-subtracted Ha using average scale factor 
#                           + one standard deviation
#     Hacs_finalminus.fits: continuum-subtracted Ha using average scale factor 
#                           - one standard deviation
#   Files
#     coords: coordinates of stars in R image from daofind
#     psel : iraf output of R photometry used to select 25 brightest stars
#     pdump.file and pdump_Ha.file : coordinates and fluxes of the stars 
#                used to find the scale factor 

#Author: Ryan Muther, Union College, July 2013
# Updated May 8, 2014 by R. Koopmann with changes:
#   - corrected value of scale factor written to header to avgScaleFactor
#   - additional documentation


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

#Aligns the R and Ha images to the R image using photometry data from the R 
#image
def alignImages(imageName,imageName_Ha):
  #Read the fwhm and seeing from the image file
  imageHeader = fits.open(imageName+".fits")[0]
  fwhm = imageHeader.header['SEEING']
  annul = 5.0 * fwhm
  aper=3.0 * fwhm
  sigma = imageHeader.header['SKYSIGMA']
  iraf.daophot(_doprint=0)
  #Do 'daofind' on the image to locate the stars
  print("1. Find stars using 'daofind'")
  
  #Configure 'datapars', 'findpars', and 'daofind'
  
  iraf.datapars.fwhmpsf=fwhm
  iraf.datapars.sigma=sigma
  iraf.datapars.datamin=-10
  iraf.datapars.ccdread="RDNOISE"
  iraf.datapars.gain="GAIN"
  iraf.datapars.exposure="EXPTIME"
  iraf.datapars.airmass="AIRMASS"
  iraf.datapars.filter="FILTER"
  iraf.datapars.obstime="TIME-OBS"

  iraf.findpars.threshold=20*sigma
  iraf.findpars.sharplo=0.2
  iraf.findpars.sharphi=1.0
  iraf.findpars.roundlo=-1.0
  iraf.findpars.roundhi=1.0
  iraf.daofind.verify='no'
  iraf.daofind.verbose='no'

  #Delete exisiting coordinate,output files of 'phot',file containg data of 'good' stars, old images
  for f in("coord","mag","psel","R_final.fits","Ha_final.fits"):
    silentDelete(f)

  iraf.daofind(imageName,'coord')

  print("File containing the coordinates of the stars is coord")


  print(" ")
  #Configure 'centerpars', 'fitskypars','photpars'
  iraf.datapars.datamax=150000
  iraf.centerpars.calgorithm="centroid"
  iraf.centerpars.cbox=16
  iraf.centerpars.maxshift=3

  iraf.fitskypars.salgorithm="mode"
  iraf.fitskypars.annulus=annul
  iraf.fitskypars.dannulus=10

  iraf.photpars.apertures=aper

  iraf.phot.verify='no'
  iraf.phot.verbose='no'

  print("2. Obtain data of stars using 'phot'")
  #Call 'phot' to get the data of the stars 
  iraf.phot (imageName,'coord','mag')

  #sort in order of increasing magnitude of stars
  iraf.psort('mag',"mag")
  boundsig=sigma+2
  boexpr="CIER==0 && PIER==0 && STDEV <="+str(boundsig)  
  print("File containing the data of the stars in order of decreasing brightness is mag")

  print (" ")

  print("3. Select stars with low error, no bad pixels")
  #Select stars that have no centering error, skyerror <sig+2 and no bad pixels
  iraf.pselect ("mag" ,"psel" ,boexpr)

  print("File containing stars with low sky error,low centering error is psel")
  print(" ")

  #Renumber the ID number of the stars in order of increasing magnitude
  iraf.prenumber ("psel")


  #Delete existing files
  for f in ("pdump.file","stars25",'alR.fits','rIn.fits','haIn.fits'):
    silentDelete(f)

  print("4. Select the 25 brightest stars")
  iraf.pselect ("psel","stars25", "ID <=25")

  print("File containing the brightest 25 'good' stars is stars25")
  print(" ")

  #Pick out only the required data of stars from the .25 file
  sys.stdout=open("pdump.file","w")
  iraf.pdump ("stars25","xcenter,ycenter,flux","yes")
  sys.stdout = sys.__stdout__
  
  print("The coordinates and flux are stored in pdump.file")
  
  #Align images
  iraf.imcopy(imageName,'rIn')
  iraf.imcopy(imageName_Ha,'haIn')
  print("Aligning images")
  #iraf.imalign.verbose='no'
  iraf.imalign("rIn,haIn", "rIn","pdump.file","R_final,Ha_final")

#Scales and subtracts the R image from the Ha image, leaving the continuum
#subtracted Ha image behind.
def continuumReduce(imageName):
  
  #Read the fwhm from the Ha image file
  imageHeader = fits.open(imageName+".fits")[0]
  fwhm = imageHeader.header['SEEING']
  annul = 5.0*fwhm
  aper=3.0*fwhm
  
  #Configure 'datapars', 'findpars', 'phot' and 'daofind'
  iraf.datapars.datamax=150000
  iraf.centerpars.calgorithm="centroid"
  iraf.centerpars.cbox=16
  iraf.centerpars.maxshift=3

  iraf.fitskypars.salgorithm="mode"
  iraf.fitskypars.annulus=annul
  iraf.fitskypars.dannulus=10

  iraf.photpars.apertures=aper

  iraf.phot.verify='no'
  iraf.phot.verbose='no'
  iraf.phot.verbose='no'
  
 #Delete old photometry and old Ha images
  for f in ('mag_Ha','psel_Ha',"Hacs_final.fits","Hacs_finalplus.fits","Hacs_finalminus.fits",imageName+"_scaledminus.fits",imageName+"_scaledplus.fits",imageName+"_scaled.fits"):
    silentDelete(f)
  
  #Call 'phot' to get the data of the stars from the aligned Ha image
  iraf.phot("Ha_final.fits",'pdump.file','mag_Ha')

  #Pick out only the required data of stars from the .25 file
  sys.stdout=open("pdump_Ha.file","w")
  iraf.pdump ("mag_Ha","xcenter,ycenter,flux","yes")
  sys.stdout = sys.__stdout__
  
  print("The Ha coordinates and flux are stored in pdump_Ha.file")
  
  #Read the data from the pdump files
  r_Data = open('pdump.file')
  r_DataList = r_Data.readlines()
  r_Data.close()
  Ha_Data = open('pdump_Ha.file')
  Ha_DataList = Ha_Data.readlines()
  Ha_Data.close()
  
  #Compute the average scale factor between the R and Ha images
  scaleFactors=[]
  print('Scale factors: ')
  for index in range(0,len(r_DataList)):
    rLine = r_DataList[index].split('  ')
    rFlux = float(rLine[2])
    haLine = Ha_DataList[index].split('  ')
    haFlux = float(haLine[2])
    scaleFactor = haFlux/rFlux
    print(str(scaleFactor))
    scaleFactors.append(scaleFactor)
  
  print(" ")
  scaleFactors=np.array(scaleFactors)
  
  #Compute the avg scale factor and standard deviation and ask if it is acceptable, changing the value if it is not.
  avgScaleFactor = np.mean(scaleFactors)
  stDev = np.std(scaleFactors)
  print("The scale factor is "+str(avgScaleFactor))
  decision = raw_input('Is this acceptable (y/n)? ')
  avgChanged=False
  if (not((decision=='yes') or (decision=='y'))):
    avgScaleFactor=float(raw_input('Enter a new scalefactor: '))
    avgChanged=True
    print("Average scale factor is now: "+str(avgScaleFactor))
  
  #configure hedit
  iraf.hedit.add='yes'
  iraf.hedit.verify='no'
  
  if (avgChanged):
    #No adjustment
    iraf.imarith('R_final',"*",avgScaleFactor,imageName+"_scaled")
    iraf.imarith('Ha_final',"-",imageName+"_scaled","Hacs_final")
    iraf.hedit("Hacs_final",'Rscale',avgScaleFactor)
    #plus .01
    iraf.imarith('R_final',"*",avgScaleFactor+.01,imageName+"_scaledplus")
    iraf.imarith('Ha_final',"-",imageName+"_scaledplus","Hacs_finalplus")
    iraf.hedit("Hacs_finalplus",'Rscale',avgScaleFactor+.01)
    #minus .01
    iraf.imarith('R_final',"*",avgScaleFactor-.01,imageName+"_scaledminus")
    iraf.imarith('Ha_final',"-",imageName+"_scaledminus","Hacs_finalminus")
    iraf.hedit("Hacs_finalminus",'Rscale',avgScaleFactor-.01)
  else:
    #No stdev adjustment
    iraf.imarith('R_final',"*",avgScaleFactor,imageName+"_scaled")
    iraf.imarith('Ha_final',"-",imageName+"_scaled","Hacs_final")
    iraf.hedit("Hacs_final",'Rscale',avgScaleFactor)
    #plus one stdev
    iraf.imarith('R_final',"*",avgScaleFactor+stDev,imageName+"_scaledplus")
    iraf.imarith('Ha_final',"-",imageName+"_scaledplus","Hacs_finalplus")
    iraf.hedit("Hacs_finalplus",'Rscale',avgScaleFactor+stDev)
    #minus one stdev
    iraf.imarith('R_final',"*",avgScaleFactor-stDev,imageName+"_scaledminus")
    iraf.imarith('Ha_final',"-",imageName+"_scaledminus","Hacs_finalminus")
    iraf.hedit("Hacs_finalminus",'Rscale',avgScaleFactor-stDev)
  
  #clear up superfluous images
  for f in (imageName+"_scaledminus.fits",imageName+"_scaledplus.fits",imageName+"_scaled.fits",'rIn.fits','haIn.fits'):
    silentDelete(f)

  print('The continuum subtracted images are Hacs_final.fits, Hacs_finalminus.fits, and Hacs_finalplus.fits')

#This bit just takes the arguments from the command line amd runs the above 
#two functions on them.
filename_R = sys.argv[1]
filename_Ha = sys.argv[2]

iraf.digiphot()
print('Aligning images...')
alignImages(filename_R,filename_Ha)
print('Subtracting continuum...')
continuumReduce(filename_R)
