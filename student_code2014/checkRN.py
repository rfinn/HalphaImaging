# Checks the date in a fits image header, and if the year is 2011 
# sets the RDNOISE key to 26.
# The fact that this program needs to exist saddens me.
#
# Author: Ryan Muther, Union College, July 2013
#
# This will likely only work on linux.  So if it breaks on something else, 
#check the function getFitsFiles().  It concatenates the directory the 
# program is currently working in a linux-specific manner

from astropy.io import fits
from pyraf import iraf
import os

def checkRN(filename):
  imageHeader = fits.open(filename)[0]
  hasDateObsKey = True
  dateString = ""
  
  try:
    dateString = imageHeader.header['DATE-OBS']
  except KeyError:
    hasDateObsKey = False
      

  year = dateString[:4]
  
  if (hasDateObsKey and year=='2011'):
    iraf.hedit(filename,'RDNOISE',26,add='yes',verify='no')
    print(" ")
  
  

#Returns a list of all the fits files in a directory
#The a and arguments is completely irrelevent.
def getFitsFiles(a,dirname,files):
  for entry in files:
    if entry.endswith('.fits'):
      checkRN(dirname+"/"+entry)

os.path.walk(os.getcwd(),getFitsFiles,None)
print('Done!')
