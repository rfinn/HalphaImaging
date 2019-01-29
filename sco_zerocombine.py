#!/usr/bin/env python

'''
 
  GOAL:
  Make lists of files that contain
  - darks with same filename suffix
  Then combine the darks

  PROCEDURE:
  - Use gethead to pull relevant header data
  - Combine darks according to exposure time 

  EXAMPLE:
     In the directory containing all darks type in the command line:
     
    %run ~/github/HalphaImaging/sco_darkcombine.py --filestring 'dark-'
  
  INPUT/OUTPUT:
  Input: 'dark' 
  Output: For combined darks --> 'darkD(exposure-time).fits' 


  REQUIRED MODULES:
  pyraf

  NOTES:
  in junkfile ftr flats still show. We changed the gethead requirements to only   bring in files that start with tr but the ftr files will not go away! =(
  
  WRITTEN BY:
  Rose A. Finn
  2019-Jan-26  
  
'''
import glob
import os
import numpy as np

import ccdproc
import astropy.units as u
import argparse

'''
from pyraf import iraf
iraf.noao()
iraf.imred()
iraf.ccdred()    
'''

parser = argparse.ArgumentParser(description ='Groups images by filter and creates flatfield images')
parser.add_argument('--filestring', dest='filestring', default='bias', help='match string for input files (default =  bias, which gets bias*.fits)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))

zeros = ccdproc.combine(fimages,method='average',sigma_clip=True,unit=u.adu)
print('writing fits')
zeros.write(args.filestring+'-combined.fits',overwrite=True)
    



                

