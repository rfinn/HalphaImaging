#!/usr/bin/env python

'''
BASIC INFORMATION ABOUT THIS CODE:
  -This is the FIFTH program you need to run in order to perform the reduction o  f HDI images.
  -This code will fix all the headers to contain basic WCS information so that i  n the future we can feed the fits images into stacking programs like SCamp and  SWarp. 
  -Updates CMMTOBS --> FILTER
        RASTRNG --> CRVAL1
        DECSTRNG --> CRVAL2
  -Adds CRPIX1, CRPIX2, CD1_1, CD2_2, CTYPE1, CTYPE2

BEFORE RUNNING THIS CODE:
  -Be sure to type ur_setup in the terminal everytime you open a new terminal     window so that ureka is activated.
  -Ensure that pyraf is still activated by retyping the commands listed in the c  omments of the FIRST program titled "uat_HDIgroupflatfiles.py".

PROCEDURE:
  -This code uses the python task .rename and .append to the list of header information so that the header contains basic WCS information that can be easily read by astronomy programs. 
  -In addition, we used different arguments using args.parse to call on adding the different header information.
     
GOAL:
    The goal of this code is to successfully add header information that is not included in the HDI raw data frames.

EXAMPLE:
   In the directory containing all flattened objects with incorrect headers type in the command line:
      /home/share/research/pythonCode/uat_HDIfixheader.py

      (or whatever the path is to where this program is stored)
   
INPUT/OUPUT:
    Input --> all ftr*.fits in directory
    Output --> hftr*.fits

REQUIRED MODULES:
    argparse
    glob
    astropy
    
EXTRA NOTES:
can be used on dome flattened images. To do so type '--filestring "dtr*.fits"' after the command

WRITTEN BY:
Rose Finn

EDITED BY:
Research Team Summer 2016 --> Grant Boughton, Natasha Collova, Sandy Spicer

'''

import argparse
import glob
#from astropy import coordinates as coord
from astropy import units as u
import ccdproc
#from astropy.coordinates import SkyCoord
from astropy.io import fits

parser = argparse.ArgumentParser(description ='Remove cosmic rays using LAcosmic')
parser.add_argument('--filestring', dest='filestring', default='ftr*.fits', help='match string for input files (default =  ftr*.fits)')
#parser.add_argument('--', dest='pixelscalex', default='0.00011808', help='pixel scale in x (default = 0.00011808)')
#parser.add_argument('--pixscaley', dest='pixelscaley', default='0.00011808', help='pixel scale in y (default = 0.00011808)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring))
nfiles=len(files)
i=1
for f in files:
    print 'ZAPPING COSMIC RAYS FOR FILE %i OF %i'%(i,nfiles)
    data = ccdproc.CCDData.read(f,unit = u.adu)
    cr_cleaned = ccdproc.cosmicray_lacosmic(data, sigclip=5)
    #print 'WRITING UPDATED FILE'
    cr_cleaned.write('c'+f,clobber=True)
    i += 1
    print '\n'

    
    


                   
