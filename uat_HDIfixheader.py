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
  -This code uses the python task .rename and .append to the list of header info  rmation so that the header contains basic WCS information that can be easily r  ead by astronomy programs. 
  -In addition, we used different arguments using args.parse to call on adding t   he different header information.
     
GOAL:
  -The goal of this code is to successfully add the correct header information t  hat previously was not there before.

EXAMPLE:
   In the directory containing all flattened objects with incorrect headers type in the command line:
      '/home/share/research/pythonCode/uat_HDIfixheader.py'(or whatever the path is to where this program is stored)
   
INPUT/OUPUT:
Input --> all ftr*.fits in directory
Output --> hftr*.fits

REQUIRED MODULES:
-pyraf

EXTRA NOTES:
can be used on dome flattened images. To do so type '--filestring "dtr*.fits"' after the command

WRITTEN BY:
Dr. Rose Finn
EDITED BY:
Research Team Summer 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood, Kaitlyn Hoag, Kelly Whalen

'''

import argparse
import glob
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits

parser = argparse.ArgumentParser(description ='Edit image headers to include basic WCS information to the HDI image headers')
parser.add_argument('--filestring', dest='filestring', default='ftr*.fits', help='match string for input files (default =  ftr*.fits)')
parser.add_argument('--pixscalex', dest='pixelscalex', default='0.00011808', help='pixel scale in x (default = 0.00011808)')
parser.add_argument('--pixscaley', dest='pixelscaley', default='0.00011808', help='pixel scale in y (default = 0.00011808)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring))
nfiles=len(files)
i=1
for f in files:
    print 'FIXING HEADER FOR FILE %i OF %i'%(i,nfiles)
    data, header = fits.getdata(f,header=True)
    header.rename_key('FILTER1','FWHEEL1')
    header.rename_key('FILTER2','FWHEEL2')

    FILTER = header['CMMTOBS']
    header.append(card=('FILTER',FILTER,'FILTER'))


    RASTRNG = header['RASTRNG']
    RA = coord.Angle(RASTRNG,unit=u.hour)
    header.append(card=('CRVAL1',RA.degree,'RA of reference point'))

    DECSTRNG = header['DECSTRNG']
    DEC = coord.Angle(DECSTRNG,unit=u.degree)

    EQUINOX = header['EQUINOX']
    header.append(card=('CRVAL2',DEC.degree,'DEC of reference point'))

    header.append(card=('CRPIX1','2048.','X reference pixel'))
    header.append(card=('CRPIX2','2048.','Y reference pixel'))
    header.append(card=('CD1_1',args.pixelscalex,'Pixel scale in X'))
    header.append(card=('CD2_2',args.pixelscaley,'Pixel scale in Y'))
    header.append(card=('CTYPE1','RA---TAN-SIP',''))
    header.append(card=('CTYPE2','DEC--TAN-SIP',''))
    header.append(card=('GAIN','1.3','gain (e-/ADU)'))
    header.append(card=('TELRA',RA.degree,'RA of reference point'))
    header.append(card=('TELDEC',DEC.degree,'DEC of reference point'))
    header.append(card=('TELEQUIN',EQUINOX,'Epoch (years)'))    
    print 'WRITING UPDATED FILE'
    fits.writeto('h'+f,data,header,clobber=True)
    i += 1
    print '\n'

    
    


                   
