#!/usr/bin/env python

'''
GOAL:
   This code will remove cosmic rays.

EXAMPLE:
   In the directory containing all flattened objects with incorrect headers type in the command line:
      /home/share/research/pythonCode/uat_HDIfixheader.py

INPUT/OUPUT:
    Input --> all tr*.fits  -- best to run on the trimmed files
    Output --> ztr*.fits

REQUIRED MODULES:
    
EXTRA NOTES:


WRITTEN BY:
Rose Finn

EDITED BY:


'''

import argparse
import glob
#from astropy import coordinates as coord
from astropy import units as u
import ccdproc
#from astropy.coordinates import SkyCoord
from astropy.io import fits

parser = argparse.ArgumentParser(description ='Remove cosmic rays using LAcosmic')
parser.add_argument('--filestring', dest='filestring', default='tr*.fits', help='match string for input files (default =  tr*.fits)')
#parser.add_argument('--', dest='pixelscalex', default='0.00011808', help='pixel scale in x (default = 0.00011808)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring))
nfiles=len(files)
i=1
for f in files:
    print 'ZAPPING COSMIC RAYS FOR FILE %i OF %i'%(i,nfiles)
    data = ccdproc.CCDData.read(f,unit = u.adu)
    cr_cleaned = ccdproc.cosmicray_lacosmic(data, sigclip=5)
    cr_cleaned.write('z'+f,clobber=True)
    i += 1
    print '\n'

    
    


                   
