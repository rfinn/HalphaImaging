#!/usr/bin/env python

'''
GOAL:
   This code will remove cosmic rays.

EXAMPLE:
   In the directory containing all flattened objects with incorrect headers type in the command line:

   from within ipython

   %run ~/github/HalphaImaging/uat_zipcosmic.py --filestring tr
   
INPUT/OUPUT:
    Input --> all tr*.fits  -- best to run on the trimmed files
    Output --> ztr*.fits

REQUIRED MODULES:
    

WRITTEN BY:
Rose Finn

EDITED BY:


'''

import argparse
import glob
#from astropy import coordinates as coord
from astropy import units as u
import ccdproc
from astropy.nddata import CCDData
#from astropy.coordinates import SkyCoord
from astropy.io import fits

parser = argparse.ArgumentParser(description ='Remove cosmic rays using LAcosmic')
parser.add_argument('--filestring', dest='filestring', default='tr*o00', help='match string for input files (default =  tr*o00, which gets tr*o00*.fits)')
parser.add_argument('--gain', dest ='gain', default= 1.3, help = 'gain in e-/ADU.  default is 1.3, which applies to HDI camera.  Siena STL11000M gain = 0.8e-/ADU.')
parser.add_argument('--rdnoise', dest = 'rdnoise', default= 7.3, help = 'readnoise in e-.  default is 1.3, which applies to HDI camera.  Siena STL11000M = 11 e-')


#parser.add_argument('--', dest='pixelscalex', default='0.00011808', help='pixel scale in x (default = 0.00011808)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)
i=1
for f in files:
    print ('ZAPPING COSMIC RAYS FOR FILE %i OF %i'%(i,nfiles))
    with fits.open(f) as hdu1:
        print ('working on ',f)

        # convert data to CCDData format and save header
        ccd = CCDData(hdu1[0].data, unit=u.adu)
        header = hdu1[0].header
        crimage = ccdproc.cosmicray_lacosmic(ccd, gain = float(args.gain), readnoise = float(args.rdnoise))
        header['HISTORY'] = 'Cosmic rays rejected using ccdproc.cosmicray_lacosmic '
        fits.writeto('z'+f,crimage,header)
        hdu1.close()
    i += 1
    print ('\n')

    
    


                   
