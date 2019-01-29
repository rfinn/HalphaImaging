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




# figure out which dark and flat to use for each file

# assume darks are in DARKS directory
# assume flats are in FLATS directory

# if calibration frames are named methodically, we should be able
# to construct the name of the flat and dark from the FILTER and EXPTIME


os.system('gethead '+args.filestring+'*.fits FILTER EXPTIME > temp-process-images')

files = open('temp-process-images','r')
i=1
nfiles = len(files)
for t in files:
    f, filt, exptime = t.split()
    
    print 'ZAPPING COSMIC RAYS FOR FILE %i OF %i'%(i,nfiles)


    
    with fits.open(f) as hdu1:
        print 'working on ',f
        outfile=f
        # convert data to CCDData format and save header
        ccd = CCDData(hdu1[0].data, unit=u.adu)
        header = hdu1[0].header

        if args.zap:
            
            zccd = ccdproc.cosmicray_lacosmic(ccd, gain = float(args.gain), readnoise = float(args.rdnoise))
            outfile = 'z'+outfile

        # subtract dark
        if args.dark:
            if args.zap:
                infile=zccd
            else:
                infile = ccd
            dccd = ccdproc.subtract_dark(infile, dark)
            outfile = 'd'+outfile

        # flatten
        if args.flat:
            # flatten
            if args.dark:
                infile = dccd
            else:
                if args.zap:
                    infile = zccd
                else:
                    infile = ccd
            fccd = ccdproc.flat_correct(infile, gain = float(args.gain), readnoise = float(args.rdnoise),add_keywork='flat_corrected')
            outfile = 'f'+outfile
            fccd.write(outfile)
    i += 1
    print '\n'

                   
