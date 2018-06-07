#!/usr/bin/env python

'''
 
  GOAL:
  Make lists of files that contain
  - dome flats with same filter
  - sky flats with same filter
  Then combine and normalize the flats so they can be used to flatten image

  PROCEDURE:
  -User should move junk files to a junk subdirectory before starting
   - Junk files include initial bias frames pointing and other garbage frames
  - Use gethead to pull relevant header data
  - Overscan subtract and trim images (we assume these image names begin with 't  r')
  - Combine flats according to flat type (dome vs sky) and filter

  EXAMPLE:
     In the directory containing all flats type in the command line:
     '/home/share/research/pythonCode/uat_HDIgroupflatfiles.py'(or whatever the      path is to where this program is stored)

  
  INPUT/OUTPUT:
  Input: 'skyflat(type of filter)' or 'domeflat(type of filter)'
         -These contain trimmed image files grouped by filter grabbed from the i         mage headers seen at the beginning of this code. 
  Output: For combined flats --> 'cskyflat(type of filter).fits' or 'cdomeflat(t          ype of filter).fits'.
          For normalized flats --> 'nskyflat(type of filter).fits' or 'ndomeflat          (type of filter).fits'.

  REQUIRED MODULES:
  pyraf

  NOTES:
  in junkfile ftr flats still show. We changed the gethead requirements to only   bring in files that start with tr but the ftr files will not go away! =(
  
  WRITTEN BY:
  Rose A. Finn
  EDITED BY:
  Natasha Collova, Tiffany Flood, and Kaitlyn Hoag 5/29/15
  Grant Boughton, Natasha Collova, and Sandy Spicer 6/3/16
  UPDATES:
  Now combines and normalizes the grouped flat files 
    Input ex. = skyflatR
    Output ex. = cskyflatR.fits (combined sky flats in R band) & nskyflatR.fits     (normalized sky flats in R band)
  
'''
import glob
import os
import numpy as np

import ccdproc
from astropy.io import fits
import argparse

'''
from pyraf import iraf
iraf.noao()
iraf.imred()
iraf.ccdred()    
'''

parser = argparse.ArgumentParser(description ='Groups images by filter and creates flatfield images')
parser.add_argument('--filestring', dest='filestring', default='ztr', help='match string for input files (default =  ztr, which gets ztr*.fits)')
#parser.add_argument('--', dest='pixelscalex', default='0.00011808', help='pixel scale in x (default = 0.00011808)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)

os.system('gethead '+args.filestring+'*f00.fits CMMTOBS > tempflats')
# tempflats is the name of a "junk file" that contains the gethead information from all the flat images.
# This file will be deleted after the information is read out in the future.

# We assume that the flat images are trimmed and the file name starts with 'ztr'
infile=open('tempflats','r')
fnames=[]
filter=[]
ftype=[]   #skyflat or domeflat

for line in infile:
    t=line.split()
    fnames.append(t[0])
    ftype.append(t[1]+t[2])
    filter.append(t[3].rstrip('\n'))
infile.close()
set_filter=set(filter)
set_ftype=set(ftype)
array_ftype=np.array(ftype)
array_filter=np.array(filter)

for f in set_ftype:
    print "flat type=",f
    for element in set_filter:
        ftype_filter = str(f)+str(element)
        print 'grouping files for filter type = ',element
        indices=np.where((array_ftype == f)&(array_filter == element))
        if len(indices[0]) > 0:
            outfile = open(ftype_filter,'w')
            for i in indices[0]:
                outfile.write(fnames[i]+'\n')
            outfile.close()
os.remove('tempflats')            
flats = glob.glob('*flat*')
flats=set(flats)-set(glob.glob('c*.fits'))-set(glob.glob('n*.fits'))     # doesn't include flat files that have already been combined or normalized in the following loop
for f in flats:
    flatimages = []
    filelist = open(f,'r')
    for fname in filelist:
        fname = fname.rstrip()
        data,header = fits.getdata(fname, header=True)
        data = data / np.median(data)
        flatimages.append(data)
    # combine flat images using median combine
    med_flat = np.median(flatimages, axis=0)
    # normalize flat image by dividing by mean
    norm_med_flat = med_flat / np.mean(med_flat)
    header['HISTORY'] = 'Combined and normalized flat field'
    fits.writeto('nc'+f,norm_med_flat,header)


                

