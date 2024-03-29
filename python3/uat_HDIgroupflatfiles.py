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
  astropy

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
import sys
import numpy as np

import ccdproc
from astropy.io import fits
import astropy.units as u
import argparse

'''
from pyraf import iraf
iraf.noao()
iraf.imred()
iraf.ccdred()    
'''

parser = argparse.ArgumentParser(description ='Groups images by filter and creates flatfield images')
parser.add_argument('--filestring', dest='filestring', default='ztr', help='match string for input files (default =  ztr, which gets ztr*.fits)')
parser.add_argument('--siena', dest='siena', default=False,action='store_true', help='set this if reducing data from Siena STL11000M CCD')
parser.add_argument('--verbose', dest='verbose', default=False,action='store_true', help='print extra messages for troubleshooting')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)

if args.siena:
    print('running on Siena data - woo hoo!')
    os.system('gethead '+args.filestring+'*.fits IMAGETYP FILTER > tempflats')
else:
    print('running on HDI data - woo hoo!')
    print(('gethead '+args.filestring+'*f00.fits CMMTOBS > tempflats'))
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
    if args.verbose:
        print(f"CMMTOBS = {line}")
    fnames.append(t[0])
    if "sky" in line:
        ftype.append("skyflat")
    elif "dome" in line:
        ftype.append("domeflat")
    else:
        ftype.append('domeflat')


    try:
        if (line.find('6620') > -1) | (line.find('ha4') > -1) | (line.find('Ha4') > -1) :
            filter.append('ha4')

        elif (line.find('6660') > -1) |(line.find('ha8') > -1) | (line.find('Ha8') > -1) :
            filter.append('ha8')
        elif (line.find('6700') > -1) |(line.find('ha12') > -1)| (line.find('Ha12') > -1) :
            filter.append('ha12')
        elif (line.find('6740') > -1) |(line.find('ha16') > -1) | (line.find('Ha16') > -1) :
            filter.append('ha16')
        elif line.find('R') > -1:
            filter.append('R')
        #elif line.find('r') > -1: # doesn't work for 2014 data - ugh!!!, also 'r' is in the filename, so need to remove filename from line
        #    filter.append('r')
        elif line.find('V') > -1:
            filter.append('V')
        elif line.find('B') > -1:
            filter.append('B')
        
        #if len(t) == 2:# adding some additional cases to accomodate 2014-04-25 what has only filter
        #    filter.append(t[1].rstrip('\n'))
        #    
        elif len(t)> 2:
        #    else:
            print('problem with determing filter!!!')
            print('probably got a multi-word entry for CMMTOBS')
            print("I'm storing the second word...")
            print('filter = ',t[4].rstrip('\n'))
            filter.append(t[4].rstrip('\n'))
        else:
            filter.append(t[3].rstrip('\n'))
    except:
        print("Problem getting filter from CMMTOBS = ",line)
        sys.exit()
    if args.verbose:
        print(f"filter = {filter[-1]}, ftype = {ftype[-1]}")
infile.close()
set_filter=set(filter)
set_ftype=set(ftype)
array_ftype=np.array(ftype)
array_filter=np.array(filter)


# create files that contain all flats taken in same filter
flat_filelist = []
for f in set_ftype:
    print('####################################')
    print("flat type=",f)
    print('####################################')
    for element in set_filter:
        if args.verbose:
            print(f"working on {f}, {element}")
        ftype_filter = str(f)+str(element)
        flat_filelist.append(ftype_filter)
        print('grouping files for filter type = ',element)
        indices=np.where((array_ftype == f)&(array_filter == element))
        if len(indices[0]) > 0:
            outfile = open(ftype_filter,'w')
            for i in indices[0]:
                outfile.write(fnames[i]+'\n')
            outfile.close()
            if args.verbose:
                print(f"\ngot {len(indices[0])} in {f},{element}\n")

# for testing
#sys.exit()

for f in flat_filelist:
    print('filelist = ',f)
    flatimages = []
    try:
        filelist = open(f,'r')
    except IOError:
        print(('Problem opening file ',f))
        print('Hope that is ok...')
        continue
    for q in filelist: flatimages.append(q.rstrip())
    if args.verbose:
        print(f"flatimages = {flatimages}")
    if len(flatimages) < 3:
        print('problem combining images from ',f)
        continue
    if args.verbose:
        print(f"\nrunning ccdproc.combine on {f}\n")
        
    # combine flat images using average combine, scale by median, sigma clip
    flat = ccdproc.combine(flatimages,scale=np.median,method='average',sigma_clip=True,unit=u.adu)
    #med_flat = ccdproc.combine(flatimages, method='median')
    # normalize flat image by dividing by mean

    # adding the np.array part to get rid of error on
    # "MaskedArray.tofile() not implemented yet."
    norm_flat = np.array(flat/ np.mean(flat))
    print(f"writing fits {'n'+f+'.fits'}")
    fits.writeto('n'+f+'.fits', norm_flat, overwrite=True)




                
# clean up
os.remove('tempflats')            
