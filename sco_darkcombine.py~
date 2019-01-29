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
parser.add_argument('--filestring', dest='filestring', default='dark', help='match string for input files (default =  dark, which gets dark*.fits)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)

os.system('gethead '+args.filestring+'*.fits EXPTIME > tempdarks')

# tempflats is the name of a "junk file" that contains the gethead information from all the flat images.
# This file will be deleted after the information is read out in the future.

# We assume that the flat images are trimmed and the file name starts with 'ztr'
infile=open('tempdarks','r')
fnames=[]
exptime=[]
ftype=[]


for line in infile:
    t=line.split()
    fnames.append(t[0])
    s = t[0].split('-')
    s2 = s[1].split('.')[0]
    ftype.append(s2[4:].replace('d','')) # pull out exptime number
    exptime.append(t[1])
infile.close()
set_ftype=set(ftype)
set_exptime=set(exptime)
array_ftype=np.array(ftype)
array_exptime=np.array(exptime)

filelist = []
for f in set_ftype:
    print "dark exposure time=",f
    ftype_group = 'dark'+str(f)
    filelist.append(ftype_group)
    indices=np.where((array_ftype == f))
    if len(indices[0]) > 0:
        outfile = open(ftype_group,'w')
        for i in indices[0]:
            outfile.write(fnames[i]+'\n')
        outfile.close()
#os.remove('tempflats')            
for f in filelist:
    print 'filelist = ',f
    fimages = []

    try:
        filelist = open(f,'r')
    except IOError:
        print('Problem opening file ',f)
        print('Hope that is ok...')
        continue
    for q in filelist: fimages.append(q.rstrip())
        
    # combine flat images using average combine, scale by median, sigma clip
    dark = ccdproc.combine(fimages,method='average',sigma_clip=True,unit=u.adu)
    print 'writing fits'
    dark.write(f+'.fits',overwrite=True)
    



                

