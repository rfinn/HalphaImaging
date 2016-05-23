#!/usr/bin/env python

'''
USE:
Changes certain header keywords to the standard.
Updates CMMTOBS --> FILTER
        RASTRNG --> CRVAL1
        DECSTRNG --> CRVAL2
Adds CRPIX1, CRPIX2, CD1_1, CD2_2, CTYPE1, CTYPE2

INPUT: all ftr*.fits in directory
OUTPUT: hftr*.fits



'''

import argparse
import glob
import os
import subprocess

parser = argparse.ArgumentParser(description ='Edit image headers to include basic WCS information to the HDI image headers')
parser.add_argument('--filestring', dest='filestring', default='ftr*.fits', help='match string for input files (default =  ftr*.fits)')
parser.add_argument('--minexptime', dest='minexptime', default=61., help='min exposure time of science frames (default = 61 sec)')
parser.add_argument('--subdir', dest='subdir', default='short_exptime', help='subdirectory for short exptime images (default = short_exptime)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring))
nfiles=len(files)
i=1
try:
    os.mkdir(args.subdir)
except OSError:
    print args.subdir,' already exists'
print 'subdirectory for short exposure time images is ',args.subdir
for f in files:
    print 'CHECKING EXPTIME FOR FILE %i OF %i'%(i,nfiles)
    read_exptime = 'gethead ' + f + ' EXPTIME'
    exptime = subprocess.check_output(read_exptime,shell=True)
    exptime = exptime.rstrip()
    if float(exptime) < args.minexptime:
        os.rename(f,args.subdir+'/'+f)
    i += 1
    print '\n'

    
    


                   
