#!/usr/bin/env python
'''
GOAL
update image header with information in .head file that is produced by scamp

RATIONALE
I am creating this code so that I can run getzp.py on images of photometric standards.  typically, we are taking only one image of these standards, so we don't want to create a mosaic using swarp.  We do want the astrometric solution found by scamp though, and so I want to merge the original image with the header information produced by scamp.

INPUT
image prefix


OUTPUT
prepend 's' to image name to signify 'scamp header'

'''
from astropy.io import fits
import argparse
import glob
import os

parser = argparse.ArgumentParser(description ='Edit image headers to include header produced by scamp')
parser.add_argument('--filestring', dest='filestring', default='hd', help='match string for input files (default =  hd, which operates on all hd*.fits files)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)
i=1
for f in files:

    prefix = f.split('.fits')[0]
    if not(os.path.exists(prefix+'.head')):
        print('no header file found for ',f)
        print('moving to next image')
        continue
    im, header = fits.getdata(f,header=True)
    t = open(prefix+'.head','r')
    scamp = t.readlines()
    t.close()
    for h in scamp:
        #print(h)
        if h.find('HISTORY') > -1:
            header.set('HISTORY',h.rstrip().split('HISTORY')[1])
        elif h.find('COMMENT') > -1:
            header.set('COMMENT',h.rstrip().split('COMMENT')[1])
        elif h.find('END') > -1:
            #print('reached end of header file for ',f)
            print('writing updated header for ',f)
            print('')
            break
        else:
            key,junk = h.split('=')
            key = key.strip()
            #print(junk)
            test = junk.rstrip().split('/')
            if len(test) > 2:
                val = test[0]
                unit = test[1]+'/'+test[2] 
            elif len(test) == 2:
                val, unit = junk.rstrip().split('/')
            #print(key,val,unit)
            try:
                val = float(val)

            except ValueError:
                #print('not a floating point')
                val = val.strip().strip("'")
            #try:
            #    print('BEFORE: ',key,header[key])
            #except KeyError:
            #    print('adding new keyword ',key)
            header.set(key,val, unit)
            #print('AFTER:  ',header[key])
    fits.writeto('s'+f, im, header, overwrite=True)
