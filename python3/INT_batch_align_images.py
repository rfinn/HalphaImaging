#!/usr/bin/env python

'''
GOAL: 
- go into each subdirectory and run getzp on coadded images
- delete any intermediate images to save space

Run this from, e.g. /home/rfinn/data/reduced/scratch-int-feb2019
- this directory has a subdirectory for each pointing


'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib
matplotlib.use("Qt5agg")
homedir = os.getenv("HOME")
telescope = 'INT'
# get list of current directory
flist1 = glob.glob('VF-*r.fits')
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True

# just doing this for rogue pointings 022 and 026
#flist1a = glob.glob('VF-*p022*r.fits')
#flist1b = glob.glob('VF-*p026*r.fits')
#flist1 = flist1a+flist1b
flist1.sort()
#print(flist1)


# shift r to match Halpha

for rimage in flist1: # loop through list


    print('##########################################')
    print('##########################################')        
    print('WORKING ON IMAGE: ',rimage)
    print('##########################################')
    print('##########################################')
    
    # read in r-band images
    # find matching Halpha image

    # grab other coadds
    rootname = rimage.split('-r')[0]
    rweightimage = rootname+'-r.weight.fits'
    coadds = glob.glob(rootname+'*.fits')
    #print(coadds)
    haimage = None
    for c in coadds:
        if (c.find('-Ha') > -1) & (c.find('weight') < 0):
            haimage = c
            print('matching ha image: ',haimage)
    if haimage is not None:
        command_string = 'python ~/github/HalphaImaging/python3/INT_align_images.py --image1 {} --image2 {} --weight2 {}'.format(haimage,rimage,rweightimage)
        # check to see if shifted r-band image exists.  if 
        try:
            print('running : ',command_string)
            os.system(command_string)
        except:
            print('##########################################')
            print('WARNING: problem running align_images for ',rimage)
            print('##########################################')

    # just running on one directory for testing purposes
    break



