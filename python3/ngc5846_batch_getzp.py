#!/usr/bin/env python

'''
GOAL: 
- go into each subdirectory and run getzp on coadded images
- delete any intermediate images to save space

Run this from, /home/rfinn/data/reduced/NGC5846/
- this directory has a subdirectory for each pointing


'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib
# get list of current directory
flist1 = glob.glob('NGC*.fits')
# overwrite output files if they exist
overwrite = True
flist1.sort()
for f in flist1: # loop through list
    print('##########################################')
    print('##########################################')        
    print('WORKING ON IMAGE: ',f)
    print('##########################################')
    print('##########################################')
            
    try:
        if f.find('Ha') > -1:
            filter='ha'
        else:
            filter='r'
        os.system('python ~/github/HalphaImaging/python3/getzp.py --instrument h --nexptime --image {} --filter {}'.format(f,filter))
    except:
        print('##########################################')
        print('WARNING: problem running getzp for ',f)
        print('##########################################')

        


