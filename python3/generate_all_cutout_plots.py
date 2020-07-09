#!/usr/bin/env python

"""
GOAL:
- generate cutouts for images with halpha cutouts

NOTES:
- run in the cutouts directory

- it will go into each subdir and run the script to generate cutouts!

"""
import glob
import os

import sys
homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/HalphaImaging/python3/')

dirlist = glob.glob('VF*')

for d in dirlist:
    os.chdir(d)
    imlist = glob.glob(d+'*-R.fits')
    if len(imlist) == 0:
        imlist = glob.glob(d+'*-r.fits')
    if len(imlist) == 0:
        continue
    for im in imlist:
        os.system('python '+homedir+'/github/HalphaImaging/python3/plot_cutouts_ha.py --r '+im+' --plotall')
        
        
    os.chdir('..')