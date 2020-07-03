#!/usr/bin/env python

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
    for im in inlist:
        os.system('python '+homedir+'/github/HalphaImaging/python3/plot_cutouts_ha.py --r '+im+' --plotall')
    os.chdir('..')
