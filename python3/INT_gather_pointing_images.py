#!/usr/bin/env python

'''

PROCEDURE:
* get list of current file
* combine images from the same pointing into one directory
* I had already split by pointing and filter, but I want to run scamp and swarp in the directory that contains all images from a pointing


'''

import os
import shutil
from astropy.io import fits
import glob

# get list of current directory
flist1 = os.listdir()
flist1.sort()

# get list of pointings - cut off filter from dir name
pointings = []
for f in flist1:
    if (f.find('pointing') > -1) & (f.find('-') > -1):
        pointings.append(f.split('-')[0])

# pointings will now have two listing for each filter b/c of r and Halpha
upointings = set(pointings) # unique pointings

subdirs = ['plots','shortexposure']
hafilters = ['Halpha','Ha6657']
for p in upointings:
    print(p)
    if not(os.path.exists(p)):
        os.mkdir(p)
    if os.path.exists(p+'-r'):
        os.system('mv '+p+'-r/* '+p+'/.')
        os.system('rmdir '+p+'-r')
    for h in hafilters:
        hadir = p+'-'+h
        if os.path.exists(hadir):
            os.system('mv {}/* {}/.'.format(hadir,p))
            for s in subdirs:
                subdirpath = '{}-{}/{}'.format(p,h,s)
                if os.path.exists(subdirpath):
                    command = 'mv {}/* {}/{}/.'.format(subdirpath,p,s)
                    print(command)
                    os.system(command)
                os.system('rmdir '+subdirpath)
            os.system('rmdir {}'.format(hadir))




