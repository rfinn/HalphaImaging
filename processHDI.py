#!/usr/bin/env python

'''
this is a wrapper script to run the many programs
that we use to reduce the HDI data.

the user should create a junkfile in each directory that
contains a list of files that should not be included in the
reduction (pointing check, focus, saturated sky flat, guided failed, etc.)



'''

import os
import sys

gitpath = os.getenv('HOME')+'/github/HalphaImaging/python3/'
sys.path.append('~/github/HalphaImaging/')

current_dir = os.getcwd()
# trim and overscan correct

os.system('python '+gitpath+'uat_overscantrim.py --filestring c')
# zap cosmic rays
os.system('python '+gitpath+'uat_zapcosmicrays.py --filestring tr')

# group flat files
os.system('python '+gitpath+'uat_HDIgroupflatfiles.py --filestring ztr')

# flatten science frames with dome flats
os.system('python '+gitpath+'uat_HDIflattendwithdome.py --filestring ztr')

# fix the HDI header
os.system('python '+gitpath+' uat_HDIfixheader.py --filestring d')

# run sextractor to create object lists
os.system('python '+gitpath+'uat_astr_mosaic.py --s --filestring h')

# run scamp to solve for astrometry
os.system('python '+gitpath+'uat_astr_mosaic.py --scamp --filestring h')

# sort objects by field
os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring h')

# run swarp?
#os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring h')
