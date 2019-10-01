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

trim = False
zap = False
group_flat = False
dflat = True
fixheader=True
astr = True

if trim:
    os.system('python '+gitpath+'uat_overscantrim.py --filestring c')
os.mkdir('ORIGINALS')
os.system('mv c*.fits ORIGINALS/.')

# zap cosmic rays
if zap:
    os.system('python '+gitpath+'uat_zapcosmicrays.py --filestring tr')
os.mkdir('TRIMMED')
os.system('mv tr*.fits TRIMMED/.')

# group flat files
if group_flat:
    os.system('python '+gitpath+'uat_HDIgroupflatfiles.py --filestring ztr')

# flatten science frames with dome flats
if dflat:
    os.system('python '+gitpath+'uat_HDIflattenwithdome.py --filestring ztr')
    os.mkdir('ZAPPED')
    os.system('mv z*.fits ZAPPED/.')


# fix the HDI header
if fixheader:
    os.system('python '+gitpath+' uat_HDIfixheader.py --filestring d')

# run sextractor to create object lists
if astr:
    os.system('python '+gitpath+'uat_astr_mosaic.py --s --filestring h')

    # run scamp to solve for astrometry
    os.system('python '+gitpath+'uat_astr_mosaic.py --scamp --filestring h')

    # sort objects by field
    os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring h')

# run swarp?
#os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring h')
