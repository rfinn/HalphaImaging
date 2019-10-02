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
import argparse

gitpath = os.getenv('HOME')+'/github/HalphaImaging/python3/'
sys.path.append('~/github/HalphaImaging/')

current_dir = os.getcwd()
# trim and overscan correct

parser = argparse.ArgumentParser(description ='Reduce HDI data, trim through scamp')



parser.add_argument('--trim', dest='trim', default=False,action='store_true', help='trim images')
parser.add_argument('--zap', dest='zap', default=False,action='store_true', help='run cosmic ray reject (this takes a long time)')
parser.add_argument('--groupflat', dest='groupflat', default=False,action='store_true', help='group files by filter and create normalized flats')
parser.add_argument('--flatwdome', dest='flatwdome', default=False,action='store_true', help='flatten images with dome flats')
parser.add_argument('--fixheader', dest='fixheader', default=False,action='store_true', help='fix headeer to get files ready for scamp.')
parser.add_argument('--astr', dest='astr', default=False,action='store_true', help='run sextractor, scamp, and group files.  best to do this on all files from an observing run at once.')
args = parser.parse_args()

trim = args.trim
zap = args.zap
group_flat = args.groupflat
dflat = args.flatwdome
fixheader=args.fixheader
astr = args.astr

if trim:
    os.system('python '+gitpath+'uat_overscantrim.py --filestring c')
mylist = ['ORIGINALS','c']
if not(os.path.exists(mylist[0])):
    os.mkdir(mylist[0])
    os.system('mv '+mylist[1]+'*.fits '+mylist[0]+'/.')

# zap cosmic rays
if zap:
    os.system('python '+gitpath+'uat_zapcosmicrays.py --filestring tr')
    mylist = ['TRIMMED','tr']
    if not(os.path.exists(mylist[0])):
        os.mkdir(mylist[0])
        os.system('mv '+mylist[1]+'*.fits '+mylist[0]+'/.')

# group flat files
if group_flat:
    os.system('python '+gitpath+'uat_HDIgroupflatfiles.py --filestring ztr')

# flatten science frames with dome flats
if dflat:
    os.system('python '+gitpath+'uat_HDIflattenwithdome.py')
    mylist = ['ZAPPED','z']
    if not(os.path.exists(mylist[0])):
        os.mkdir(mylist[0])
        os.system('mv '+mylist[1]+'*.fits '+mylist[0]+'/.')


# fix the HDI header
if fixheader:
    os.system('python '+gitpath+'uat_HDIfixheader.py --filestring d')
    mylist = ['FLATTENED','d']
    if not(os.path.exists(mylist[0])):
        os.mkdir(mylist[0])
        os.system('mv '+mylist[1]+'*.fits '+mylist[0]+'/.')


# best approach is to create a directory of all flattened images
# from a given run.
# then run scamp in this folder.
# this will automatically scale images from different nights
# and you can swarp images with different exptime and
# images taken on different nights. :)
if astr:
    # run sextractor to create object lists
    os.system('python '+gitpath+'uat_astr_mosaic.py --s --filestring h')

    # run scamp to solve for astrometry
    os.system('python '+gitpath+'uat_astr_mosaic.py --scamp --filestring h')

    # sort objects by field
    os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring h')

# run swarp?
#os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring h')
