#!/usr/bin/env python

'''
This is a wrapper script to run the many programs
that we use to reduce the HDI data.

Before running this program, 
the user should create a junkfile in each directory that
contains a list of files that should not be included in the
reduction (pointing check, focus, saturated sky flat, guided failed, etc.)

USAGE:

processHDI.py --trim --zap --groupflat --flatwdome --fixheader

python ~/github/HalphaImaging/processHDI.py --trim --zap --groupflat --flatwdome --fixheader

'''

import os
import sys
import argparse
import glob

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
parser.add_argument('--filestring', dest = 'filestring', default = 'h', help = 'string to use to get input files (default = "h", which grabs all of the files "h*o00.fits")')

parser.add_argument('--se', dest='se', default=False,action='store_true', help='run sextractor.  best to do this on all files from an observing run at once.')
parser.add_argument('--scamp', dest='scamp', default=False,action='store_true', help='run scamp, and group files.  best to do this on all files from an observing run at once.')
parser.add_argument('--submed', dest='submed', default=False,action='store_true', help='run subtract_median to remove the background before running swarp.')
parser.add_argument('--swarp', dest='swarp', default=False,action='store_true', help='run swarp to make coadded images.')
parser.add_argument('--filelist', dest='filelist', default='swarp_input', help='list of image sets to run swarp on.  the file should contain the list of all Rband groups, for example: ls pointing*_R > swarp_input.  This will look for the corresponding list of halpha images.')
parser.add_argument('--zp', dest='zp', default=False,action='store_true', help='run getzp.py on all *coadd.fits images')
parser.add_argument('--uat', dest='uat', default=False,action='store_true', help='set when running on UAT groups data that uses a bunch of different halpha filters')
args = parser.parse_args()

trim = args.trim
zap = args.zap
group_flat = args.groupflat
dflat = args.flatwdome
fixheader=args.fixheader
#astr = args.astr

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
    os.system('python '+gitpath+'uat_HDIgroupflatfiles.py --filestring ztr --verbose')

# flatten science frames with dome flats
if dflat:
    os.system('python '+gitpath+'uat_HDIflattenwithdome.py')
    mylist = ['ZAPPED','z']
    if not(os.path.exists(mylist[0])):
        os.mkdir(mylist[0])
        os.system('mv '+mylist[1]+'*.fits '+mylist[0]+'/.')


# fix the HDI header
if fixheader:
    if args.uat:
        os.system('python '+gitpath+'uat_HDIfixheader.py --filestring d --uat')
    else:
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
if args.se:
    # run sextractor to create object lists
    os.system('python '+gitpath+'uat_astr_mosaic.py --s --filestring '+args.filestring)
if args.scamp:
    # run scamp to solve for astrometry
    os.system('python '+gitpath+'uat_astr_mosaic.py --scamp --filestring '+args.filestring)

    # sort objects by field
    os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring '+args.filestring)

if args.submed:
    # subtract median from images before running swarp
    os.system('python '+gitpath+'subtract_median.py --filestring '+args.filestring)

    # sort objects by field
    os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring '+args.filestring)


# run swarp?
#os.system('python '+gitpath+'uat_HDIsortobjects.py --filestring h')
if args.swarp:
    infile = open(args.filelist,'r')
    print('running swarp')
    i = 0
    for f in infile:
        f = f.rstrip()
        print('PROCESSING INPUT FILE: ',f)
        if f.find(r'_R') > -1:
            print('R filter')
            rootname = f.split('_R')[0]
        elif f.find(r'_r') > -1:
            print('sdss r filter')
            rootname = f.split('_r')[0]
        else:
            print('unexpected file name format')
            sys.exit()
        # get set of halpha images
        print('rootname = ',rootname)
        fnames = glob.glob(rootname+'_h*coadd.fits')
        print('halpha file = ',fnames)
        if len(fnames) > 1:
            print('got more than one Halpha image - crazy!')
            print('hope this is ok...')
            multiha = True
        else:
            multiha = False
            halist = fnames[0]
        # get name of R-band coadd
        rcoadd_image = f+'.noback.coadd.fits'
        # run swarp on r images
        print('python '+gitpath+'uat_astr_mosaic.py --swarp --l '+f)
        os.system('python '+gitpath+'uat_astr_mosaic.py --swarp --l '+f+' --noback')
        

        if multiha:
            for h in fnames:
                # run swarp on halpha, with r as reference image
                os.system('python '+gitpath+'uat_astr_mosaic.py --swarp --l '+h+' --refimage '+rcoadd_image+' --noback')
        else:
            # run swarp on r, with r as reference image
            os.system('python '+gitpath+'uat_astr_mosaic.py --swarp --l '+halist+' --refimage '+rcoadd_image+' --noback')
        # run swarp on r, with r as reference image
        os.system('python '+gitpath+'uat_astr_mosaic.py --swarp --l '+f+' --refimage '+rcoadd_image+' --noback')
        #break
        
    infile.close()

# solve for photometric zp
if args.zp:
    filelist = glob.glob('*coadd.fits')
    for f in filelist:
        if f.find('_h') > -1:
            photfilter = 'r'
        elif f.find('-h') > -1:
            photfilter = 'r'
        elif f.find('_r') > -1:
            photfilter = 'r'
        elif f.find('-r') > -1:
            photfilter = 'r'
        elif f.find('_R') > -1:
            photfilter = 'R'
        elif f.find('-R') > -1:
            photfilter = 'R'
        os.system('python '+gitpath+'getzp.py --image '+f+' --filter '+photfilter+' --instrument i ')
        #break
