#!/usr/bin/env python

'''
GOAL:
  The goal of this program is to have scamp compute the astrometric 
  solution and then have swarp create image mosaics.
  This assumes that images have been reduced through flatfielding.
  If using HDI data, you need to run a separate program to add convert
  some header keywords to standard values and to add a rough WCS solution.

PROCEDURE:
  - copy default setup files to current directory
  - Run sextractor on each image
  - Make a list containing all .cat files
  - Run scamp
  - Run swarp

EXAMPLE:
   In the directory containing all flattened objects with fixed headers to run sextractor type in the command line:
      '/Users/alfalfa/Github/HalphaImaging/uat_astr_mosaic.py --s'(or whatever the path is to where this program is stored)

   To get swarp to create aligned images in multiple bans (e.g. Halpha and R-band), do the following
    uat_astr_mosaic.py --swarp --l A1367-h02_R

    uat_astr_mosaic.py --swarp --l A1367-h02_ha12 --refimage 'A1367-h02_R.coadd.fits'

    uat_astr_mosaic.py --swarp --l A1367-h02_R --refimage 'A1367-h02_R.coadd.fits'



WHAT THIS CODE DOES:
INPUT/OUPUT:
REQUIRED MODULES:
EXTRA NOTES:
WRITTEN BY:
Rose Finn

EDITED BY:
Research Team 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood, Kaitlyn Hoag, and Kelly Whalen
'''
import glob
import os
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import subprocess

parser = argparse.ArgumentParser(description ='Run sextractor, scamp, and swarp to determine WCS solution and make mosaics')
parser.add_argument('--filestring', dest = 'filestring', default = 'h', help = 'string to use to get input files (default = "h", which grabs all of the files "h*o00.fits")')
parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Run sextractor to create object catalogs')
parser.add_argument('--scamp', dest = 'scamp', default = False, action = 'store_true', help = 'Run scamp')

parser.add_argument('--siena', dest = 'siena', default = False, action = 'store_true', help = 'set this when running Siena data')
parser.add_argument('--pisces', dest = 'pisces', default = False, action = 'store_true', help = 'set this when running PISCES data')
parser.add_argument('--swarp', dest = 'swarp', default = False, action = 'store_true', help = 'Run swarp to create mosaics')
parser.add_argument('--l', dest = 'l', default = False, help = 'List of images to input to swarp')
parser.add_argument('--d',dest = 'd', default ='~/github/HalphaImaging/astromatic', help = 'Locates path of default config files.  Default is ~/github/HalphaImaging/astromatic')
parser.add_argument('--refimage',dest = 'refimage', default = None,  help = 'use a reference image to set center and size of output mosaic')
parser.add_argument('--m',dest = 'm', default = False, action = 'store_true', help = 'set if running for mosaic data')
parser.add_argument('--int',dest = 'int', default = False, action = 'store_true', help = 'set if running on INT data')
parser.add_argument('--noback',dest = 'noback', default = False, action = 'store_true', help = 'set to disable background subtraction in swarp')

args = parser.parse_args()

print(('testing ',args.refimage))
if args.refimage:
    print((args.refimage,' has a value of True'))

# get input files
#print 'cp ' +args.d + '/default.* .'
os.system('cp ' +args.d + '/default.* .')
files = sorted(glob.glob(args.filestring+"*.fit*"))

nfiles = len(files)
print(('number of files = ',nfiles))
i = 1


if args.s:
    for f in files:
        #read_exptime = 'gethead ' + f + ' EXPTIME'
        #exptime = subprocess.check_output(read_exptime,shell=True)
        #exptime = exptime.rstrip()
        #print f, exptime
        #if float(exptime) > 61.:
        print(('RUNNING SEXTRACTOR ON FILE %i OF %i'%(i,nfiles)))
        t = f.split('.fits')
        froot = t[0]
        if args.int:
            os.system('sex ' + f + ' -c default.sex.INT -CATALOG_NAME ' + froot + '.cat')
        elif args.siena:
            os.system('sex ' + f + ' -c default.sex.siena -CATALOG_NAME ' + froot + '.cat')
        elif args.pisces:
            os.system('sex ' + f + ' -c default.sex.siena -CATALOG_NAME ' + froot + '.cat -DETECT_MINAREA 10 -DETECT_THRESH 2 -ANALYSIS_THRESH 2 -BACK_SIZE 128 -GAIN 100')
            print('HEY! PISCES!!!!!')
        else:
            os.system('sex ' + f + ' -c default.sex.HDI -CATALOG_NAME ' + froot + '.cat')
        #os.rename('check.fits', froot + 'check.fits')
            
        i += 1
        
if args.scamp:
    os.system('ls '+args.filestring+'*.cat > scamp_input_cats')
    print('RUNNING SCAMP')
    if args.siena:
        os.system('scamp @scamp_input_cats -c default.scamp.siena')
    if args.int:
        os.system('scamp @scamp_input_cats -c default.scamp.INT')
    elif args.pisces:
        os.system('scamp @scamp_input_cats -c default.scamp.pisces ')
    else:
        os.system('scamp @scamp_input_cats -c default.scamp')
    print('DONE')
    
if args.swarp:
    #HDI data
    pixel_scale = 0.425
    defaultswarp = 'default.swarp'
    if args.int:
        pixel_scale = 0.331
        defaultswarp = 'default.swarp.INT'
    if args.noback:
        outimage = '.noback.coadd.fits'
        weightimage = '.noback.coadd.weight.fits'        
    else:
        outimage = '.coadd.fits'
        weightimage = '.coadd.weight.fits'                
    if args.refimage:
        data,header = fits.getdata(args.refimage,header=True)
        w = WCS(header)
        image_size = data.shape

        ra,dec = w.wcs_pix2world(image_size[0]/2.,image_size[1]/2.,1)
        center = str(ra)+','+str(dec)
        mosaic_image_size = str(image_size[1])+','+str(image_size[0])
        #print(('output mosaic image size = ',mosaic_image_size))
        #print(('center of mosaic = ',center))
    if not(args.l):
        print('No file list provided for swarp')
    else:
        
        print('RUNNING SWARP')
        if args.m: 
            infile = open(args.l,'r')
            outfile = open('masks','w')
            for line in infile:
                t=line.split('.fits')
                outfile.write(t[0]+'_bpm.pl \n')
            outfile.close
            os.system('swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage+' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE @masks')
        else:
            # CENTER_TYPE MANUAL
            #CENTER RA,DEC
            #PIXEL_SCALE 0.425
            if args.refimage:
                print('using reference image with swarp')
                commandstring = 'swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage+' -CENTER_TYPE MANUAL -CENTER '+center+' -PIXEL_SCALE '+str(pixel_scale)+' -IMAGE_SIZE '+mosaic_image_size 
            else:
                commandstring='swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage
            if args.noback:
                commandstring = commandstring+' -SUBTRACT_BACK N'
            print('running the following command:')
            print(commandstring)
            os.system(commandstring)

        print('DONE')

