#!/usr/bin/env python

'''
BASIC INFORMATION ABOUT THIS CODE:
  -This is the SIXTH program you need to run in order to perform the data reduct  ion of HDI images. 
  -This program has SCamp compute the astrometric solution and then has SWarp c   reate image mosaics. 


BEFORE RUNNNING THIS CODE:
  -Be sure to type ur_setup in the terminal everytime you open a new terminal window so that ureka is activated.
  -Ensure that pyraf is still activated by retyping the commands listed in the c  omments of the FIRST program titled "uat_HDIgroupflatfiles.py".


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
parser.add_argument('--filestring', dest = 'filestring', default = 'hcftr*o00.fits', help = 'string to use to get input files (default = "hcftr*o00.fits")')
parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Run sextractor to create object catalogs')
parser.add_argument('--scamp', dest = 'scamp', default = False, action = 'store_true', help = 'Run scamp')
parser.add_argument('--swarp', dest = 'swarp', default = False, action = 'store_true', help = 'Run swarp to create mosaics')
parser.add_argument('--l', dest = 'l', default = False, help = 'List of images to input to swarp')
parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of default config files')
parser.add_argument('--refimage',dest = 'refimage', default = None,  help = 'use a reference image to set center and size of output mosaic')
parser.add_argument('--m',dest = 'm', default = False, action = 'store_true', help = 'set if running for mosaic data')

args = parser.parse_args()

print 'testing ',args.refimage
if args.refimage:
    print args.refimage,' has a value of True'

# get input files
#print 'cp ' +args.d + '/default.* .'
os.system('cp ' +args.d + '/default.* .')
files = sorted(glob.glob(args.filestring))

nfiles = len(files)
i = 1


if args.s:
    for f in files:
        #read_exptime = 'gethead ' + f + ' EXPTIME'
        #exptime = subprocess.check_output(read_exptime,shell=True)
        #exptime = exptime.rstrip()
        #print f, exptime
        #if float(exptime) > 61.:
        print 'RUNNING SEXTRACTOR ON FILE %i OF %i'%(i,nfiles)
        t = f.split('.fits')
        froot = t[0]
        os.system('sex ' + f + ' -c default.sex.hdi -CATALOG_NAME ' + froot + '.cat')
        os.rename('check.fits', froot + 'check.fits')
        i += 1
        
if args.scamp:
    os.system('ls *.cat > scamp_input_cats')
    print 'RUNNING SCAMP'
    os.system('scamp @scamp_input_cats -c default.scamp')
    print 'DONE'
    
if args.swarp:
    #HDI data
    pixel_scale = 0.425
    if args.refimage:
        data,header = fits.getdata(args.refimage,header=True)
        w = WCS(header)
        image_size = data.shape

        ra,dec = w.wcs_pix2world(image_size[0]/2.,image_size[1]/2.,1)
        center = str(ra)+','+str(dec)
        mosaic_image_size = str(image_size[1])+','+str(image_size[0])
        print 'output mosaic image size = ',mosaic_image_size
        print 'center of mosaic = ',center     
    if not(args.l):
        print 'No file list provided for swarp'
    else:
        
        print 'RUNNING SWARP'
        if args.m:
            infile = open(args.l,'r')
            outfile = open('masks','w')
            for line in infile:
                t=line.split('.fits')
                outfile.write(t[0]+'_bpm.pl \n')
            outfile.close
            os.system('swarp @' + args.l + ' -c default.swarp -IMAGEOUT_NAME ' + args.l + '.coadd.fits -WEIGHTOUT_NAME ' + args.l + '.coadd.weight.fits -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE @masks')
        else:
            # CENTER_TYPE MANUAL
            #CENTER RA,DEC
            #PIXEL_SCALE 0.425
            if args.refimage:
                print 'using reference image with swarp'
                os.system('swarp @' + args.l + ' -c default.swarp -IMAGEOUT_NAME ' + args.l + '.coadd.fits -WEIGHTOUT_NAME ' + args.l + '.coadd.weight.fits -CENTER_TYPE MANUAL -CENTER '+center+' -PIXEL_SCALE '+str(pixel_scale)+' -IMAGE_SIZE '+mosaic_image_size)
            else:
                os.system('swarp @' + args.l + ' -c default.swarp -IMAGEOUT_NAME ' + args.l + '.coadd.fits -WEIGHTOUT_NAME ' + args.l + '.coadd.weight.fits ')
        print 'DONE'
         
