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

Summer 2025 - yes, that's 2025!

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
parser.add_argument('--instrument', dest='instrument', default='h', const='h', nargs='?', choices=['h', 'i', 'm'], help='instrument.  options are h=HDI, m=mosaic, i=INT')
parser.add_argument('--m',dest = 'm', default = False, action = 'store_true', help = 'set if running for mosaic data')
parser.add_argument('--int',dest = 'int', default = False, action = 'store_true', help = 'set if running on INT data')
parser.add_argument('--noback',dest = 'noback', default = False, action = 'store_true', help = 'set to disable background subtraction in swarp')
parser.add_argument('--psfex',dest = 'psfex', default = False, action = 'store_true', help = 'use psfex config files when running sextractor')

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

def update_coadd_header(coadd, input_image_list):
    '''
    GOAL:
    scamp is not passing the header keywords through to the coadded image,
    so the goal of this routine is to add the keywords manually.  This

    INPUT:
    * coadd - coadded image to add header keywords too
    * input_image_list - list of images to get header keywords from
    '''

    # read in coadd header and data
    hdu = fits.open(coadd)



    # write out coadd with new header
    header_fields = ['OBSDATE','AIRMASS','EXPTIME','MEDSUB','SKYMED','SKYSTD']# List of FITS keywords to propagate
    infile = open(input_image_list,'r')
    input_images = []
    for line in infile:
        input_images.append(line.rstrip())
    input_images.sort()
    print("\ngathering header info from these input images: \n",input_images)
    infile.close()

    
    mosaic_flag = False
    h1 = fits.getheader(input_images[0])
    try:
        instrument = h1['INSTRUME']
        if 'Mos' in instrument:
            mosaic_flag = True
            header_fields[0] = 'DATE-OBS'
            print(f"updating header fields for Mosaic to: {header_fields}")
    except:
        print(f"WARNING: could not get INSTRUME in the header of image {input_images[0]}")
    
    # add ccdnoise
    f = 'CCDNOISE'
    hdu[0].header.set(f,7.3)
    

    nimage = 1    
    for im in input_images:
        # record the image name
        hdu[0].header.set(f'INIMAG{nimage}',im)

        # 
        iheader = fits.getheader(im)
        # check if this is mosaic

        
        for f in header_fields:
            newfield = f"{f}{nimage}"
            if len(newfield) > 8:
                newfield = newfield.replace('-','')
            #print(newfield)

            # for int, we are getting 10s of images into a coadd, so we are over the 8 char limit
            # shortening the header field names in this case
            # AIRMASS -> ARMASS
            # EXPTIME -> EXPT
            if (len(newfield) > 8) & ('AIRMASS' in newfield): 
                newfield = newfield.replace('I','')
            elif (len(newfield) > 8) & ('EXPTIME' in newfield):
                newfield = newfield.replace('IME','')
               
            try:
                if f in ['MEDSUB','SKYMED','SKYSTD']:
                    hstring = f"{iheader[f]:.3f}"
                else:
                    hstring = iheader[f]
                hdu[0].header.set(newfield,hstring)
            except KeyError:
                print(f"WARNING: Keyword {f} not found")
        nimage += 1


    # write out image
    hdu.writeto(coadd,overwrite=True)

    
if args.s:
    for f in files:
        #read_exptime = 'gethead ' + f + ' EXPTIME'
        #exptime = subprocess.check_output(read_exptime,shell=True)
        #exptime = exptime.rstrip()
        #print f, exptime
        #if float(exptime) > 61.:

        # do not run source extractor if there is a 'weight' in the filename
        if f.find('.weight') > -1:
            continue

        # do not run source extractor if there is a 'bpm' in the filename
        if f.find('.bpm') > -1:
            continue
        
        print(('RUNNING SEXTRACTOR ON FILE %i OF %i'%(i,nfiles)))
        t = f.split('.fits')
        froot = t[0]
        extra_commands = ""
        if args.instrument == 'i':
            config_file = "default.sex.INT"
        elif args.instrument == 'h':
            config_file = "default.sex.HDI"
        elif args.instrument == 'm':
            config_file = "default.sex.MOS" 
        elif args.siena | args.pisces:
            config_file = "default.sex.siena"
        elif args.pisces:
            config_file = "default.sex.siena"
            extra_commands = " -DETECT_MINAREA 10 -DETECT_THRESH 2 -ANALYSIS_THRESH 2 -BACK_SIZE 128 -GAIN 100"
            print('HEY! PISCES!!!!!')
        else:
            config_file = "default.sex.HDI"

        if args.psfex:
            config_file = "default.sex.HDI.psfex"
            

        ## RUN SOURCE EXTRACTOR
        se_call = f"sex {f} -c {config_file} -CATALOG_NAME {froot}.cat {extra_commands}"

        print(f"running source extractor as: {se_call}")
        os.system(se_call)
        
        #os.rename('check.fits', froot + 'check.fits')
            
        i += 1
        
if args.scamp:
    os.system('ls '+args.filestring+'*.cat > scamp_input_cats')
    print('RUNNING SCAMP')
    if args.siena:
        os.system('scamp @scamp_input_cats -c default.scamp.siena')
    if args.instrument == 'i':
        os.system('scamp @scamp_input_cats -c default.scamp.INT')
    elif args.pisces:
        os.system('scamp @scamp_input_cats -c default.scamp.pisces ')
    else:
        os.system('scamp @scamp_input_cats -c default.scamp')
    print('DONE')
    
if args.swarp:
    #HDI data

    defaultswarp = 'default.swarp' # for HDI
    if args.instrument == 'i':
        pixel_scale = 0.331
        defaultswarp = 'default.swarp.INT'
    elif args.instrument == 'm':
        defaultswarp = 'default.swarp.MOS'
        pixel_scale = 0.4255        
    elif args.instrument == 'h':
        defaultswarp = 'default.swarp.HDI'
        pixel_scale = 0.4255

    # TODONE - do I need a case for BOK? No, I have a separate program for post-pipeline (NOAO pipeline) processing of BOK images: BOK_run_swarp.py
    if args.noback:
        outimage = '.noback.coadd.fits'
        weightimage = '.noback.coadd.weight.fits'        
    else:
        outimage = '.coadd.fits'
        weightimage = '.coadd.weight.fits'                
    if args.refimage:
        data,header = fits.getdata(args.refimage,header=True)
        refwcs = WCS(header)
        image_size = data.shape

        #ra,dec = refwcs.wcs_pix2world(image_size[0]/2.,image_size[1]/2.,1)
        # why not use CRVAL1 and CRVAL2 from reference image???
        ra = header['CRVAL1']
        dec = header['CRVAL2']

        center = str(ra)+','+str(dec)        


        
        mosaic_image_size = str(image_size[1])+','+str(image_size[0])

        # get pixel scale from the ref image rather than using something fixed for all images
        # it could be that the different halpha images will have different
        #pscale = wcs.utils.proj_plane_pixel_scales(refwcs) # in deg -> arcsec
        #pixel_scale = round(pscale[0]*3600,5) # in arcsec per pixel
        
        #print(('output mosaic image size = ',mosaic_image_size))
        #print(('center of mosaic = ',center))
    if not(args.l):
        print('No file list provided for swarp')
    else:
        
        print('RUNNING SWARP')
        # this is the basic command for running swarp
        commandstring = 'swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage +' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+str(pixel_scale)
        if args.m:

            # this create a list of bpm files that tell swarp to ignore the chip gaps, etc.
            infile = open(args.l,'r')
            outfile = open('masks','w')
            for line in infile:
                t=line.split('.fits')
                outfile.write(t[0]+'_bpm.pl \n')
            outfile.close
            # may have used this command for virgo mosaic images
            #os.system('swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage+' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE @masks')
            # updating for uat Halpha groups, mosaic data
            #os.system('swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage+' -WEIGHT_TYPE MAP_WEIGHT -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+str(pixel_scale))
            commandstring += ' -WEIGHT_TYPE MAP_WEIGHT'
        # CENTER_TYPE MANUAL
        #CENTER RA,DEC
        #PIXEL_SCALE 0.425
        if args.refimage:
            print('using reference image with swarp')
            #commandstring = 'swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage+' -CENTER_TYPE MANUAL -CENTER_TYPE MANUAL -CENTER '+center+' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+str(pixel_scale)+' -IMAGE_SIZE '+mosaic_image_size#+' -RESAMPLE N'
            commandstring += ' -CENTER_TYPE MANUAL -CENTER '+center+' -IMAGE_SIZE '+mosaic_image_size
        #else:
        #    commandstring='swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage +' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+str(pixel_scale)
        if args.noback:
            commandstring += ' -SUBTRACT_BACK N'
        print('running the following command:')
        print(commandstring)
        os.system(commandstring)


        # add function to add header keywords from individual images into the coadd name
        update_coadd_header(args.l + outimage,  args.l)
    
        print('DONE')

