#!/usr/bin/env python


'''

GOAL: create coadd images 

OVERVIEW:
* designed to complete final stages of coaddition for INT WFC images
* assumes processing is complete through astrometric correction.  
  - code will look for .head file

REFERENCES:
https://photutils.readthedocs.io/en/stable/background.html

https://reproject.readthedocs.io/en/stable/mosaicking.html

ccdproc wcs projection
https://ccdproc.readthedocs.io/en/latest/image_combination.html
'''

#from photutils import make_source_mask
import os
homedir = os.getenv("HOME")
import sys
print(homedir)

sys.path.append('/home/rfinn/github/HalphaImaging/python3/')

#import imutils
import ha_imutils as imutils

from reproject import reproject_interp
#from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs

import ccdproc as ccdp
#from ccdproc import wcs_project
#from ccdproc import Combiner

from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import mad_std

import numpy as np

##########################################################
### SET UP MASKS
##########################################################

# read in WFC ccd masks
mask1 = CCDData.read(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits',unit='adu')
mask2 = CCDData.read(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits',unit='adu')
mask3 = CCDData.read(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits',unit='adu')
mask4 = CCDData.read(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits',unit='adu')
#mask1.data = mask1.data.astype('bool')
#mask2.data = mask2.data.astype('bool')
#mask3.data = mask3.data.astype('bool')
#mask4.data = mask4.data.astype('bool')

# dictionary to translate between instrument in header and mask file
maskdict = {'INTWFC1':mask1,\
            'INTWFC2':mask2,\
            'INTWFC3':mask3,\
            'INTWFC4':mask4}

##########################################################
### SET UP DIRECTORY AND HEADER FIELDS TO BRING IN
##########################################################
imagedir = '/home/rfinn/data/reduced/scratch-int-feb2019/attempt2/pointing149/'
keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']

##########################################################
### GET RID OF SHORT EXPOSURES
##########################################################
ic = ccdp.ImageFileCollection(imagedir, keywords=keys, glob_include='WFC.r.*4PA.fits',glob_exclude='*coadd*.fits')

# get rid of 15 sec r-band exposure, b/c can't figure out how to do the weighting!!!
if not os.path.exists('shortexposure'):
    os.mkdir('shortexposure')
expflag = ic.summary['exptime'] < 20
if sum(expflag) > 0:
    quick_files = ic.files[expflag]
    for f in quick_files:
        os.rename(f,'shortexposure/'+f)

##########################################################
### CREATE NEW IMAGE COLLECTION WITHOUT SHORT EXPOSURE
##########################################################
ic = ccdp.ImageFileCollection(imagedir, keywords=keys, glob_include='WFC.*PA.fits',glob_exclude='*coadd*.fits')

expt = ic.summary['exptime']
myweights = np.array(1/expt,'f')



##########################################################
### ADD MASKS TO IMAGES
##########################################################
addmask=False
if addmask:
    # add masks and save images as CCDData objects
    ccddata = []
    for i,h in enumerate(myhdus):
        ccddata.append(CCDData(h.data,mask=maskdict[h.header['INSTRMNT']].data,unit='adu',wcs=WCS(h.header)))


##########################################################
### SUBTRACT MEDIAN FROM IMAGE
##########################################################

#print('image statistics before subtracting median')
#for h in myhdus:
#    print(np.median(h.data))
# subtract median for all images
print('subtracting median from images')
for hdu in ic.hdus:
    # background subtraction
    hdu.data = imutils.subtract_median_sky(hdu.data)

#print('image statistics after subtracting median')
#for h in myhdus:
#    print(np.median(h.data))

##########################################################
### FIND BEST WCS FOR ALL IMAGES
##########################################################

# find wcs for all images
print('finding optimal output wcs')
wcs_out, shape_out = find_optimal_celestial_wcs(ic.hdus)


##########################################################
### REPROJECTING IMAGES
##########################################################
#print('Reprojecting Images, please be patient...')
#array, footprint = reproject_and_coadd(ccddata,wcs_out,shape_out=shape_out,\
#                                       reproject_function=reproject_interp,combine_function='mean') # change to exact once this works
#outim = CCDData(array,wcs=wcs_out,unit='adu')
#outweight = CCDData(footprint,wcs=wcs_out,unit='adu')
#outim.write('testcoadd.fits',overwrite=True)
#outim.write('testcoadd.weight.fits',overwrite=True)

print('Starting loop to reproject images...')
filters = set(ic.summary['filter'])

if not os.path.exists('reprojected'):
    os.mkdir("reprojected")
    
for f in filters:
    reprojected  = []
    allwcs = []
    allweight = np.zeros(shape_out)
    allim = np.zeros(shape_out)
    for i,imagename in enumerate(ic.files):
        if ic.summary['filter'][i] == f:
            print('processing ',imagename)

            w = WCS(myhdus[i].header)
            
            # change to reproject_exact once code works
            reprodata, footprint = reproject_interp(ccddata[i], wcs_out,shape_out=shape_out)
            
            reprodata = wcs_project(cdata,wcs_out)
            newimg = CCDData(reprodata, wcs=wcs_out, unit='adu',mask=exposuremask)                
            reprojected.append(newimg)
            allweight += footprint


def combine_images(dir):
    workingdir = os.getcwd()
    os.chdir("reprojected")
    ic = ccdp.ImageFileCollection(imagedir, keywords=keys, glob_include=args.filestring+'*.fits',glob_exclude='*coadd*.fits')
    filters = set(ic.summary['filter'])
    for f in filters:
        images = ic.files_filtered(filter=f)
        combined_image = ccdp.combine(images,method='average',\
                                  sigma_clip=True,sigma_clip_low_thresh=2,\
                                  sigma_clip_high_thresh=5,sigma_clip_func=np.ma.median,\
                                  sigma_clip_dev_func=mad_std,minmax_clip=True,nlow=1,nhigh=0)
        combined_image.write(f+'coadd.fits',overwrite=True)
        weightim = CCDData(allweight,wcs=wcs_out,unit='adu')
        weightim.write(f+'coadd.weight.fits',overwrite=True)    

