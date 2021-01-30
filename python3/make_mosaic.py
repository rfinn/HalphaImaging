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

from photutils import make_source_mask
import os
homedir = os.getenv("HOME")
import sys
print(homedir)

sys.path.append('/home/rfinn/github/HalphaImaging/python3/')
import imutils

#from reproject import reproject_interp
#from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs
import ccdproc as ccdp
from ccdproc import wcs_project
from ccdproc import Combiner
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import CCDData

# read in WFC ccd masks
mask1 = fits.getdata(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits')
mask2 = fits.getdata(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits')
mask3 = fits.getdata(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits')
mask4 = fits.getdata(homedir+'/data/WFCINTmasks/WFCINTccdmask1.fits')


# dictionary to translate between instrument in header and mask file
maskdict = {'INTWFC1':mask1,\
            'INTWFC2':mask2,\
            'INTWFC3':mask3,\
            'INTWFC4':mask4}

imagedir = '/home/rfinn/data/reduced/scratch-int-feb2019/attempt2/pointing149/'
keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']
ic = ccdp.ImageFileCollection(imagedir, keywords=keys, glob_include='WFC*.fits',glob_exclude='*coadd*.fits')

# find wcs for all images
wcs_out, shape_out = find_optimal_celestial_wcs(ic.files)


# for each image

filters = set(ic.summary['filter'])

for f in filters:
    allimages  = []
    allwcs = []

    for imagename in ic.files_filtered(filter=f):
        data = fits.getdata(imagename)
        # check to see if header file exists
        filebasename = imagename.split('.fits')[0]
        altheader = filebasename+'.head'
        if os.path.exists(altheader):
            wcsinfo = imutils.convert_headfile_header(altheader)
            header = fits.getheader(imagename)
            for key in wcsinfo:
                if key.startswith('PV'):
                    continue
                header.set(key,value=wcsinfo[key])
        else:
            header = fits.getheader(imagename)
        w = WCS(header)
        # background subtraction
        data = imutils.subtract_median_sky(data)
        # convert data to CCDdata and add mask
        cdata = CCDData(data,mask=maskdict[header['INSTRMNT']],unit='adu',wcs=w)
        #cdata.mask = maskdict[header['INSTRMNT']]
        reprodata = wcs_project(cdata,wcs_out)
        allimages.append(reprodata)
        # reproject individual images onto combined wcs
    combiner = Combiner(allimages)
    #combiner.sigma_clipping(low_thresh=2, high_thresh=5, func=np.ma.median)
    stacked_image = combiner.average_combine()
    stacked_image.write(f+'coadd.fits')
# mosaic




# coaddition by filter
# apply badpixel masks
