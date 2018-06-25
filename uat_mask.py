#!/usr/bin/env python

'''
PURPOSE:

The goal of the program is to create a mask for a galaxy image to mask
out other objects within the cutout area.

USAGE:

from within ipython type:

   %run ~/github/HalphaImaging/uat_mask.py --image 'A1367-140231-R.fits'

you just need to run this on R-band images.

Interacting with the display is finicky.  This works fine when running
within ipython - not so much when running from the command line.  When running
from the command line, I am not able to interact with the figure.  This may
have something to do with setting block=False in show().


PROCEDURE:


REQUIRED MODULES:
   os
   astropy
   numpy
   argsparse
   matplotlib
   scipy


'''

import os
import sys
from astropy.io import fits
import numpy as np
import argparse
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

defaultcat='default.sex.HDI.mask'


parser = argparse.ArgumentParser(description ='Create a mask for extraneous objects in field')
parser.add_argument('--image', dest = 'image', default = None, help = 'image to mask')
parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of default config files')
parser.add_argument('--threshold', dest = 'threshold', default = .005, help = "sextractor DEBLEND_MINCONT: 0=lots of deblending; 1=none (default = .005)",action="store")
parser.add_argument('--snr', dest = 'snr', default = .005, help = "sextractor DEBLEND_MINCONT: 0=lots of deblending; 1=none (default = .005)",action="store")
parser.add_argument('--cmap', dest = 'cmap', default = 'gist_heat_r', help = "color map to use when displaying image mask.  default is gist_heat_r.") 
args = parser.parse_args()

sextractor_files=['default.sex.HDI.mask','default.param','default.conv','default.nnw']
for file in sextractor_files:
    os.system('ln -s '+args.d+'/'+file+' .')

def read_se_cat():
    sexout=fits.getdata('test.cat')
    xsex=sexout['XWIN_IMAGE']
    ysex=sexout['YWIN_IMAGE']
    dist=np.sqrt((yc-ysex)**2+(xc-xsex)**2)
    #   find object ID
    objIndex=np.where(dist == min(dist))
    objNumber=sexout['NUMBER'][objIndex]
    return objNumber[0] # not sure why above line returns a list

def runse():
    print 'using a deblending threshold = ',args.threshold
    
    os.system('sex %s -c %s -CATALOG_NAME test.cat -CATALOG_TYPE FITS_1.0 -DEBLEND_MINCONT %f'%(args.image,defaultcat,args.threshold))
    maskdat = fits.getdata('segmentation.fits')
    center_object = read_se_cat()
    maskdat[maskdat == center_object] = 0
    return maskdat

adjust_mask = True
figure_size = (6,3)

# create name for output mask file
t = args.image.split('-')
mask_image=t[0]+'-'+t[1]+'-'+t[2]+'-mask.fits'
print 'saving image as: ',mask_image

# read in image and define center coords
image, imheader = fits.getdata(args.image,header = True)
yc,xc = image.shape
xc = xc/2.
yc = yc/2.

v1,v2=scoreatpercentile(image,[5.,99.5])

# run sextractor on input image
# return segmentation image with central object removed
maskdat = runse()

while adjust_mask:
    plt.close('all')
    plt.figure(1,figsize=figure_size)
    plt.clf()
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.subplot(1,2,1)
    plt.imshow(image,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
    plt.title('image')
    plt.subplot(1,2,2)
    #plt.imshow(maskdat,cmap='gray_r',origin='lower')
    plt.imshow(maskdat,cmap=args.cmap,origin='lower')
    plt.title('mask')
    plt.gca().set_yticks(())
    #plt.draw()
    plt.show(block=False)
    
    
    t=raw_input('enter:\n   pixel value to remove object in mask;\n   t to adjust SE threshold (0=no deblend, 1=lots); \n   w to write output and quit; \n   q to quit without saving\n')
    try:
        objID = int(t)
        maskdat[maskdat == objID] = 0.
    except ValueError:
        adjust_scale = False
        if t.find('t') > -1:
            t = raw_input('enter new threshold')
            args.threshold = float(t)
            runse()
        if t.find('q') > -1:
            sys.exit()
        if t.find('w') > -1:
            newfile = fits.PrimaryHDU()
            newfile.data = maskdat
            newfile.header = imheader
            fits.writeto(mask_image, newfile.data, header = newfile.header, clobber=True)
            adjust_mask = False

# clean up
#sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
for file in sextractor_files:
    os.system('unlink '+file)


    
