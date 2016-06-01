#!/usr/bin/env python

'''
PURPOSE:

The goal of the program is to create a mask for a galaxy image so that 

USAGE:

PROCEDURE:

'''

import os
import sys
from astropy.io import fits
import numpy as np
from pyraf import iraf
import argparse


defaultcat='default.sex.sdss'

mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    sedir='/Users/rfinn/research/LocalClusters/sextractor/'
elif mypath.find('home') > -1:
    print "Running on coma"
    sedir='/home/share/research/LocalClusters/sextractor/'

parser=argparse.ArgumentParser()
parser.add_argument("agcnumber",help='agcnumber')
#parser.add_argument("-t",'--threshold',help="sextractor DEBLEND_MINCONT: 0=lots of deblending; 1=none (default = .005)",action="store")
#parser.add_argument("-c",'--catalog', help="optional input, sextractor default.sex file",action="store")
#parser.add_argument("-d",'--display', help="display result in ds9",action="store_true")
args = parser.parse_args()


#
# assume galaxy is at the center of the image
fdulist = fits.open(image)
t=fdulist[0].data
n2,n1=t.shape
xc=n1/2.
yc=n2/2.
fdulist.close()

# run sextractor to generate a list of objects in the image
# generate 'segmentation image'
#try:
#    os.system('ln -s /home/share/research/LocalClusters/sextractor/* .')
#except:
#    print 'could not link to /home/share/research/LocalClusters/sextractor'
#    print defaultcat,' should be in the current directory'
sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
for file in sextractor_files:
    os.system('ln -s '+sedir+file+' .')

if args.threshold:
    os.system('sex %s -c %s -CATALOG_NAME test.cat -CATALOG_TYPE ASCII -DEBLEND_MINCONT %f'%(image,defaultcat,deblend))
else:
    os.system('sex %s -c %s -CATALOG_NAME test.cat -CATALOG_TYPE ASCII'%(image,defaultcat))

#   parse sextractor output to get x,y coords of objects        

sexout=np.loadtxt('test.cat',usecols=[0,1,2],comments='#')
sexnumber=sexout[:,0]
xsex=sexout[:,1]
ysex=sexout[:,2]
dist=np.sqrt((yc-ysex)**2+(xc-xsex)**2)

#   find object ID
objIndex=np.where(dist == min(dist))
objNumber=sexnumber[objIndex]
objNumber=objNumber[0] # not sure why above line returns a list
mask_image=image.split('.fits')[0]+'mask.fits'


fdulist = fits.open('segmentation.fits',mode='update')
t=fdulist[0].data
# replace object ID values with zero
f = (t == objNumber)
replace_values = np.zeros(t.shape)
t=t*(~f) + replace_values*f
# set all bad values to 1
t[t> .1]=t[t > .1]/t[t > .1]
fdulist.flush()
# write out mask image
if os.path.isfile(mask_image):
    os.remove(mask_image)

hdu=fits.PrimaryHDU(t)
hdu.writeto(mask_image)

fdulist.close()

#   use iraf imexam to replace object ID values with zero
#iraf.imreplace(working_dir+mask_image,value=0,lower=objNumber-.5,upper=objNumber+.5)
# convert segmentation image to object mask by replacing the object ID of target galaxy with zeros

# clean up
sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
for file in sextractor_files:
    os.system('unlink '+file)

if args.display:
    try:
        d.set('frame delete all')
    except:
        d=ds9.ds9()
        d.set('frame delete all')
    d.set('file new '+image)
    d.set('zscale')
    d.set('file new '+mask_image)
    d.set('zscale')
