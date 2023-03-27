#!/usr/bin/env python
"""
GOAL: fix the wcs headers in the 90prime data

PROCEDURE:

* According to an email from Taran LeRoy Esplin on 3/15/2021
######################################################################### 
Starting with the RA and DEC written by the the telescope to the headers: 
amplifiers 1-16 (extensions 2-17 of a raw file) have the following 
approximate offsets from the (central) RA and DEC (in decimal degrees):



Order of headers:
0 = primary header, this has telescope RA and DEC

remaining headers are in a weird order:

IM4
IM3
IM2
IM1
IM8
IM7
IM6
IM5
IM9
IM10
IM11
IM12
IM13
IM14
IM15
IM16

######################################################################### 

"""
#from ccdproc import ImageFileCollection
import os
import glob
import numpy as np
from astropy.coordinates import Angle
from astropy.io import fits
from astropy import units as u
import argparse

# use argparse to get the file prefix

parser = argparse.ArgumentParser()
parser.add_argument('--filestring', default='ut', help='Prefix of image files.  The default is ut.  All images with ut*.fits will be operated on.')
parser.add_argument('--nooverwrite', action='store_true', default=False, help='Default behavior is to overwrite existing files.  Set this to save output as w+input_filename.')
args = parser.parse_args()

RA_offsets = np.array([0.1536,0.1543,0.4118,0.4127,\
                       0.1500,0.1513,0.4066,0.4089,\
                       -0.1508,-0.1512,-0.4093,-0.4096,\
                       -0.1537,-0.1541,-0.4114,-0.4107])

DEC_offsets_extensions = np.array([0.1346,0.3886,0.133,0.38606,\
                                   -0.3912,-0.1377,-0.3922,-0.1387,\
                                   0.392,0.138,0.3911,0.1376,\
                                   -0.134,-0.388,-0.1335,-0.387])
CD1_1 = 0

CD1_2 = np.array([+1.257e-4,+1.257e-4,-1.257e-4,-1.257e-4,\
                  +1.257e-4,+1.257e-4,-1.257e-4,-1.257e-4,\
                  -1.257e-4,-1.257e-4,+1.257e-4,+1.257e-4,\
                  -1.257e-4,-1.257e-4,+1.257e-4,+1.257e-4])

CD2_1 = np.array([+1.257e-4,-1.257e-4,+1.257e-4,-1.257e-4,\
                  +1.257e-4,-1.257e-4,+1.257e-4,-1.257e-4,\
                  -1.257e-4,+1.257e-4,-1.257e-4,+1.257e-4,\
                  -1.257e-4,+1.257e-4,-1.257e-4,+1.257e-4])

CD2_2 = 0

other_fields = ['SIMPLE','BITPIX',\
                'NAXIS','NAXIS1','NAXIS2',\
                'CRPIX1','CRPIX2',\
                'CTYPE1','CTYPE2','CUNIT1','CUNIT2']

field_values = ['T',-32,\
                2,2016,2048,\
                1008,1024,\
                'RA---TAN','DEC--TAN','deg','deg']


# From Zou+2017 https://iopscience.iop.org/article/10.3847/1538-3881/aa72d9
#Table 2 
#Gain and Readout Noise for Each CCD
#CCD	CCD #1				CCD #2				CCD #3				CCD #4				
#HDU	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	
#Gain	1.495	1.620	1.524	1.519	1.546	1.656	1.608	1.590	1.399	1.503	1.425	1.553	1.505	1.535	1.521	1.584	
#Noise	7.353	7.288	7.313	6.597	9.299	9.509	9.302	8.006	10.122	8.897	9.474	9.019	9.527	7.940	6.903	7.716	
gain = np.array([1.495,	1.620,	1.524,	1.519,\
                 1.546,	1.656,	1.608,	1.590,\
                 1.399,	1.503,	1.425,	1.553,\
                 1.505,	1.535,	1.521, 1.584])	

# get list of images
filenames = glob.glob(args.filestring+'*.fits')

# loop over images
for f in filenames:

    hdu = fits.open(f)

    # get RA DEC from primary hdu
    RA = hdu[0].header['RA']
    DEC = hdu[0].header['DEC']
    #print(f,RA,DEC)
    # convert hexigesimal to degrees
    ra_deg = Angle(RA,unit=u.hour).deg
    dec_deg = Angle(DEC,unit=u.deg).deg    
    #print(ra_deg,dec_deg)

    # correct for cosine of declination
    RA_offsets_extensions = RA_offsets/np.cos(dec_deg*np.pi/180)

    # loop through extensions
    for i in range(1,len(hdu)):
        # update wcs
        # update CRVAL1
        hdu[i].header['CRVAL1'] = ra_deg + RA_offsets_extensions[i-1]
        # update CRVAL2
        hdu[i].header['CRVAL2'] = dec_deg + DEC_offsets_extensions[i-1]
        # update CD1_1
        hdu[i].header['CD1_1'] = CD1_1      
        # update CD1_2
        hdu[i].header['CD1_2'] = CD1_2[i-1]
        # update CD2_1
        hdu[i].header['CD2_1'] = CD2_1[i-1]
        # update CD2_2
        hdu[i].header['CD2_2'] = CD2_2

        # what about the reference pixels???
        # leaving if for now and will see what it looks like
        for k,field in enumerate(other_fields):
            #print(field,field_values[k])
            hdu[i].header[field] = field_values[k]
        # update the GAIN in each image
        hdu[i].header.set('GAIN',gain[i-1],'GAIN in e/ADU Zou 2017')

    # add comment to primary header to note that wcs is corrected
    hdu[0].header.set('fixWCS',1,'rough fix to 90prime WCS')

    # write image
    print(f"updating WCS header info for {f} ({RA},{DEC})")
    if args.nooverwrite:
        hdu.writeto('w'+f,overwrite=True)
    else:
        hdu.writeto(f,overwrite=True)
