#!/usr/bin/env python

"""
GOAL: 

to adjust the bias levels in flats and targets based on how the bias drifts with time

"""


'''
GOAL:


PROCEDURE:

- run theli through splitting all images

- create the average bias for frames taken at the beginning of the night (before UT day = 0.3) 
- fit the line to each amp using the following, where the line is anchored by the bias counts 
  at the start of a given day.  For each night, we derive an equation for the linear drift 
  in bias versus time, using the slope that was fitted with linear regression to all the 
  bias data from the April 2022 run.

  slope(amp) = (yf - yi)/(xf-xi)

  bias_obs = slope(amp)*(UT_obs - UT_bias) + med_bias

- for each science flat

  bias_subtracted = science - bias_obs


USAGE:

- run this from the root directory for each night, 

- after:
  - fix wcs
  - sort files: the data have been sorted into BIAS, FLAT, TARGET-r, TARGET-Ha4nm
  - BIAS frames have been combined by theli




REFERENCES:

https://ccdproc.readthedocs.io/en/latest/reduction_toolbox.html

'''
#import ccdproc
from ccdproc import ImageFileCollection
import os
from astropy.io import fits
from astropy.time import Time
import numpy as np
import glob

# these are the slopes that we fit to the April 2022 bias data
# fitting counts as a function of UT hr/24.  X range of data is from 0 - hour < 0.6 days
# Ryan Keenan took data at the beginning and end of the night on several nights
biasslope = {1:-54.464, 2:-56.099, 3:-24.081, 4:-31.598, 5:9.965, 6:-40.599, 7:16.967, 8:-14.124, 9:-10.698, 10:-8.141, 11:-26.000, 12:-79.479, 13:-22.814, 14:-32.337, 15:14.889, 16:-46.241, }

# using data from the beginning of the run only

biasslope = {1:-40.699, 2:-41.567, 3:-19.300, 4:-23.278, 5:3.528, 6:-31.717, 7:11.185, 8:-12.076, 9:-8.032, 10:-5.718, 11:-17.136, 12:-57.562, 13:-18.726, 14:-26.345, 15:11.223, 16:-34.378, }
# number of amps in 90prime
namps=16

# save top directory name
basedir = os.getcwd()


def get_bias_UT():
    # pick an image
    # check if theli preserves the UT time of observation in the bias - it does not!!!
    # need to get the UT for the bias frame
    os.chdir('BIAS')
    bias_files = glob.glob('*1P.fits')
    dateobs = []
    for b in bias_files:
        header = fits.getheader(b)
        dateobs.append(Time(header['DATE-OBS']).mjd)
    dateobs = np.array(dateobs)
    # get the mean UT of the bias frames
    bias_mjd = np.mean(dateobs)
    os.chdir(basedir)
    # return modified julian date
    return bias_mjd
                    
def subtract_bias(ic,amp=1,bias_mjd=0.):

    # open the bias frame
    biasfile = basedir+f"/BIAS/BIAS_{amp}.fits"
    bhdu = fits.open(biasfile)
    bias = bhdu[0].data
    bias_fracday = bias_mjd%1
    bhdu.close()

    # get median of bias from, so set intercept of linear scaling relation with time
    # BTW, this all made sense to us at some point...

    medbias = np.median(bias)

    # loop through images and subtract bias
    for f in ic.values('file'):
        f = os.path.basename(f)
        hdu = fits.open(f)
        tobs = Time(hdu[0].header['DATE-OBS'])
        fracday = tobs.mjd%1
        
        # get scale factor for bias
        bias_offset = biasslope[amp]*(fracday - bias_fracday) 

        # subtract the bias

        hdu[0].data += (bias + bias_offset)

        # add card to head to record that bias has been subtracted
        hdu[0].header.set('BIASCOR','vfproc','no overscan im')
        hdu[0].header.set('BIASLOPE',f"{biasslope[amp]:.2f}",'slope of bias drift vs UThr/24')
        hdu[0].header.set('BIASOFF',f"{bias_offset:.2f}",'bias offset for this image/amp')
        
        # write out bias-subtracted file
        # writing to a new file while I am testing
        hdu.writeto(f,overwrite=True)            
        hdu.close()


bias_mjd = get_bias_UT()
# subtract bias from flatfield images

filters = ['r','Ha4nm']
subdir_names = ['FLAT-','target-']
# FOR TESTING, JUST USING R AND TARGET
#filters = ['r']
#subdir_names = ['target-']
for filt in filters:
    for sd in subdir_names:
        subdir_name = sd+filt
        print(f"working on directory {subdir_name}\n")
        os.chdir(subdir_name)
        for i in range(namps):
            filestring = f"*_{i+1}P.fits"
            ic = ImageFileCollection(os.getcwd(),keywords='*',glob_include=filestring)    
    
            subtract_bias(ic,amp=i+1,bias_mjd=bias_mjd)
        os.chdir(basedir)
    
