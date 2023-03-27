#!/usr/bin/env python

"""

GOAL:
create a supersky flat by combining science frames taken through a given filter

PROCEDURE:


USAGE:
%run ~/github/HalphaImaging/python3/BOK_03_sky_flat.py Ha4nm

where the string after the program name is the filter, like r or Ha4nm

NOTE:
This is REALLY, REALLY slow.  
It takes about 2-3 minutes per amp, so 45 min for all 16 amps.  Yuck!

Should run this on an HPC instead of my laptop!
"""


import glob
import os
import numpy as np

from astropy.io import fits
from astropy.stats import mad_std
from astropy import units as u
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS

import timeit

import ccdproc as ccdp
from photutils import make_source_mask
import multiprocessing as mp

homedir = os.getenv("HOME")
import sys
sys.path.append(os.path.join(homedir,'/github/HalpaImaging/python3/'))
import imutils



usemp = True

filter = sys.argv[1]
amp = int(sys.argv[2])
medsky_results = []
def collect_results_medsky(result):

    global results
    medsky_results.append(result)


def inv_median(a):
    return 1/np.ma.median(a)

def get_masked_image(filename):
    hdu = fits.open(filename)
    mask = make_source_mask(hdu[0].data,nsigma=3,npixels=5,dilate_size=20)
    masked_data = np.ma.array(hdu[0].data,mask=mask)
    #clipped_array = sigma_clip(masked_data,cenfunc=np.ma.mean)

    mean,median,std = sigma_clipped_stats(masked_data,sigma=3.0,cenfunc=np.ma.mean)
    
    # kluging this because I can't figure out how to return image to mp.pool
    this_median.append(median)
    masked_images.append(CCDData(hdu[0].data,unit=u.electron,mask=mask))

    hdu.close()

def get_bright_star_flag(filelist):
    staronimage = np.zeros(len(filelist),'bool')
    brightstar = Table.read('/home/rfinn/research/legacy/gaia-mask-dr9-bright-Virgo-g-lt-8.fits')
    starcoord = SkyCoord(brightstar['ra'],brightstar['dec'],frame='icrs',unit='deg')
    for i,im in enumerate(filelist):
        h = fits.getheader(im)
        w = WCS(h)
        x,y = w.world_to_pixel(starcoord)
        flag = (x > 0) & (x < h['NAXIS1']) & (y>0) & (y < h['NAXIS2'])        
        if np.sum(flag) > 0:
            staronimage[i] = True
    return staronimage

if not os.path.exists('SKYFLAT'):
    os.mkdir('SKYFLAT')

# get list of flattene images from each subdir

subfolders = [f.path for f in os.scandir(os.getcwd()) if f.is_dir()]

filelist = []
for f in subfolders:
    if f'target-{filter}' in f:
        filelist += glob.glob(f+'/*.fits')

#print(filelist)
amps = np.arange(1,17)
all_skyflats = []
all_med = []
skyfiles = []
amps = [amp]
for a in amps:

    a1 = []
    # cheating to help 
    this_median = []
    masked_images = []
    for f in filelist:
        if f'_{a}PA.fits' in f:
            a1.append(f)
    t_0 = timeit.default_timer()
    print("#################################")
    print(f"working on amp {a}")
    print("\tlooking for bright stars...")
    # look for bright stars and remove frames with bright stars
    bright_star_flag = get_bright_star_flag(a1)
    # print update
    print(f"\tremoving {np.sum(bright_star_flag)} images due to bright stars")
    print()
    # remove images with bright stars
    a1 = list(np.array(a1)[~bright_star_flag])
    t_1 = timeit.default_timer()
    print(f"\telapsed time: {round((t_1-t_0),3)} sec")
    
    print(f"\tprocessing {len(a1)} images for amp {a} skyflat")    
    # combine files
    # scale by median
    # method = median
    print("\tgetting object masks")
    if usemp:
        medsky_pool = mp.Pool(mp.cpu_count())
        myresults = [medsky_pool.apply_async(get_masked_image(im),callback=collect_results_medsky) for im in a1]
    
        medsky_pool.close()
        medsky_pool.join()
    else:
        for im in a1:
            get_masked_image(im)

    t_2 = timeit.default_timer()
    print(f"\telapsed time: {round((t_2-t_1),3)} sec")

    print("\tcombining masked images")    
    skyflat = ccdp.combine(masked_images,
                           method='median', scale=1/np.array(this_median),
                           sigma_clip=True, sigma_clip_low_thresh=2.5, sigma_clip_high_thresh=2.5,
                           sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,unit=u.electron)

    # testing clip_extrema instead of sigma_clip
    # try throwing out the top 1/3 of the images so that objects/stars are gone
    #skyflat = ccdp.combine(masked_images,
    #                       method='median', #scale=1/np.array(this_median),
    #                       clip_extreme=True, nlow=int(len(a1)/3), nhigh=int(len(a1)/3),unit=u.electron)
    t_3 = timeit.default_timer()
    print(f"\telapsed time: {round((t_3-t_2),3)} sec")

    all_skyflats.append(skyflat)
    all_med.append(np.median(np.array(this_median)))


    skyflat.write(f'SKYFLAT/skyflat{a:02d}.fits',overwrite=True)
    skyfiles.append(f'SKYFLAT/skyflat{a:02d}.fits')

    medvalues = open(f'SKYFLAT/medvalue-{a}.dat','w')
    for m in all_med:
        medvalues.write(f"{m:.3f}\n")
    medvalues.close()
    t_4 = timeit.default_timer()
    print(f"\ttime for this amp: {round((t_3-t_0),3)} sec")
    

'''
all_med = np.array(all_med)
mean_sky = np.mean(all_med)

norm_factor = mean_sky/all_med


for i,f in enumerate(all_skyflats):
    hdu = fits.open(skyfiles[i])
    hdu[0].data = hdu[0].data*norm_factor[i]
    hdu[0].header.set('skyflat', True)    
    hdu[0].header.set('scale', norm_factor[i])
    hdu.writeto(f'SKYFLAT/skyflat{i+1:02d}.fits',overwrite=True)

'''
