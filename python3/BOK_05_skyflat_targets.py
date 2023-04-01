#!/usr/bin/env python
"""
flatten science images with sky flat

"""
import os
from astropy.io import fits
import numpy as np
import sys
import glob

filter = sys.argv[1]
# get list of flattened images from each subdir
subfolders = [f.path for f in os.scandir(os.getcwd()) if f.is_dir()]

filelist = []
for f in subfolders:
    if f'target-{filter}' in f:
        filelist += glob.glob(f+'/*.fits')

#print(filelist)
amps = np.arange(1,17)
for a in amps:

    a1 = []
    # cheating to help 
    this_median = []
    masked_images = []
    for f in filelist:
        if f'_{a}PA.fits' in f:
            a1.append(f)

    print("#################################")
    print(f"working on amp {a}")
    
    print(f"\tprocessing {len(a1)} images for amp {a} skyflat")    
    # combine files
    # scale by median
    # method = median

    skyflat=f'SKYFLAT-{filter}/nskyflat{a:02d}.fits'
    flat = fits.open(skyflat)
    fdata = flat[0].data
    flat.close()
    for f in a1:
        hdu = fits.open(f)
        hdu[0].data = hdu[0].data/fdata
        hdu[0].header.set("skyflat",1,"skyflat applied")
        #basename = os.path.basename(f)
        #fout = f.replace(basename,'f'+basename)
        hdu.writeto(f,overwrite=True,output_verify='ignore')
        hdu.close()
