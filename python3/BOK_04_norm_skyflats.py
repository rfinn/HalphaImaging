#!/usr/bin/env python

"""
GOAL: 

normalize the skyflats for the entire mosaic

USAGE:

python ~/github/HalpaImaging/python3/BOK_04_norm_skyflats.py r

"""
import numpy as np
from astropy.io import fits
import sys

filter = sys.argv[1]
# get median values
med = []
for i in range(16):
    d = np.loadtxt(f'SKYFLAT-{filter}/medvalue-{i+1}.dat')
    #print(d)
    med.append(d)
med = np.array(med)
# zeros are place holders for the count levels in
# amps with bright stars
maskdata = med < 0.1

# create a masked array of median sky levels
mmed = np.ma.array(med,mask=maskdata)

# normalize each image column by the mean counts in each image
norm_mmed = mmed/np.ma.mean(mmed,axis=0)

# now take the median in each amp
ampvals = np.ma.median(norm_mmed,axis=1)


## track images where no chip is rejected due to bright stars
#imflag = np.sum(maskdata,axis=0) == 0
#
## calculate statistics of amps for the images where all amps are good
#ampvals = np.median(med[:,imflag],axis=1) 

# norm flats
ave = np.mean(ampvals)
print("ave = ",ave)

# number to scale each image by
# images will be divided by this flat, so need to do med/ave
# this way, images with higher sky values will get divided by bigger number
scale = ampvals/ave


for i in range(16):
    imname = f'SKYFLAT-{filter}/skyflat{i+1:02d}.fits'
    print(f"{imname}, {ampvals[i]:.2f}, {scale[i]:.5f}")
    hdu = fits.open(imname)
    hdu[0].data = hdu[0].data * scale[i]

    # write out scaled data
    hdu.writeto(f'SKYFLAT-{filter}/nskyflat{i+1:02d}.fits',overwrite=True)

