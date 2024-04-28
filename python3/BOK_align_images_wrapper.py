#!/usr/bin/env python

"""
GOAL: 

wrapper to INT_align_images.py where you can pass in just one argument so that it's easy to run with gnu parallel

USAGE:

python ~/github/HalphaImaging/python3/BOK_align_images.py  VF-172.650+35.671-BOK-20210417-VFID2484-Ha4.fits


NOTES:

On draco, 

cd /data-pool/Halpha/coadds/BOK-coadds-v0
ls *Ha4.fits > Ha4_filelist
wc -l Ha4_filelist

should have 66 files


Then to test with one galaxy

python ~/github/HalphaImaging/python3/BOK_align_images_wrapper.py  VF-172.650+35.671-BOK-20210417-VFID2484-Ha4.fits


and then run using parallel:

parallel --eta python ~/github/HalphaImaging/python3/BOK_align_images_wrapper.py  :::: Ha4_filelist

"""
import os
import sys
import glob

himage = sys.argv[1]
rimage = himage.replace('Ha4.fits','r.fits')
weightimage = rimage.replace('r.fits','r.weight.fits')

rshifted = rimage.replace('r.fits','r-shifted.fits')
if not os.path.exists(rimage):
    print("warning: default r image name is not working for ",himage, rimage)
    t = himage.split('-')
    ######################################################
    # check to see if r-band was taken on a different date
    ######################################################

    # get date from filestring
    obsdate = t[-3]

    # replace date with *
    search_name = rimage.replace(obsdate,'*')

    # look for files using glob
    filelist = glob.glob(search_name)

    # if a match was found, redefine names
    if len(filelist) > 0:
        rimage = filelist[0]
        weightimage = rimage.replace('r.fits','r.weight.fits')
        rshifted = rimage.replace('r.fits','r-shifted.fits')
if os.path.exists(rshifted):
    print("shifted image found - skipping ",himage)
    sys.exit()
os.system(f"python ~/github/HalphaImaging/python3/INT_align_images.py --image1 {himage} --image2 {rimage} --weight2 {weightimage}")
