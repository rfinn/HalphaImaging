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

python ~/github/HalphaImaging/python3/BOK_align_images.py  VF-172.650+35.671-BOK-20210417-VFID2484-Ha4.fits


and then run using parallel:

parallel --eta python ~/github/HalphaImaging/python3/BOK_align_images.py  :::: Ha4_filelist

"""
import os
import sys

himage = sys.argv[1]
rimage = haimage.replace('Ha4.fits','r.fits')
weightimage = haimage.replace('r.fits','r.weight.fits')

os.system(f"python ~/github/HalphaImaging/python3/INT_align_images.py --image1 {himage} --image2 {rimage} --weight2 {weightimage}")
