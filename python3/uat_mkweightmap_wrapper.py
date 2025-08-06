#!/usr/bin/env python

"""
GOAL:
- run uat_MOSmkweightmap.py on images with the same prefix

USAGE:

python ~/github/HalphaImaging/python3/uat_MOSmkweightmap_wrapper.py --filestring iwf

this will find all files that match iwf*.fits, and for each file, it will call uat_mkweightmap.py
"""


import argparse
from astropy.io import fits
import glob
import os

homedir = os.getenv("HOME")

parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--filestring', dest = 'filestring', default = 'h', help = 'string to use to get input files (default = "h" which grabs all files "h*.fits")')

args = parser.parse_args()


filelist = glob.glob(args.filestring+"*.fits")

for f in filelist:
    os.system(f"{homedir}/github/HalphaImaging/python3/mkweightmap.py --filename {f} --badval 0")
