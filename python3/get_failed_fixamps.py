#!/usr/bin/env python

"""
GOAL: print out files that have sky subtraction but fixamps failed.

"""
import glob
import os

mfiles = glob.glob('mksb*ooi*v1.fits')
mfiles.sort()
for m in mfiles:
    zfile = 'z'+m
    if not os.path.exists(zfile):
        print(m)
