#!/usr/bin/env python

"""
GOAL: print out files that have sky subtraction but fixamps failed.

"""
import glob
import os

mfiles = glob.glob('mksb*ooi*v1.fits')
mfiles.sort()
outdir = 'fixamp-failures'
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
for m in mfiles:
    zfile = 'z'+m
    hfile = m.replace('.fits','.head')
    if not os.path.exists(zfile) and os.path.exists(hfile):
        print(m)
        os.rename(m,outdir+'/'+m)
    elif not os.path.exists(zfile) and not os.path.exists(hfile):
        print(f"{m} missing {hfile}")
    
