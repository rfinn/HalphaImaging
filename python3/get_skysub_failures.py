#!/usr/bin/env python

"""
GOAL: print out files that have sky subtraction but fixamps failed.

"""
import glob
import os

kfiles = glob.glob('ksb*ooi*v1.fits')
kfiles.sort()
outdir = 'skysub-failures'
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
for k in kfiles:
    mfile = 'm'+k
    if not os.path.exists(mfile):
        print(k)
        os.rename(k,outdir+'/'+k)
    
