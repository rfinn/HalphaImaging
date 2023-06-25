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
    hfile = k.replace('.fits','.head')
    if not os.path.exists(mfile) and os.path.exists(hfile):
        print(k)
        os.rename(k,outdir+'/'+k)
    elif not os.path.exists(mfile) and not os.path.exists(hfile):
        print(f"{k} missing {hfile}")
    
