#!/usr/bin/env python

"""
GOAL: print out files that have sky subtraction but fixamps failed.

"""
import glob
import os

mfiles = glob.glob('ksb*ooi*v1.fits')
mfiles.sort()
outdir = 'skysub-failures'
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
for k in kfiles:
    zfile = 'm'+k
    hfile = k.replace('.fits','.head')
    if not os.path.exists(mfile) and os.path.exists(hfile):
        print(m)
        os.rename(m,outdir+'/'+m)
    elif not os.path.exists(mfile) and not os.path.exists(hfile):
        print(f"{m} missing {hfile}")
    
