#!/usr/bin/env python

"""
move fVF*.fits from /data-pool/Halpha/coadds/virgo-coadds-MOSAIC/ to VF*.fits in /data-pool/Halpha/coadds/all-virgo-coadds/

"""
import shutil
import glob
import os

indir = '/data-pool/Halpha/coadds/virgo-coadds-MOSAIC/'
outdir = '/data-pool/Halpha/coadds/all-virgo-coadds/'

flist = glob.glob(indir+'fVF*.fits')
flist.sort()

for f in flist:
    outname = os.path.join(outdir,os.path.basename(f).replace('fVF','VF'))
    print(f"{f} -> {outname}")
    #shutil.copyfile(f,outname)


print('\n\n')
wlist = glob.glob(indir+'VF*weight.fits')
wlist.sort()

for w in wlist:
    outname = os.path.join(outdir,os.path.basename(w))    
    print(f"{w} -> {outname}")
    #shutil.copyfile(f,outname)
