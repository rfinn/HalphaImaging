#!/usr/bin/env python

"""
run from the directory where the VFXXX mosaics are

"""
import glob
import os
import shutil

infiles = glob.glob('nVF*.fits')
infiles.sort()

print("Good news!")
print()
print(f"I found {len(infiles)} files to move!!!")
print()
# copy nVF to VF in the coadd directory


outdir = '/media/rfinn/hdata/coadds/BOK2021pipeline/'

for f in infiles:

    # copy and rename the coadded image
    outname = outdir+f.replace('nVF','VF')
    print(f"moving {f} -> {outname}")        
    shutil.copy(f,outname)

    # copy the weight file
    weightname = f.replace('nVF','VF').replace('.fits','.weight.fits')
    shutil.copy(weightname,outdir+weightname)
