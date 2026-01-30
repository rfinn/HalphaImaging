#!/usr/bin/env python
"""
* want to rerun getzp on all coadds b/c I've change how I am fitting ZP
* I have incorporated color correction from Matteo Fossati

* want to run all again to make sure all is done in a uniform way for all coadds

* this assumes that the coadds have the most recent naming convention, which includes the telescope and filter

OVERVIEW:

This is for running the first round of getzp with flattening on the coadds, after reconstructing
them in Dec 2025/Jan 2026.

the files are organized by the 3 observing runs, and the coadds are in subdirectories according to the pointing.

goal is to:
- flatten images in their current location, and 
- (separate program) then use script to grab the flattened files and copy them to the main coadd directory, while renaming them with the proper VFS coadd name.

USAGE:
- create an input file that contains a list of all the pointing subdirectory names
- then run in parallel

parallel 

"""
import glob
import os
import sys
import numpy as np

    
def runone(f):

        
    iinstrument = 'i'
    

    if 'Halpha' in f: #int 197
        ffilter = 'ha197'
    elif 'Ha6657' in f:
        ffilter = 'ha227'
    elif 'r.coadd' in f:
        ffilter = 'r'
    elif 'R.coadd' in f:
        ffilter = 'R'
    else:
        print('WARNING: did not recognize filter ',filter, f)
        sys.exit()


    
    # construct string
    getzpstring = f"python ~/github/HalphaImaging/python3/getzp.py --image {f} --instrument {iinstrument} --filter {ffilter}"
    #if instrument == 'BOK':
    #    getzpstring += ' --fixbok'
    print()
    print(f"Running getzp.py for file {f}")
    #print(getzpstring)

    # run getzp
    os.system(getzpstring)

# run one galaxy at a time, in parallel
# don't do this - it corrupted a bunch of images - need to figure out what I'm doing wrong in getzp.py
# might be that a store image as a temp file, and these get overwritten/accessed by multiple programs?
if len(sys.argv) > 1:
    topdir = os.getcwd()
    
    subdirname = sys.argv[1]
    os.chdir(subdirname)

    # get coadds
    filelist = glob.glob('*coadd*.fits')
    # call runone for each coadd
    for f in filelist:
        if f.find('weight') > -1:
            continue

        # skip CS images
        if f.find('CS.fits') > -1:
            continue
        if f.find('CS-ZP.fits') > -1:
            continue
        #print("calling runone for ",filename)
        runone(f)
    os.chdir(topdir)





