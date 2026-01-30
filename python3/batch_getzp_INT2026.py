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
    if f.find('weight') > -1:
        return

    # skip CS images
    if f.find('CS.fits') > -1:
        return
    if f.find('CS-ZP.fits') > -1:
        return
    # how do we skip the unshifted files?
    # this only seems to apply to INT data (why did I need to shift these, anyway???)
    # so if file has INT and -r
    if ('INT' in f) and ('r.fits' in f):
        return
    
    if ('BOK' in f) and ('r.fits' in f):
        return
    
    # determine telescope and filter from filename
    t = f.split('-')
    #print(t)
    # check to see if declination is negative
    # if it is, then index of instrument will be off by one
    if f[10] == '+':
        instrument = t[2]
        filter = t[5].replace('.fits','')
    else:
        instrument = t[3]
        filter = t[6].replace('.fits','')
        

    halpha_filters = ['ha4','Halpha','Ha6657','Ha4']
    if filter in ['ha4','Ha4','Ha+4nm']:
        ffilter = 'ha4'
    elif filter == 'Halpha': #int 197
        ffilter = 'ha197'
    elif filter == 'Ha6657'
        ffilter = 'ha227'
    elif filter == 'r':
        ffilter = 'r'
    elif filter == 'R':
        ffilter = 'R'
    else:
        print('WARNING: did not recognize filter ',filter, f)
        sys.exit()

    # get instrument
    if instrument == 'INT':
        iinstrument = 'i'
    elif instrument == 'BOK':
        iinstrument = 'b'
    elif instrument == 'HDI':
        iinstrument = 'h'
    elif instrument == 'MOSAIC':
        iinstrument = 'm'
    elif instrument == 'MOS': # didn't have this case for mosaic images
        iinstrument = 'm'
    else:
        print('WARNING: did not recognize instrument ',instrument, f)
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
    for fname in filelist:
        if f.find('weight') > -1:
            continue

        # skip CS images
        if f.find('CS.fits') > -1:
            continue
        if f.find('CS-ZP.fits') > -1:
            continue
        print("calling runone for ",filename)
        runone(filename)
    os.chdir(topdir)





