#!/usr/bin/env python
"""

OVERVIEW:

This is for running the first round of getzp with flattening on the coadds, after reconstructing
them in Dec 2025/Jan 2026.

the files are organized by the 3 observing runs (2019 Feb, 2019 May, 2022 May), and the coadds are in subdirectories according to the pointing.

goal is to:
- flatten images in their current location, and 
- (separate program) then use script to grab the flattened files and copy them to the main coadd directory, while renaming them with the proper VFS coadd name.

X I am not running this in parallel because it corrupted a bunch of images - not sure why...  

USAGE:

python ~/github/HalphaImaging/python3/batch_getzp_INT2026.py 2019febINT-all

parallel --eta python ~/github/HalphaImaging/python3/batch_getzp_INT2026.py :::: dirlist


https://github.com/rfinn/HalphaImaging/wiki/06---INT-Data-and-Theli#run-getzp-to-flatten-images

"""
import glob
import os
import sys
import numpy as np

    
def runone(f):

        
    iinstrument = 'i'

    # updating in March 2026 - I don't understand what happened and why some of the resulting flatfields looked so bad.
    # this is particularly troublesome for Halpha
    #
    # when I rerun it now, the residual surfaces do not always match what is on the website - puzzling!
    # tested halpha, and one round of flattening with 3rd order spline works pretty well
    # I also used a minmag=14.5 for halpha (checked m=14, and this also works well
    # takes about ~60 sec to run one halpha field, but I can omit the diagnositic plot of putting residuals on image b/c
    # this is slow.

    # saving values from Feb 2026 run
    nflatten = 1
    usespline = True # using spline for r-band images too
    norder = 3
    if 'Halpha' in f: #int 197
        ffilter = 'ha197'
        nflatten = 2
        usespline = True
        norder = 4
    elif 'Ha6657' in f:
        ffilter = 'ha227'
        nflatten = 2
        usespline = True
        norder = 4
    elif 'r.coadd' in f:
        ffilter = 'r'
    elif 'R.coadd' in f:
        ffilter = 'R'
    else:
        print('WARNING: did not recognize filter ',filter, f)
        sys.exit()

    # Values for March 28
    nflatten = 1
    usespline = True # using spline for r-band images too
    norder = 3
    if 'Halpha' in f: #int 197
        ffilter = 'ha197'
        minmag = 14.5
    elif 'Ha6657' in f:
        ffilter = 'ha227'
        minmag = 14.5
    elif 'r.coadd' in f:
        ffilter = 'r'
        minmag = 15.
    elif 'R.coadd' in f:
        ffilter = 'R'
        minmag = 15.        
    else:
        print('WARNING: did not recognize filter ',filter, f)
        sys.exit()

        

    
    # construct string
    getzpstring = f"python ~/github/HalphaImaging/python3/getzp.py --image {f} --instrument {iinstrument} --filter {ffilter} --flatten {nflatten} --useastropy --minmag {minmag}"
    if usespline:
        getzpstring += " --spline --spline_order {norder}"
    #if instrument == 'BOK':
    #    getzpstring += ' --fixbok'
    print()
    print(f"Running getzp.py for file {f}")
    print(getzpstring)

    # run getzp
    os.system(getzpstring)

def runone_directory(subdirname):
    topdir = os.getcwd()
    print("working on subdir ",subdirname)
    os.chdir(subdirname)

    # get coadds
    if '2019' in os.getcwd():
        filelist1 = glob.glob('p*coadd*.fits')
        filelist2 = glob.glob('lmp*coadd*.fits')
        filelist = filelist1 + filelist2
    else:
        filelist = glob.glob('t*coadd*.fits')   
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

        # when rerunning second time, I just need to rerun on halpha - skip r?
        # going to run on both, just for completeness, symmetry, etc...
        
        runone(f)
    os.chdir(topdir)


# run one galaxy at a time, in parallel
# don't do this - it corrupted a bunch of images - need to figure out what I'm doing wrong in getzp.py
# might be that a store image as a temp file, and these get overwritten/accessed by multiple programs?
if len(sys.argv) > 1:
    topdir = os.getcwd()
    
    subdirname = sys.argv[1]
    runone_directory(subdirname)

else:  # assume we are starting in /data-pool/Halpha/processed
    rootdir = os.getcwd()
    rundirs = ["2019febINT-all", "2019nelvy-all", "2022MayINT-allfiles-v2"]

    for rundir in rundirs:
        os.chdir(os.path.join(rootdir, rundir))

        # get list of subdirectories
        subdirs = [d for d in os.listdir('.') if os.path.isdir(d)]

        for subdirname in subdirs:
            if subdirname.startswith(("lmpointing", "pointing", "target_")):
                print(f"Running {subdirname}")

                os.chdir(subdirname)

                try:
                    runone_directory(subdirname)
                except Exception as e:
                    print(f"ERROR in {subdirname}: {e}")

                # move back to rundir
                os.chdir('..')

    os.chdir(rootdir)
    




