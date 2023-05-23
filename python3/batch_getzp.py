#!/usr/bin/env python
"""
* want to rerun getzp on all coadds b/c I've change how I am fitting ZP
* I have incorporated color correction from Matteo Fossati

* want to run all again to make sure all is done in a uniform way for all coadds

* this assumes that the coadds have the most recent naming convention, which includes the telescope and filter

USAGE:

from /media/rfinn/hdata/coadds/all-virgo-coadds

python ~/github/HalphaImaging/python3/batch_getzp.py


"""
import glob
import os
import sys

# grab all files in the directory
rfiles = glob.glob('VF*.fits')
rfiles.sort()
for i,f in enumerate(rfiles):
    # skip weight files
    if f.find('weight') > -1:
        continue

    # skip CS images
    if f.find('CS.fits') > -1:
        continue
    # how do we skip the unshifted files?
    # this only seems to apply to INT data (why did I need to shift these, anyway???)
    # so if file has INT and -r
    if ('INT' in f) and ('r.fits' in f):
        continue
    
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
        

    halpha_names = ['ha4','Halpha','Ha6657','Ha4']
    if filter in halpha_names:
        ffilter = 'ha'
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
        iinstrument == 'm'
    else:
        print('WARNING: did not recognize instrument ',instrument, f)
        sys.exit()
    
    # construct string
    getzpstring = f"python ~/github/HalphaImaging/python3/getzp.py --image {f} --instrument {iinstrument} --filter {ffilter}"
    print()
    print(f"Running getzp.py for file {i}/{len(rfiles)}")
    print(getzpstring)

    # run getzp
    os.system(getzpstring)

    # uncomment the following for testing
    #if i > 0:
    #    break
        