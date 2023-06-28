#!/usr/bin/env python
"""
* want to rerun getzp on all coadds b/c I've change how I am fitting ZP
* I have incorporated color correction from Matteo Fossati

* want to run all again to make sure all is done in a uniform way for all coadds

* this assumes that the coadds have the most recent naming convention, which includes the telescope and filter

USAGE:

from /media/rfinn/hdata/coadds/all-virgo-coadds

python ~/github/HalphaImaging/python3/batch_getzp.py

To run on a subset of files, give a match string VF*matchstring*.fits:

  python ~/github/HalphaImaging/python3/batch_getzp.py BOK

"""
import glob
import os
import sys
import numpy as np
import multiprocessing as mp

image_results = []
def collect_results(result):

    global results
    image_results.append(result)


    
def runone(i,f):
    if f.find('weight') > -1:
        return

    # skip CS images
    if f.find('CS.fits') > -1:
        return
    # how do we skip the unshifted files?
    # this only seems to apply to INT data (why did I need to shift these, anyway???)
    # so if file has INT and -r
    if ('INT' in f) and ('r.fits' in f):
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
    #if instrument == 'BOK':
    #    getzpstring += ' --fixbok'
    print()
    print(f"Running getzp.py for file {i}/{len(rfiles)}")
    #print(getzpstring)

    # run getzp
    os.system(getzpstring)

# grab all files in the directory
if len(sys.argv) > 1:
    matchstring = sys.argv[1]
    rfiles = glob.glob('VF*'+matchstring+'*.fits')
else:
    rfiles = glob.glob('VF*.fits')
rfiles.sort()
indices = np.arange(len(rfiles))

#for rimage in rimages: # loop through list
image_pool = mp.Pool(mp.cpu_count())
myresults = [image_pool.apply_async(runone,args=(i,rfiles[i]),callback=collect_results) for i in indices]
    
image_pool.close()
image_pool.join()
image_results = [r.get() for r in myresults]


