#!/usr/bin/env python

'''
GOAL:
- sort Bok 90Prime data according to target and filter
-

RELEVANT HEADER KEYWORDS:

    OBSTYPE: TARGET, SKY (skyflat), BIAS, FOCUS

    IMAGETYP: object, sky, focus, zero

    OBJECT: name that we assigned, i.e. pointing-021
    - need to strip spaces

    CATNAME:  similar to OBJECT, but capitalized? not all are the same.  use OBJECT

    WFFBAND - filter; need to strip spaces

    WFFPOS - filter position (can use this to double check)

PROCEDURE:

get list of filter, imagetype, and object names

for flats, create a directory for each filter called FLAT-HALPHA.  Move flat to appropriate directory

for bias, create BIAS, and move files there

for focus, create FOCUS directory, move files there

for science frames, sort them by filter into, e.g. SCIENCE-r or SCIENCE-Halpha
- these will get sorted by object once bias subtraction and flatfielding is done.

After basic calibration, we can sort science objects, make directory for each unique set of OBJECT-WFFBAND.  move files to appropriate directory.

'''
#import ccdproc
from ccdproc import ImageFileCollection
import os
import numpy as np
# updating 2022 to make this usable for BOK scripts, rather than making a separate

ic = ImageFileCollection(os.getcwd(),keywords='*',glob_include='data*.fits')

######################################################
### MOVE BIAS FRAMES
######################################################
if not os.path.exists('BIAS'):
    os.mkdir('BIAS')

# get bias files
bfiles = ic.files_filtered(obstype='zero')
for f in bfiles:
    print('move ',f,' -> BIAS/'+f)
    os.rename(f,'BIAS/'+f)

######################################################
### MOVE FLAT FRAMES
######################################################
# get list of filter, imagetype, and object names

# dome flats
filenames = ic.values('file')
filters = ic.values('filter')

unique_filters = set(filters)

# check for dome flats first
flattype = ['dome flat','sky flat']
flatnames = ['DOMEFLAT','SKYFLAT']
for i,d in enumerate(flatnames):
    for f in unique_filters:
        filt = f.replace(' ','-')
        dirname = d+'-'+filt
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        # get matching files
        flag = (ic.summary['obstype'] == flattype[i]) & (ic.summary['filter'] == f)
        # check for any matches with this combo of flat type and filter
        if np.sum(flag > 0):
            for ff in ic.summary['file'][flag]:
                print('move ',ff,' -> ',os.path.join(dirname,ff))# move files to subdir
                os.rename(ff,os.path.join(dirname,ff))# move files to subdir


######################################################
### MOVE OBJECT FRAMES
######################################################


# check for dome flats first
flattype = ['object']
flatnames = ['target']
for i,d in enumerate(flatnames):
    for f in unique_filters:
        filt = f.replace(' ','-')
        dirname = d+'-'+filt
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        # get matching files
        flag = (ic.summary['obstype'] == flattype[i]) & (ic.summary['filter'] == f)
        # check for any matches with this combo of flat type and filter
        if np.sum(flag > 0):
            for ff in ic.summary['file'][flag]:            
                print('move ',ff,' -> ',os.path.join(dirname,ff))
                os.rename(ff,os.path.join(dirname,ff))# move files to subdir




        




        
