#!/usr/bin/env python

'''
GOAL:
- sort INT WFC data according to target and filter
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
#from ccdproc import ImageFileCollection
import os
import glob
from astropy.io import fits
import sys
# this is not working
# seems like something is wrong with DATE-OBS card, but don't know what

#ic = ImageFileCollection(os.getcwd(),keywords=['FILTER','OBJECT'],glob_include='r*.fit*')

filters=[]
object_names=[]
infiles = glob.glob('r*.fit*')
if len(infiles) ==0:
    infiles = glob.glob('WFC*.fits')
    if len(infiles) == 0:
        print("no files found that match r*.fit* or WFC*.fits")
        sys.exit()
for f in infiles:
    print(f)
    d,h = fits.getdata(f,header=True)
    filters.append(h['FILTER'])
    object_names.append(h['OBJECT'])

    
# get list of filter, imagetype, and object names

#filters = ic.values('wffband')
#object_names = ic.values('object')
#filters = ic.values('FILTER')
#object_names = ic.values('OBJECT')

# create a unique set of object-filter
filefilterlist = ["{}-{}".format(str(i).replace(" ","").replace("-","").replace("_",""),str(j).replace(" ","")) for i,j in zip(object_names,filters)]

# find bias and change the name to BIAS (without filter)
# replace Skybackground with SKYFLAT
for i,d in enumerate(filefilterlist):
    if d.find('Bias') > -1:
        filefilterlist[i] = 'BIAS'
    if d.find('Skybackground') > -1:
        print('found it')
        filefilterlist[i] = filefilterlist[i].replace("Skybackground", "SKYFLAT")

# keep only unique directory names
dirnames = set(filefilterlist)



# make subdirectory for each object-filter combination
for d in dirnames:
    print(d)
    if d.find('--') > -1:
        continue
    if os.path.exists('../'+d):
        continue
    else:
        os.mkdir('../'+d)
        #os.mkdir('../'+d+'/MASK_IMAGES')
# move files to appropriate subdirectory
for f,d in zip(infiles,filefilterlist):
    try:
        os.rename(f,'../'+d+'/'+f)
        # move SE cat if it exists
        secat = f.split('.fits')[0]+'.cat'
        try:
            os.rename(secat,'../'+d+'/'+secat)
        except:
            pass
        scampfile = f.split('.fits')[0]+'.head'
        try:
            os.rename(scampfile,'../'+d+'/'+scampfile)
        except:
            pass
            
    except:
        print("Error moving file {} to directory {}".format(f,d))

        




        
