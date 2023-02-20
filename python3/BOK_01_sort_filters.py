#!/usr/bin/env python

'''
GOAL:
- sort Bok 90Prime data according to image type and filter
- this sets up the files as needed for theli


PROCEDURE:

get list of filter, imagetype, and object names

for bias, create BIAS, and move files there

for flats, create a directory for each filter called FLAT-HALPHA.  Move flat to appropriate directory

for science frames, sort them by filter into, e.g. target-r or target-Halpha
- these will get sorted by object once bias subtraction and flatfielding is done.

After basic calibration, we can sort science objects, 
make directory for each unique set of OBJECT-FILTER.  
move files to appropriate directory.  But theli actually does this part...

'''
#import ccdproc
from ccdproc import ImageFileCollection
import os


def move_files(subdir_name,ic_sub):
    if not os.path.exists(subdir_name):
        os.mkdir(subdir_name)
    for f in ic_sub.values('file'):
        f = os.path.basename(f)
        try:
            os.rename(f,subdir_name+'/'+f)
        except:
            print("Error moving file {} to directory {}".format(f,subdir_name))


ic = ImageFileCollection(os.getcwd(),keywords='*',glob_include='*.fits')


# move bias frames to subdirectory
ic_bias = ic.filter(imagetyp='zero')
subdir_name = 'BIAS'

move_files(subdir_name,ic_bias)


# move flats to their own directory based on filter
filters = ['r','Ha4nm']
for filt in filters:
    ic_sub = ic.filter(imagetyp='flat',filter=filt)
    subdir_name = 'FLAT-'+filt
    move_files(subdir_name,ic_sub)
   
# move object files to a target-r or target-Ha4nm directory
filters = ['r','Ha4nm']
for filt in filters:
    ic_sub = ic.filter(imagetyp='object',filter=filt)
    subdir_name = 'target-'+filt
    move_files(subdir_name,ic_sub)





        






        
