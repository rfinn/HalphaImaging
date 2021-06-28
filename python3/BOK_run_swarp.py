#!/usr/bin/env python

'''


'''

import ccdproc as ccdp
import os

# get list of images (images only)


# sort images by location and filter
# alternatively could use object name,
# but not all are correct, so need to fix names


# create file list with r-band images
# create a file list with r-band weight images


# create file list with halpha images
# create a file list with halpha weight images


# run swarp to mosaic r-band

# run swarp to mosaic halpha

# alight r and halpha imaging


def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i+1
    except UnboundLocalError:
        return 0


homedir = os.getenv("HOME")
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
runscamp=False
runswarp=True


rawdir = os.getcwd()

keys = ['OBJECT', 'RA', 'DEC', 'FILTER', 'EXPTIME']
ic = ccdp.ImageFileCollection(rawdir, keywords=keys, glob_include='*.fit*', glob_exclude='*test*.fit')

