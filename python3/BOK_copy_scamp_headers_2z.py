#!/usr/bin/env python

"""
GOAL:
copy the scamp header files so they have the 'z' prefix

this should have happened automatically in BOK_fixampoffsets.py but didn't

"""

import glob
import shutil
import os

filelist = glob.glob('mksb*.head')
filelist.sort()

for scamp_header in filelist:
    # check to see if header files exist
    # if they do, then copy to have z prefix

    if os.path.exists(scamp_header):
        shutil.copy(scamp_header,'z'+scamp_header)
