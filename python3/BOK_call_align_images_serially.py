#!/usr/bin/python

"""
this is painful, but reproject doesn't seem to work well in parallel

so calling each image, one at a time :(

pass file with Ha images 
"""


import sys
import os


filelist = sys.argv[1]

# get number of lines
t = open(filelist,'r')
nlines = len(t.readlines())
t.close()


infiles = open(filelist,'r')
for i,line in enumerate(infiles):
    print()
    print(f"good news - running {i+1} out of {nlines}")
    t = line.rstrip()
    print(t)
    s = f"python ~/github/HalphaImaging/python3/BOK_align_images_wrapper.py {t}"
    os.system(s)

infiles.close()
