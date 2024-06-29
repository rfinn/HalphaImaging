#!/usr/bin/python

"""
this is painful, but reproject doesn't seem to work well in parallel

so calling each image, one at a time :(

pass file with Ha images 
"""


import sys
import os


filelist = sys.argv[1]

infiles = open(filelist,'r')

for line in infiles:
    t = line.rstrip()
    print(t)
    s = f"python ~/github/HalphaImaging/python3/BOK_align_images_wrapper.py {t}"
    os.system(s)
    
