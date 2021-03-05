#!/usr/bin/env python

'''
I inserted a subtract median sky step between scamp and swarp, 
and renamed the median-subtracted files mh*.

the scamp files are h*.head

so swarp doesn't associate the corrected headers


'''
import glob
import os

flist = glob.glob('h*.head')

for f in flist:
    command = "cp {} m{}".format(f,f)
    os.system(command)
