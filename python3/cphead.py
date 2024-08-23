#!/usr/bin/env python

'''
I inserted a subtract median sky step between scamp and swarp, 
and renamed the median-subtracted files mh*.

the scamp files are h*.head

so swarp doesn't associate the corrected headers


'''
import glob
import os
import argparse

parser = argparse.ArgumentParser(description ='Run sextractor, scamp, and swarp to determine WCS solution and make mosaics')
parser.add_argument('--filestring', dest = 'filestring', default = 'gm', help = 'string to prepend to the scamp h*.head files (default = "gm", which assumes  "h*o00.fits")')

args = parser.parse_args()


###########################################
# copy scamp files
###########################################

flist = glob.glob('h*.head')

for f in flist:
    command = "cp {} m{}".format(f,f)
    os.system(command)


