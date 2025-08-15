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

parser = argparse.ArgumentParser(description ='Copy the header files with a new prefix to match the further processed images so that swarp correctly associates the header files with the image files.')
parser.add_argument('--filestring', dest = 'filestring', default = 'h', help = 'string to use to get input files (default = "h", which grabs all of the files "h*.head")')
parser.add_argument('--addprefix', dest = 'addprefix', default = 'gm', help = 'string to prepend to the scamp h*.head files (default = "gm", which assumes  "h*o00.fits")')

args = parser.parse_args()


###########################################
# copy scamp files
###########################################

flist = glob.glob(args.filestring+'*.head')

for f in flist:
    command = f"cp {f} {args.addprefix}{f}"
    os.system(command)


