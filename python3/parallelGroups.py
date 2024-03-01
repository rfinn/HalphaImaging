#!/usr/bin/env python

"""
GOAL: launch HDI reduction for all nights' data from a given run simultaneously.  I don't know how much is memory intensive vs processor intensive, but it can't hurt, right?



USAGE:
Feed data directory on command line, e.g.:

python ~/github/HalphaImaging/python3/parallelGroups.py 20140423


To run in parallel:

* create a file that is a list of the subdirectories

* then type:


parallel --eta python ~/github/HalphaImaging/python3/parallelGroups.py :::: dirlist.txt
"""
import os
import sys

# save the starting directory; return to this at the end
startdir = os.getcwd()

# move to directory
dirname = sys.argv[1]
os.chdir(dirname)

# run uat_mvjunk.py
os.system('python ~/github/HalphaImaging/python3/uat_mvjunk.py')

# run processHDI.py
os.system('python ~/github/HalphaImaging/processHDI.py --trim --zap --groupflat --flatwdome --fixheader')

os.chdir(startdir)
