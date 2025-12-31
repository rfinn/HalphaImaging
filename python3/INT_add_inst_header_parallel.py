#!/usr/bin/env python

"""
OVERVIEW:
- this is a wrapper for INT_add_inst_header.py
- this is intended to be used with gnu parallel

USAGE:



PROCEDURE:
directory is passed on command line

move to directory

run INT_add_inst_header.py

move back to parent directory

"""
import os
import sys

# directory is passed on command line
subdir = sys.argv[1]

# move to directory
topdir = os.getcwd()

os.chdir(subdir)

# run INT_add_inst_header.py
os.system("python ~/github/HalphaImaging/python3/INT_add_inst_header.py"
# move back to parent directory
os.chdir(topdir)
