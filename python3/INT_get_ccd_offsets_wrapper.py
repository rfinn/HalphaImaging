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

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ="This program measures the offset of each ccd relative to PanSTARRS, calculates the scaling required to bring the ccds onto a uniform scale, and adds a keyword PANSCALE to the image headers.  PANSCALE can be used as the flux scaling keyword with SWarp.")
    parser.add_argument('--dirlist', dest = 'dirlist', default = 'dirlist', help = 'file containing list of directories to run on.  default is dirlist')
    parser.add_argument('--filestring', dest = 'filestring', default = 'mWFC', help = 'Filestring to pass into INT_get_ccd_offsets.py.  Usually this is mWFC or WFC.')        
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'Will run on one directory only.')
    
    args = parser.parse_args()
    
    dirlist = open(args.dirlist, 'r')

    for d in dirlist:
        subdir = d.rstrip('\n')

        # move to directory
        topdir = os.getcwd()

        os.chdir(subdir)

        # run INT_add_inst_header.py
        os.system(f"python ~/github/HalphaImaging/python3/INT_get_ccd_offsets.py --filestring {args.filestring} --verbose")
        # move back to parent directory
        os.chdir(topdir)

        if args.testing:
            sys.exit()
