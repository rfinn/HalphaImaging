#!/usr/bin/env python

'''

PROCEDURE:
* get list of current file

I ran the INT_scamp_swarp before merging all filters from the same pointing into one folder.

I then made the coadds, and the files ended up with -r_r or -Halpha_Halpha in their names.

just writing a quick script to get rid of the double filter names
'''

import os
import shutil
from astropy.io import fits
import glob


import argparse
import subprocess


homedir = os.getenv("HOME")
# define directory for all coadds
output_dir_coadds = args.outdir
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()

working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
workingdir = os.getcwd()
for subdir in flist1:
    if os.path.isdir(subdir) & (subdir.find('pointing')>-1 ):# & (subdir.find('-') > -1):
        print("I'm in "+subdir+"!!!")
        os.chdir(subdir)
        # store pointing and filter
        #long_pointing,filter = subdir.split('-')
        # redoing for when both filters are in the same directory
        long_pointing = subdir
        pointing = long_pointing.replace('ointing','')
        

        ## GRAB THE COADDS IN THIS DIRECTORY
                

        filesuffix = ['-r_r','-Halpha_Halpha','-Ha6657_Ha6657']
        newsuffix = ['_r','_Halpha','_Ha6657']
        for i,fsuffix in enumerate(filesuffix):
            files = glob.glob('*{}*.fits'.format(fsuffix))
            for f in files:
                outname = f.replace(fsuffix,newsuffix[i])
                os.rename(f,outname)


