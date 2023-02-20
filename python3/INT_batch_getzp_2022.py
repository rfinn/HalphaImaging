#!/usr/bin/env python

'''
GOAL: 
* go into each subdirectory target-r_0/coadd-r/ and run getzp on coadded image coadd.fits
* this assumes that theli has been used through coaddition

USAGE:
* Run this from, e.g. /media/rfinn/ssd4t/rfinn/data/INT/2022-allfiles-v2
  - this directory has a subdirectory for each pointing

* to run type

  python ~/github/HalphaImaging/python3/INT_batch_getzp_2022.py

NOTES:
* this is a slight variation on INT_batch_getzp.py
* it assumes directory structure created by theli, using theli to split objects and create coadds
'''

import os
import shutil
import glob
from astropy.io import fits
import matplotlib

matplotlib.use("Qt5agg")

homedir = os.getenv("HOME")

telescope = 'INT'

# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()

subdir_rootname = 'target'
subdir_rootname = 'pointing'

#print(flist1)
#flist1 = ['pointing022','pointing026']
for subdir in flist1: # loop through list

    #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
    if os.path.isdir(subdir) & (subdir.startswith('target-r_') | subdir.startswith('target-Halpha_')):
            
        if subdir.startswith('target-r'):
            filter = 'r'
            getzp_filter = 'r'
        else:
            filter = 'Halpha'
            getzp_filter = 'ha'
        # look for coadd subdir

        coadd_subdir = subdir+'/coadd_'+filter+'/'
        print(coadd_subdir)
        if os.path.exists(coadd_subdir):
            print()
            print('##########################################')
            print('##########################################')        
            print('WORKING ON DIRECTORY: ',coadd_subdir)
            print('##########################################')
            print('##########################################')
            
            os.chdir(coadd_subdir)
            try:
                os.system('python ~/github/HalphaImaging/python3/getzp.py --instrument i --image coadd.fits --filter '+getzp_filter)
            except:
                print('##########################################')
                print('WARNING: problem running getzp for ',coadd_subdir)
                print('##########################################')

        os.chdir(working_dir)
        # just running on one directory for testing purposes
        #break



