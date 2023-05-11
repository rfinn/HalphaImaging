#!/usr/bin/env python

'''

PROCEDURE:

* rename images that have the filename format

  VF-20210315-BOK-VFID2593-r.fits

  to a format that includes ra and dec

* get list of current files

USAGE:

* running this from /media/rfinn/hdata/coadds/BOK2021pipeline

'''

import os
import shutil
from astropy.io import fits
import glob




homedir = os.getenv("HOME")
# define directory for all coadds
telescope = 'BOK'
# get list of current directory
flist1 = os.listdir()

working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
workingdir = os.getcwd()
for f in flist1:
    if os.path.isdir(f):
        #this will skp over subdirectories, etc
        continue
    h = fits.getheader(f)
    ra = float(h['CRVAL1'])
    dec = float(h['CRVAL2'])

    t,dateobs,telescope,pointing,filterwsuffix = f.split('-')
    # create string for output name
    if float(dec) < 0:
        outfile = 'VF-{:07.3f}-{:06.3f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filterwsuffix)
    else:
        outfile = 'VF-{:07.3f}+{:06.3f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filterwsuffix)


    print('renaming ',f,'->',outfile)
    os.rename(f,outfile)
    #os.chdir(workingdir)


