#!/usr/bin/env python

'''

PROCEDURE:

* rename images to a uniform format
* files have format

  VF-2017-05-23-HDI-p010-R-noback-coadd.fits

  and need to have

  VF-208.704+40.669-BOK-20210417-VFID2068-r.fits 

* get list of current files

USAGE:

* running this from /home/rfinn/data/reduced/virgo-coadds-HDI/


'''

import os
import shutil
from astropy.io import fits
import glob


outdir = '/media/rfinn/hdata/coadds/virgo-coadds-HDI/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

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

    t,year,month,day,telescope,pointing,filter, junk1,suffix = f.split('-')
    # create string for output name
    dateobs = year+month+day

    # remove coadd from suffix

    #print(suffix)
    #print(suffix.replace('coadd',''))
    
    if float(dec) < 0:
        outfile = 'VF-{:07.3f}-{:06.3f}-{:s}-{:s}-{:s}-{:s}{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter, suffix.replace('coadd',''))
    else:
        outfile = 'VF-{:07.3f}+{:06.3f}-{:s}-{:s}-{:s}-{:s}{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter,suffix.replace('coadd',''))


    print('renaming ',f,'->',outdir+outfile)
    shutil.copy(f,os.path.join(outdir,outfile))
    #os.rename(f,outfile)



