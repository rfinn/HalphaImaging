#!/usr/bin/env python

'''

PROCEDURE:

rename images that have the filename format

VF-217.9524+7.6805-BOK-20210414-VFID5541-r.fits

to a format that includes ra and dec with leading zeros
b/c I screwed it up the first time...

plus I also included too many decimals on RA and DEC.  
Now using 3 instead of 4.

USAGE:

* running this from /media/rfinn/hdata/coadds/BOK2021pipeline

python ~/github/HalphaImaging/python3/BOK_fix_coadd_names_2.py

'''

import os
#import glob
import shutil

copyflag = False
if not copyflag:
    # change this for the appropriate input/output filelist
    outdir = '/media/rfinn/hdata/coadds/virgo-coadds-INT-2019-Feb/'
    outdir = '/media/rfinn/hdata/coadds/virgo-coadds-INT-2019-Jun/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)


# get list of current directory
flist1 = os.listdir()

# overwrite output files if they exist
overwrite = True
flist1.sort()
workingdir = os.getcwd()
for f in flist1:
    if os.path.isdir(f):
        #this will skp over subdirectories, etc
        continue
    if not f.find('INT') > -1:
        continue
    print(f)
    if f.find('+') >-1: # this means dec is positive
        t,radec,telescope,dateobs,pointing = f.split('-')[0:5]
        junk,filterwsuffix = f.split(pointing)
        ra,dec = radec.split("+")
    else: #dec is negative and is joined with a - instead of a +
        t,ra,dec,telescope,dateobs,pointing,filterwsuffix = f.split('-')[0:6]
        dec = -1*float(dec)
        junk,filterwsuffix = f.split(pointing)        
    ra = float(ra)
    dec = float(dec)
    # create string for output name
    if float(dec) < 0:
        outfile = 'VF-{:07.3f}-{:06.3f}-{:s}-{:s}-{:s}{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filterwsuffix)
    else:
        outfile = 'VF-{:07.3f}+{:06.3f}-{:s}-{:s}-{:s}{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filterwsuffix)

    if copyflag:
        print('renaming ',f,'->',outdir+outfile)
    else:
        print('renaming ',f,'->',outfile)

    # comment the following until you are sure that the output names are correct

    if copyflag:
        shutil.copy(f,os.path.join(outdir,outfile))
    else:
        os.rename(f,outfile)

    #os.chdir(workingdir)


