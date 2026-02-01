#!/usr/bin/env python

'''

PROCEDURE:
* run from draco:/data-pool/Halpha/processed

* will look in ['2019febINT-all', '2019nelvy-all','2022MayINT-allfiles-v2']

* program will look in each run folder, and then in each sub directory that has a "pointing" or "target_" in the name

* it will copy that to the output_dir_coadds directory specified below

* coadds will be renamed by RA, DEC, telescope, pointing, and filter

'''

import os
import shutil
from astropy.io import fits
import glob

homedir = os.getenv("HOME")
# define directory for all coadds

telescope = 'INT'

rundirs = ['2019febINT-all', '2019nelvy-all','2022MayINT-allfiles-v2']





if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description ="This program copies coadd and weight files to output directory, renaming according to VFS convention.")
    #parser.add_argument('--filestring', dest = 'filestring', default = 'UAT', help = 'filestring to match. default is UAT, which will grab all UAT*R.fits files. Note: weight images will be renamed with their corresponding science coadd.')
    #parser.add_argument('--se', dest = 'se', default = False, action='store_true', help = 'Run source extractor')    
    #parser.add_argument('--scamp', dest = 'scamp', default = False, action='store_true', help = 'Run scamp')
    #parser.add_argument('--submedian', dest = 'submedian', default = False, action='store_true', help = 'Subtract a median sky value from each image.')
    parser.add_argument('--outdir', dest = 'outdir', default = '/data-pool/Halpha/coadds-2025DEC', action='store_true', help = 'output directory to copy renamed coadds to.  default is /data-pool/Halpha/coadds-2025DEC ')    
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'Will print matches but will not edit headers.')
    
    args = parser.parse_args()


    output_dir_coadds = args.outdir
    topdir = os.getcwd()
    for rd in rundirs:

        os.chdir(rd)

        # get list of directory, including lm/pointing/target_ subdirs
        flist1 = os.listdir()

        # overwrite output files if they exist
        overwrite = True
        flist1.sort()

        thisrundir = os.getcwd() # save for returning after each pointing

        for f1 in flist1:

            # check this is a directory and a pointing
            if (os.path.isdir(f1)) & (f1.startswith('pointing') | f1.startswith('lmpointing') | f1.startswith('target_')):
                os.chdir(f1)
                print('WORKING ON DIRECTORY: ',f1)

                # get list of coadds
                coaddlist = glob.glob('f*coadd.fits')
                coaddlist.sort()

                for coadd in coaddlist:

                    # get info from header
                    cheader = fits.getheader(coadd)

                    # construct output filename
                    ra = float(cheader['CRVAL1'])
                    dec = float(cheader['CRVAL2'])

                    thisdate = cheader['DATE-OBS']
                    t = thisdate.split('T')
                    dateobs = t[0].replace('-','')

                    if '2022' in rd:
                        pointing = cheader['OBJECT']#.replace('lm_pointing_','lmp').replace('pointing-','p')
                    else:
                        pointing = f1

                    tfilter = cheader['FILTER']
                    if float(dec) < 0:
                        outfile = os.path.join(output_dir_coadds,'VF-{:07.3f}-{:06.3f}-{:s}-{:s}-{:s}-{:s}.fits'.format(ra,abs(dec),telescope,dateobs,pointing,tfilter))
                    else:
                        outfile = os.path.join(output_dir_coadds,'VF-{:07.3f}+{:06.3f}-{:s}-{:s}-{:s}-{:s}.fits'.format(ra,abs(dec),telescope,dateobs,pointing,tfilter))

                    # check for weight file
                    weightfile = coadd.replace('coadd.fits','coadd.weight.fits')
                    weightflag = False
                    if os.path.exists(weightfile):
                        outweightfile = outfile.replace('.fits','.weight.fits')
                        weightflag = True
                        
                    # copy file to output coadd directory
                    if args.testing:
                        print(f"testing: copying {coadd} {outfile}")
                        if weightflag:
                            print(f"testing: copying {coadd} {outfile}")                            
                    else:
                        shutil.copyfile(coadd,outfile)
                        # copy weight file to output coadd directory
                        if weightflag:
                            shutil.copyfile(weightfile,outweightfile)                    
                    







