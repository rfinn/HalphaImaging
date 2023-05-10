#!/usr/bin/env python

'''

PROCEDURE:
* gather coadds for a given pointing into the same directory

* this is an updated version of INT_gathercoadds.py
  - previously, the object names were pointing20, for example, whereas in 2022 they 
    are named for the primary target according to its VFID, e.g. VFID1010

* in 2022, I used theli to split targets, and that made the naming conventions differ
  - theli makes target-r_0, target-r_1, etc

* the goal of this script is to get the coadd in 
  - target-r_0/coadd-r/coadd.fits and t
  - target-Halpha_0/coadd-Halpha/coadd.fits
  - and move them to VFID3054/


* run this from the base data directory, like ~/data/INT
  - this have subfolders arranged by date

* it then looks in targetXXX-r/coadd-r for the coadd.fits file
  - it will copy that to the output_dir_coadds directory specified below

* it will repeat for Halpha

* coadds will be renamed by RA, DEC, telescope, pointing, and filter


USAGE:
* run after INT_batch_getzp_2022.py
  - assumes that coadds will be named fcoadd.fits after flattening from getzp.py

* run in /media/rfinn/ssd4t/rfinn/data/INT/2022-allfiles-v2/


NEXT UP:

* align the r-band image to the Halpha image
  - INT_align_images.py
  - INT_batch_align_images.py

* then feed into halphgui!


'''

import os
import shutil
from astropy.io import fits
from astropy.time import Time

homedir = os.getenv("HOME")
# define directory for all coadds
output_dir_coadds ='/media/rfinn/hdata/coadds/virgo-coadds-int-2022/'
if not os.path.exists(output_dir_coadds):
    os.mkdir(output_dir_coadds)
telescope = 'INT'


# get list of current directory
# I am expecting subdirectories named target-r_0, target-Halpha_0, etc...
flist1 = os.listdir()

working_dir = os.getcwd()

# overwrite output files if they exist
overwrite = True
flist1.sort()
print(flist1)
for subdir in flist1:
    # if item is a directory and the name contains pointing, then assume it is a target
    if os.path.isdir(subdir) & subdir.startswith('target-r_'):

        # the 2022 data don't need the extra flattening so need to remove the f from the image name
        coadd_r = os.path.join(subdir,'coadd_r/fcoadd.fits')
        weight_r = os.path.join(subdir,'coadd_r/coadd.weight.fits')        
                                
        coadd_ha = os.path.join(subdir.replace('-r','-Halpha'),'coadd_Halpha/ffcoadd.fits')
        weight_ha = os.path.join(subdir.replace('-r','-Halpha'),'coadd_Halpha/coadd.weight.fits')
        
        # this is directory structure setup by theli
        # when I processed them myself, the coadds are in the main directory
        if os.path.exists(coadd_r):
            print("found ",coadd_r)
            if os.path.exists(coadd_ha):
                # create string for output name
                images = [coadd_r, coadd_ha]
                weights = [coadd_r, coadd_ha]
                filters = ['r','Halpha']
                for i,filter in enumerate(filters):
                    h = fits.getheader(images[i])
                    ra = float(h['CRVAL1'])
                    dec = float(h['CRVAL2'])
                    try:
                        full_dateobs = h['DATEOBS1']
                        dateobs = full_dateobs.split('T')[0].replace('-','')
                    except KeyError:
                        # assume the header was "fixed" to contain MJD-OBS
                        mjd_obs = Time(h['MJD-OBS'],format='mjd')
                        full_dateobs = mjd_obs.isot
                        dateobs = full_dateobs.split('T')[0].replace('-','')                        
                    pointing = h['OBJECT']
                    
                    if float(dec) < 0:
                        outfile = output_dir_coadds+'VF-{:.4f}-{:.4f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter)
                    else:
                        outfile = output_dir_coadds+'VF-{:.4f}+{:.4f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter)

                    # copy imfile to outfile
                    out1 = outfile+'.fits'
                    out2 = outfile+'.weight.fits'
                    if (not os.path.exists(out1)) or overwrite:
                        print('\t   copy ',images[i],' -> ',out1)
                        shutil.copyfile(images[i],out1)
                    else:
                        print('\t   '+out1,' already exists. set overwrite if you want to copy anyway.')
                    if (not os.path.exists(out2)) or overwrite:
                        # copy weight file
                        shutil.copyfile(weights[i],out2)                    
                        print('\t   weight file = ',weights[i])
            else:
                print("\twarning!  did not find halpha image ",coadd_ha)
                        



