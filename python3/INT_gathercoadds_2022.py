#!/usr/bin/env python

'''

PROCEDURE:
* gather coadds for a given pointing into the same directory

* this is an updated version of INT_gathercoadds.py
  - previously, the object names were pointing20, for example, whereas in 2022 they 
    are named for the primary target according to its VFID, e.g. VFID1010 (v1 catalogs)

* in 2022, I used theli to split targets, and that made the naming conventions differ
  - theli makes target-r_0, target-r_1, etc

X* the goal of this script is to get the coadd in 
X  - target-r_0/coadd-r/coadd.fits and t
X  - target-Halpha_0/coadd-Halpha/coadd.fits
X  - and move them to VFID3054/

* Updating 2023-06-28 to run on draco
  - the data were copied there under /data-pool/laptop-backup/rfinn/data/INT/2022-allfiles-v2

* run this from the base data directory, like ~/data/INT/2022-allfiles-v2/
  - this have subfolders arranged by date

* it then looks in targetXXX-r/coadd-r for the coadd.fits file
  - it will copy that to the output_dir_coadds directory specified below

* it will repeat for Halpha

* coadds will be renamed by RA, DEC, telescope, pointing, and filter
  - uses RA and DEC from Halpha image

USAGE:

* NO LONGER TRUE- I run getzp in the final mosaic directory so that all the plots for all the coadds
  are in the same place
X  run after INT_batch_getzp_2022.py
X  - assumes that coadds will be named fcoadd.fits after flattening from getzp.py

* UPDATE: on draco, data are in /data-pool/laptop-backup/rfinn/data/INT/2022-allfiles-v2
X run in /media/rfinn/ssd4t/rfinn/data/INT/2022-allfiles-v2/


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
output_dir_coadds ='/data-pool/Halpha/coadds/virgo-coadds-int-2022-v3/'
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


# creating filelist by hand because we are only missing a few galaxies
flist1 = ['target-r_5','target-r_6','target-r_7','target-r_8','target-r_9']
flist1 = ['target-r_7a','target-r_7b']
print(flist1)


for subdir in flist1:
    # if item is a directory and the name contains pointing, then assume it is a target
    if os.path.isdir(subdir) & subdir.startswith('target-r_'):
        print("got here!")
        # the 2022 data don't need the extra flattening so need to remove the f from the image name
        # looks like fcoadd.fits does not exist for the missing coadds, so changing this to coadd.fits
        # can deal with flattening once the coadds are gathered up
        coadd_r = os.path.join(subdir,'coadd_r/coadd.fits')
        weight_r = os.path.join(subdir,'coadd_r/coadd.weight.fits')        
                                
        coadd_ha = os.path.join(subdir.replace('-r','-Halpha'),'coadd_Halpha/coadd.fits')
        weight_ha = os.path.join(subdir.replace('-r','-Halpha'),'coadd_Halpha/coadd.weight.fits')
        
        # this is directory structure setup by theli
        # when I processed them myself, the coadds are in the main directory
        print(coadd_r)
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

                    ##
                    # TODO - make sure I'm using the RA and DEC from the Halpha image
                    # because the r-band image will be aligned to that and could have a different RA and DEC
                    #
                    # this is an update after moving to draco
                    ##
                    if float(dec) < 0:
                        outfile = output_dir_coadds+'VF-{:.3f}-{:.3f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter)
                    else:
                        outfile = output_dir_coadds+'VF-{:.3f}+{:.3f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filter)

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
                        



