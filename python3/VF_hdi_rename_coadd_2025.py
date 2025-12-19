#!/usr/bin/env python

'''
change names of coadds to something descriptive and unique.

DEC 2025:
- originally I named the coadds something like VF-2017-05-23-HDI-p012-R.fits, but then migrated them again to 

VF-RA-DEC-TEL-DATEOBS-POINTING-FILTER.fits

- instead of doing this in two steps like I did the first time around, I am updating this program to go right from coadd name to the final version.

USAGE:

from /mnt/qnap_home/rfinn/halpha-processed/raw/HDI/virgo-allfiles-2017/coadds-2025DEC

python ~/github/HalphaImaging/python3/VF_hdi_rename_coadd_2025.py --test --outdir /data-pool/Halpha/coadds-2025DEC/










OLD DESCRIPTION PRIOR TO DEC 2025:
The coadd names will include the year, instrument, date, pointing and filter.  This is necessary because I used the same pointing number to refer to different fields during different runs, meaning that for each run, I remade the finding charts starting pointing001, pointing002, etc.  We don't want to overwrite coadds, to the title needs to be unique.

```
VF-2017-05-23-HDI-p012-R.fits
```

The following script will rename the hdi coadds.  It gets the date from the equinox header, and then converts the equinox (which is in decimal years) to a string.  Note, the dates will be in UT.
```
python ~/github/HalphaImaging/python3/vf_hdi_rename_coadds.py
```



'''
import glob
import os
import argparse

from astropy.io import fits
from astropy.time import Time


def get_updated_coadd_name(imname, survey='VF'):
    h = fits.getheader(f)

    ############################################
    ## GET DATE
    ############################################    
    
    # store time 
    t = h['EPOCH']
    # convert to year, month,day
    t = Time(t,format='decimalyear')
    dateobs = t.iso.split()[0]

    #dateobs_time = h['OBSDATE']
    #dateobs = dateobs_time.split('T')[0]
    dateobs = dateobs.replace('-','')
    #print(dateobs)    


    ############################################
    ## GET POINTING
    ############################################    
    
    # read in object
    o = h['OBJECT'] # should split into "pointing" and "10" for example
    # try to identify format
    #print(o,len(o))
    if len(o.split()) > 1:
        split_string=' '
        #print('object names contain ',split_string)        
    elif len(o.split('-')) > 1:
        split_string='-'
        #print('object names contain ',split_string)
    elif len(o.split('_')) > 1:
        split_string='_'
        #print('object names contain ',split_string)        
    
             
             
    if (o.find('lm') > -1)| (o.find('LM') > -1):
        # low-mass pointing
        prefix = "lmp"
    else:
        prefix = "p"
    pointing = "{}{:03d}".format(prefix,int(o.split(split_string)[1]))

    ############################################
    ## GET FILTER
    ############################################        
    # read in filter
    ffilter = h['FILTER']

    ############################################
    ## NOT RELEVANT ANYMORE B/C I REMOVED NOBACK IN 2025
    ############################################        
    # no back - does this mean no background,
    # or no background subtraction???
    if f.find('noback') > -1:
        noback='noback-'
    else:
        noback=''

    ############################################
    ## CHECK NAME FOR WEIGHT
    ############################################        

    # check to see if this is a weight image
    if f.find('weight') > -1:
        suffix='coadd.weight.fits'
    else:
        suffix='coadd.fits'
    #outfile = "VF-{}-{}-{}-{}-{}".format(dateobs,instrument,pointing,ffilter,suffix)
  

    ############################################
    ## GET RA AND DEC IN DEG
    ############################################        

    #h = fits.getheader(imname)
    ra = float(h['CRVAL1'])
    dec = float(h['CRVAL2'])
    

    ############################################
    ## GET INSTRUMENT
    ############################################    
    
    
    temptel = h['INSTRUME']
    if (temptel.find('hdi') > -1) | (temptel.find('HDI') > -1):
        telescope = "HDI"
    elif (temptel.find('mos') > -1) | (temptel.find('MOS') > -1):
        telescope = "MOS"
    else:
        print("WARNING: could not parse header field INSTRUME")
        print("\t setting instrument to HDI")

    # let's not switch everything to underscores b/c it's wreaking havoc with halphagui
    p,filterwithsuffix = imname.split('_')
    #pointing = pointing.replace('-','_')
    # remove coadd from name
    filterwithsuffix = filterwithsuffix.replace('.coadd','')
    # create string for output name

    if float(dec) < 0:
        #outfile = f'UAT_{ra:07.3f}-{dec:06.3f}_{telescope}_{dateobs}_{pointing}_{filterwithsuffix}'
        outfile = f'{survey}-{ra:07.3f}-{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'        
    else:
        #outfile = f'UAT_{ra:07.3f}+{dec:06.3f}_{telescope}_{dateobs}_{pointing}_{filterwithsuffix}'
        outfile = f'{survey}-{ra:07.3f}+{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'        
    return outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Rename HDI coadds to VF format.  EX: VF-2017-05-23-HDI-p012-R-coadd.fits')


    parser.add_argument('--filestring', dest='filestring', default='pointing', help='filestring that points to files to rename.  glob will pull filestring*.fits')
    parser.add_argument('--indir', dest = 'indir', default = '.', help = 'directory of input images.  Default is current directory')
    parser.add_argument('--outdir', dest = 'outdir', default = '.', help = 'directory to write output images to.  Default is current directory')
    
    parser.add_argument('--test', dest='test', default=False, action='store_true', help='set this to print out new filenames but not actually rename the files')
    

    args = parser.parse_args()

    filelist = glob.glob(os.path.join(args.indir,"*"+args.filestring+"*.fits"))
    filelist.sort()

    for f in filelist:
        new_name = get_updated_coadd_name(f)
        new_output_image = os.path.join(args.outdir, new_name)        
        # read in header
        print('moving ',f,' -> ',new_output_image)
        if not args.test:
            os.system(f"cp {f} {new_output_image}")
