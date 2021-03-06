#!/usr/bin/env python

'''
change names of coadds to something descriptive and unique.

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


parser = argparse.ArgumentParser(description ='Rename HDI coadds to VF format.  EX: VF-2017-05-23-HDI-p012-R-coadd.fits')


parser.add_argument('--filestring', dest='filestring', default='pointing', help='list of image sets to run swarp on.  the file should contain the list of all Rband groups, for example: ls pointing*_R > swarp_input.  This will look for the corresponding list of halpha images.')
parser.add_argument('--test', dest='test', default=False, action='store_true', help='set this to print out new filenames but not actually rename the files')

args = parser.parse_args()

flist = glob.glob(args.filestring+'*.fits')
flist.sort()

for f in flist:
    # read in header
    header = fits.getheader(f)
    
    # store time 
    t = header['EPOCH']
    # convert to year, month,day
    t = Time(t,format='decimalyear')
    dateobs = t.iso.split()[0]
    
    # read in instrument
    i = header['INSTRUME']
    # convert to capitals
    instrument = i.upper()
    
    # read in object
    o = header['OBJECT'] # should split into "pointing" and "10" for example
    # try to identify format
    print(o,len(o))
    if len(o.split()) > 1:
        split_string=' '
        print('object names contain ',split_string)        
    elif len(o.split('-')) > 1:
        split_string='-'
        print('object names contain ',split_string)
    elif len(o.split('_')) > 1:
        split_string='_'
        print('object names contain ',split_string)        
    
             
             
    if (o.find('lm') > -1)| (o.find('LM') > -1):
        # low-mass pointing
        prefix = "lmp"
    else:
        prefix = "p"
    pointing = "{}{:03d}".format(prefix,int(o.split(split_string)[1]))
    # read in filter
    ffilter = header['FILTER']

    # no back
    if f.find('noback') > -1:
        noback='noback-'
    else:
        noback=''
    # check to see if this is a weight image
    if f.find('weight') > -1:
        suffix=noback+'coadd.weight.fits'
    else:
        suffix=noback+'coadd.fits'
    outfile = "VF-{}-{}-{}-{}-{}".format(dateobs,instrument,pointing,ffilter,suffix)
    print('moving ',f,' -> ',outfile)
    if not args.test:
        os.rename(f,outfile)
