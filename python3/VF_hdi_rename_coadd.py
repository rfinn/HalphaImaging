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
    if (o.find('lm') > -1)| (o.find('LM') > -1):
        # low-mass pointing
        pointing = "lmp{:03d}".format(float(o.split()[1]))
    else:
        pointing = "lmp{:03d}".format(float(o.split()[1]))
    
    # read in filter
    ffilter = header['FILTER']

    # check to see if this is a weight image
    if f.find('weight') > -1:
        suffix='coadd.weight.fits'
    else:
        suffix='coadd.fits'
    outfile = "VF-{}-{}-{}-{}-{}".format(dateobs,instrument,pointing,ffilter,suffix)
    print('moving ',f,' -> ',outfile)
    #os.rename(f,outfile)
