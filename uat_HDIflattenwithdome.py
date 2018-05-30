#!/usr/bin/env python

'''
GOAL:
- Takes all object images and flattens them using a normalized dome flat in the same filter

PROCEDURE:
- get all files that start with tr*.fits
- get filter for each image
- create a list of unique filters
- divide each image by the corresponding flat
- assumes that the flats already exist and are called ndomeflatR (normalized domeflat in R filter)

EXAMPLE:
- In the directory containing all objects to flatten type in the command line:

/home/share/research/pythonCode/uat_HDIflattenwithdome.py


INPUT/OUTPUT:
- Input --> tr*o00.fits (Trimmed fits files)
- Output --> dtr*o00.fits (Dome flattened and trimmed fits files)

REQUIRED MODULES:


EXTRA NOTES:
  -**We did not reduce our HDI data from KPNO in April 2015 using our dome flats because they were insufficient. We used sky flats instead.
  

WRITTEN BY:
  Dr. Rose Finn
  
EDITED BY:
  Research Team Summer 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood, Kaitlyn Hoag, and Kelly Whalen
  Finn May 2018 to use python and not pyraf
'''
import glob
import os
import numpy as np
from astropy.io import fits


os.system('gethead tr*o00.fits CMMTOBS > junkfile1')
infile=open('junkfile1','r')
filenames=[]
filefilter=[]
for line in infile:
    t=line.split('.fits')
    filenames.append(t[0]+'.fits')
    filefilter.append(t[1].rstrip('\n').replace(' ',''))
infile.close()
filters=set(filefilter)
print 'these are the filters I got: ',filters

filefilter=np.array(filefilter) # make into character array

## for f in filters:
##     objectgroup='object_'+f
##     fobjectgroup='fobject_'+f
##     print 'objectgroup = ',objectgroup
##     indices=np.where(filefilter == f)
##     if len(indices[0]) > 0:
##         outfile = open(objectgroup,'w')
##         outfile2 = open(fobjectgroup,'w')
##         for i in indices[0]:
##             outfile.write(filenames[i]+'\n')
##             outfile2.write('d'+fnames[i]+'\n')
##         outfile.close()
##         outfile2.close()
##     iraf.imarith(operand1 = "@"+objectgroup, op = "/", operand2="ndomeflat"+f, result = "@" + fobjectgroup)
## os.remove('junkfile1')

# rewriting to use python and not pyraf


ffilter=np.array(filefilter) # make into character array
for f in filters:
    print 'PROCESSING IMAGES IN FILTER ',f
    objectgroup='object_'+f
    fobjectgroup='fobject_'+f
    print 'objectgroup = ',objectgroup
    indices=np.where(filefilter == f)
    if len(indices[0]) > 0:
        # open flat file
        if not(os.path.exists("ndomeflat"+f+".fits")):
            print "can't find "+"ndomeflat"+f
            print "skipping to next filter"
            continue
        else:
            flatdata = fits.getdata("ndomeflat"+f+".fits")
            for i in indices[0]:
                data,header = fits.getdata(filenames[i],header=True)
                dataout = data / flatdata
                header['HISTORY'] = 'Flattened using ndomeflat'+f
                fits.writeto('d'+filenames[i],dataout,header,clobber=True)
os.remove('junkfile1')


