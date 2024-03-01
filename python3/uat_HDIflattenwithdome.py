#!/usr/bin/env python

'''
GOAL:
- Takes all object images and flattens them using a normalized dome flat in the same filter

PROCEDURE:
- get all files that start with ztr*.fits
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
  Rose Finn
  
EDITED BY:
  Research Team Summer 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood, Kaitlyn Hoag, and Kelly Whalen
  Finn May 2018 to use python and not pyraf

  updated 7/27/18 by RAF
  
'''
import glob
import os
import numpy as np
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(description ='Flatten images with dome flat taken in corresponding filter')
parser.add_argument('--filestring', dest='filestring', default='ztr', help='match string for input files (default =  ztr, which gets ztr*.fits)')
parser.add_argument('--verbose', dest='verbose', default=False,action='store_true', help='print extra messages for troubleshooting')
args = parser.parse_args()



os.system('gethead '+args.filestring+'*o00.fits CMMTOBS > junkfile1')
infile=open('junkfile1','r')
filenames=[]
filefilter=[]
for line in infile:
    t=line.split('.fits')
    filenames.append(t[0]+'.fits')
    if (line.find('6620') > -1) | (line.find('ha4') > -1) | (line.find('Ha4') > -1) :
        filefilter.append('ha4')
    elif (line.find('6660') > -1) |(line.find('ha8') > -1) | (line.find('Ha4') > -1) :
        filefilter.append('ha8')
    elif (line.find('6700') > -1) |(line.find('ha12') > -1) | (line.find('Ha4') > -1) :
        filefilter.append('ha12')
    elif (line.find('6740') > -1) |(line.find('ha16') > -1) | (line.find('Ha4') > -1) :
        filefilter.append('ha16')
    
    elif t[1].find('R') > -1:
        filefilter.append('R')
    elif t[1].find('r') > -1:
        filefilter.append('r')
    elif t[1].find(' r ') > -1:
        filefilter.append('r')
    elif t[1].find('V') > -1:
        filefilter.append('V')
    else:
        print("WARNING: no filter found for ",t[1])
    #filefilter.append(t[1].rstrip('\n').replace(' ',''))
infile.close()
filters=set(filefilter)
print('these are the filters I got: ',filters)

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
    print('PROCESSING IMAGES IN FILTER ',f)
    objectgroup='object_'+f
    fobjectgroup='fobject_'+f
    print('objectgroup = ',objectgroup)
    if f == 'ha4':
        indices=np.where((filefilter == f) | (filefilter == '6620'))
    else:
        indices=np.where(filefilter == f)
    if len(indices[0]) > 0:
        # open flat file
        if not(os.path.exists("ndomeflat"+f+".fits")):
            print("can't find "+"ndomeflat"+f)
            if f == 'ha4':
                print('looking for ndomeflat6620.fits')
                f = 'ndomeflat6620.fits'
                print(f)
                print(os.path.exists(f))
                if os.path.exists(f):
                    flatfile = f
                else:
                    print("skipping to next filter")
                    continue
            else:   
                print("skipping to next filter")
                continue
        else:
            flatfile = "ndomeflat"+f+".fits"
        flatdata = fits.getdata(flatfile)
        for i in indices[0]:
            data,header = fits.getdata(filenames[i],header=True)
            dataout = data / flatdata
            header['HISTORY'] = 'Flattened using ndomeflat'+f
            fits.writeto('d'+filenames[i],dataout,header,overwrite=True)
os.remove('junkfile1')


