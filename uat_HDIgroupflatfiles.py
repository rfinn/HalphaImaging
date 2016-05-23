#!/usr/bin/env python

'''
  GOAL:
  Make lists of files that contain
  - dome flats with same filter
  - sky flats with same filter
  Then combine and normalize the flats so they can be used to flatten image

  PROCEDURE:
  -User should move junk files to a junk subdirectory before starting
   - Junk files include initial bias frames pointing and other garbage frames
  - Use gethead to pull relevant header data
  - Overscan subtract and trim images (we assume these image names begin with 'tr'
  - Combine flats according to flat type (dome vs sky) and filter

  EXAMPLE:
     In the directory containing all flats type in the command line:
      '/home/share/research/pythonCode/uat_HDIgroupflatfiles.py'(or whatever the path is to where this program is stored)

  
  REQUIRED MODULES:
  pyraf

  NOTES:
  in junkfile ftr flats still show. We changed the gethead requirements to only bring in files that start with tr but the ftr files will not go away! =(
  
  WRITTEN BY:
  Rose A. Finn
  EDITED BY:
  Natasha Collova, Tiffany Flood, and Kaitlyn Hoag 5/29/15

  UPDATES:
  
'''
import glob
import os
import numpy as np
from pyraf import iraf
    
os.system('gethead tr*f00.fits CMMTOBS > junkfile')
#we assume that the flat images are trimmed and the file name starts with 'tr'
infile=open('junkfile','r')
fnames=[]
filter=[]
ftype=[]   #skyflat or domeflat
for line in infile:
    t=line.split()
    fnames.append(t[0])
    ftype.append(t[1]+t[2])
    filter.append(t[3].rstrip('\n'))
infile.close()
set_filter=set(filter)
set_ftype=set(ftype)
array_ftype=np.array(ftype)
array_filter=np.array(filter)

for f in set_ftype:
    print "flat type=",f
    for element in set_filter:
        ftype_filter = str(f)+str(element)
        print 'filter type = ',element
        indices=np.where((array_ftype == f)&(array_filter == element))
        print indices, len(indices[0])
        if len(indices[0]) > 0:
            outfile = open(ftype_filter,'w')
            for i in indices[0]:
                outfile.write(fnames[i]+'\n')
            outfile.close()
os.remove('junkfile')            
flats = glob.glob('*flat*')
for f in flats:
    iraf.flatcombine(input='@'+f, output='c'+f, combine=median)
    mean = iraf.imstat(input=f, fields='mean', lower='INDEF', format=0, Stdout=1)
    iraf.imarith(operand1='@c'+f, op='/', operand2=mean, result='n'+f

                

