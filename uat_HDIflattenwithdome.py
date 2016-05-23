#!/usr/bin/env python

'''
  GOAL:
   - Takes all object images and flattens them using a normalized dome flat in the same filter

EXAMPLE:
   In the directory containing all objects to flatten type in the command line:
      '/home/share/research/pythonCode/uat_HDIflattenwithdome.py'(or whatever the path is to where this program is stored)
  
'''
import glob
import os
import numpy as np
from pyraf import iraf


os.system('gethead tr*o00.fits CMMTOBS > junkfile1')
infile=open('junkfile1','r')
fnames=[]
ffilter=[]
for line in infile:
    t=line.split('.fits')
    fnames.append(t[0]+'.fits')
    ffilter.append(t[1].rstrip('\n').replace(' ',''))
infile.close()
filters=set(ffilter)

ffilter=np.array(ffilter) # make into character array
for f in filters:
    objectgroup='object_'+f
    fobjectgroup='fobject_'+f
    print 'objectgroup = ',objectgroup
    indices=np.where(ffilter == f)
    if len(indices[0]) > 0:
        outfile = open(objectgroup,'w')
        outfile2 = open(fobjectgroup,'w')
        for i in indices[0]:
            outfile.write(fnames[i]+'\n')
            outfile2.write('d'+fnames[i]+'\n')
        outfile.close()
        outfile2.close()
    iraf.imarith(operand1 = "@"+objectgroup, op = "/", operand2="ndomeflat"+f, result = "@" + fobjectgroup)
os.remove('junkfile1')


