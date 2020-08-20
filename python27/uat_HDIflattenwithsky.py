#!/usr/bin/env python

'''
  BASIC INFORMATION ABOUT THIS CODE:
  -This is the THIRD program you need to run in order to perform the data recuti  on process of HDI images.
  -This code will flatten the combined and normalized sky flat images. This is a  very important part of the data reduction process because we need to take into  account of the multiplicative effect of temperature fluctutations within the C  CD (which sky f  lats can help correct for).
  -Sky flats correct for pixel to pixel variations within the CCD. Sky flats are   a type of flatfielding correction in the data reduction process. 

  BEFORE RUNNING THIS CODE:
  -Be sure to type ur_setup in the terminal everytime you open a new terminal wi  ndow so that ureka is activated.
  -Ensure that pyraf is still activated by retyping the commands listed in the c  omments of the FIRST program titled "uat_HDIgroupflatfiles.py".

  PROCEDURE:
  -The FIRST program grouped the flats and corrected for additive effects within  the CCD, but now we have to use the flatfielding technique (just like we did i  n Program 2 with the dome flats) to help correct for the multiplicative effect  of pixel to pixel variations.
  -The specific task that is performed to flatfield the images is the task "imar  ith" which is the line of code that actually divides the images by the normali  zed pixel values which will correct for the multiplicative effect.  

  GOAL:
  - Takes all object images and flattens them using a normalized sky flat in th   e same filter
  - It does this through division because these flats are correcting for the mul  tiplicative effect of pixel to pixel variations. 

  EXAMPLE:
   In the directory containing all objects to flatten type in the command line:
      '/home/share/research/pythonCode/uat_HDIflattenwithsky.py'(or whatever the path is to where this program is stored)

  WHAT THIS CODE DOES:
  -Mainly, this code uses the created organized folders containing normalized fl  at files and grabs them in a for loop and everytime it loops over it combines   the flats and performs the flatfield division thereafter. 

  INPUT/OUTPUT:
  Input --> tr*o00.fits (trimmed images within the created flat folders from pro                        gram 1 titled "uat_HDIgroupflatfiles.py")
  Output --> ftr*o00.fits (f stands for flattened after it gets flattened with t                          hese sky flats). 

  REQUIRED MODULES:
  -pyraf

  EXTRA NOTES:
  -We ended up only flattening our HDI images using sky flats rather than our do  me flats because they were not sufficient enough. 

  WRITTEN BY:
  Dr. Rose Finn
  EDITED BY: 
  Research Team Summer 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood,   Kaitlyn Hoag, and Kelly Whalen.  
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
            outfile2.write('f'+fnames[i]+'\n')
        outfile.close()
        outfile2.close()
    iraf.imarith(operand1 = "@"+objectgroup, op = "/", operand2="nskyflat"+f, result = "@" + fobjectgroup)
os.remove('junkfile1')


