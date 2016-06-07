#!/usr/bin/env python

'''
  BASIC INFORMATION ABOUT THIS CODE:
  -This is the SECOND program you need to run in order to perform the data reduc  tion process of HDI images.
  -This code will flatten the combined and normalized dome flat images. This is   an important part of the data reduction process because we need to take into a  ccount of the effect of temperature fluctuations within the CCD (which dome fl  ats will correct for).
  
  BEFORE RUNNING THIS CODE:
  -Be sure to type ur_setup in the terminal everytime you open a new terminal wi  ndow so that ure  ka is activated.
  -Ensure that pyraf is activated still by retyping the commands listed in the c  omments of the FIRST program titled "uat_HDIgroupflatfiles.py"

  PROCEDURE:
  -The FIRST code only grouped the flat files and organized them and combined th  em, but only did the math behind normalizing the images. In this way, we had t  o grab all the grouped together files of flats from the previous code and flat  ten them by subtracting out the additive noise from the CCD.
    -More specifically, after we performed the said above using a for loop by ac    cessing the newly created files from the FIRST code, and performing the flat    tening by dividing the normalized flat by the original image so that the exc    ess noise is removed. 

  GOAL:
      - Takes all object images and flattens them using a normalized dome flat i      n the same filter
      - This program is necessary to perform the flattening of the fits files us      ing the dome flats. We need to perform this step so that the noise within       the CCD can be removed using dome flats which correct for temperature fluc      tuations within the CCD. 

  EXAMPLE:
      In the directory containing all objects to flatten type in the command lin      e:
        '/home/share/research/pythonCode/uat_HDIflattenwithdome.py'(or whatever         the path is to where this program is stored)

  WHAT THIS CODE DOES:
  -Mainly, this code performs the actual division of the dome flats so that the   noise from the CCD can be taken out. 

  INPUT/OUTPUT:
  Input --> tr*o00.fits (Trimmed fits files)
  Output --> dtr*o00.fits (Dome flattened and trimmed fits files)

  REQUIRED MODULES:
  -pyraf
  
  EXTRA NOTES:
  -**We did not reduce our HDI data from KPNO in April 2015 using our dome flats  because they were insufficient. We used sky flats instead, only to account fo   r the additive effect of temperature within the CCD with increasing exposure    time. The sky flats were more sufficient in testing linearity and seeing wher   e the counts leveled off and were used more so throughout the observing proce   ss.  

  WRITTEN BY:
  Dr. Rose Finn
  EDITED BY:
  Research Team Summer 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood,   Kaitlyn Hoag, and Kelly Whalen
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


