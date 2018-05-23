#!/usr/bin/env python

'''

BASIC INFORMATION ABOUT THIS CODE:
  -This is the FOURTH program you need to perform the data reduction process of   HDI images.
  -This code gets rid of short exposure times (less than 61 seconds) so that the  only exposure frames we will use to analyze our images are exposures that are   more clear and less likely to have abnormalities in the image itself. 

BEFORE RUNNING THIS CODE:
  -Be sure to type ur_setup in the terminal everytime you open a new terminal window so that ureka is activated.

PROCEDURE:
  -This code creates a subdirectory called short_exptime where all short exptimes will go to so they are not included in the data reduction process.
  -This code uses args.parse to also parse the code so that it was able to call on the directory when the for loop found images shorter than the limited exposure time (< 61 seconds).
  -Usually when observing, the majority of object images taken with such short o  f exposure times are used for pointing stars, focus frames, standard stars, etc. Those of which we do not want to include in our reduction. 

GOAL:
  -Mainly, the goal of this code is to elminate short exposure times. It does th   is so that we only have necessary exposure times for the data reduction. Short exposure times are usually not object frames and/or object images that we d   on't want in our data. There was a reason that whoever was observing didn't take a longer exposure time on that object, perhaps a focusing issue, pointing issue, etc. 

EXAMPLE:
  -In the directory, containing all fits files that are less than the 61 second exposure time limit, type in the ocmmand line:
   '/home/share/research/pythonCode/uat_mvshortexptime.py'

WHAT THIS CODE DOES:
  -Mainly, this code moves the short exposure times to the subdirectory titled "   short_exptime". 

INPUT/OUTPUT:
  Input --> all ftr*.fits in directory
  Output --> None, just moves those ftr*.fits files that have a small exposure t             ime into a subdirectory so we don't use these images.

REQUIRED MODULES:
  -pyraf

EXTRA NOTES:
  -We wrote this code during the summer of 2015 and wrote it in this way to prac  tice using args.parse in a program. Of course, there are many different ways o  f writing a program like this.

WRITTEN BY:
Dr. Rose Finn
EDITED BY:
Research Team Summer 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood, Kaitlyn Hoag, and Kelly Whalen  



'''

import argparse
import glob
import os
import subprocess

parser = argparse.ArgumentParser(description ='Move images with short exposure times to a subdirectory')
parser.add_argument('--filestring', dest='filestring', default='ftr', help='match string for input files (default =  ftr, which looks for ftr*.fits files)')
parser.add_argument('--minexptime', dest='minexptime', default=61., help='min exposure time of science frames (default = 61 sec)')
parser.add_argument('--subdir', dest='subdir', default='short_exptime', help='subdirectory for short exptime images (default = short_exptime)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)
i=1
try:
    os.mkdir(args.subdir)
except OSError:
    print args.subdir,' already exists'
print 'subdirectory for short exposure time images is ',args.subdir
for f in files:
    print 'CHECKING EXPTIME FOR FILE %i OF %i'%(i,nfiles)
    read_exptime = 'gethead ' + f + ' EXPTIME'
    exptime = subprocess.check_output(read_exptime,shell=True)
    exptime = exptime.rstrip()
    if float(exptime) < args.minexptime:
        os.rename(f,args.subdir+'/'+f)
    i += 1
    print '\n'

    
    


                   
