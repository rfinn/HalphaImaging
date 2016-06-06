#!/usr/bin/env python

'''
  BASIC INFORMATION ABOUT THIS CODE:
  -This is the FIRST program you need to run in order to perform the data reduct  ion of HDI images.
  -This code will create trimmed files of your original images which will remove  overscan regions in addition to combining the flats by grouping together the f  ilter and flat type and normalizing those combined images.
 
  BEFORE RUNNING THIS CODE:
  -Type "ur_setup" into the terminal to enable ureka
  -Make sure that you are able to run these codes even if your images are in a s  eperate directory. In order to do so, you must type the following into terminal 
          emacs -nw .profile
          export PATH=$PATH:~/directory_of_your_python_code
      -Check to see if this added path was successfully added by typing: 
          echo $PATH
  -Create a junk director to put pointing tests, focus frames, and bad exposures  into within the directory where all your images are stored so that they are el  iminated from the data reduction process.
  -IN ORDER TO START PYRAF --> type the following into terminal
          pyraf
          imred
          ccdred
          epar ccdproc --> use parameters that are listed on our authorea websit          e titled "Halpha Imaging of Nearby Galaxy Groups" written by Dr. Rose           Finn, Grant Boughton, Natasha Collova, Sandy Spicer, and Kelly Whalen. 

  GOAL:
  Make lists of files that contain
  - dome flats with same filter
  - sky flats with same filter
  Then combine and normalize the flats so they can be used to flatten image

  PROCEDURE:
  -User should move junk files to a junk subdirectory before starting
   - Junk files include initial bias frames pointing and other garbage frames
  - Use gethead to pull relevant header data
  - Overscan subtract and trim images (we assume these image names begin with 't  r')
  - Combine flats according to flat type (dome vs sky) and filter

  EXAMPLE:
     In the directory containing all flats type in the command line:
     '/home/share/research/pythonCode/uat_HDIgroupflatfiles.py'(or whatever the      path is to where this program is stored)

  WHAT THIS CODE DOES:
  -This program combines and stores the information of the normalized flats of t  hat same filter so that this additive effect could be subtracted from our data  This code also subtracts the overscan regions.
  -This program does the above be grouping all the flat files by filter and flat   type so that they can be combined accordingly and normalized thereafter.  
  
  INPUT/OUTPUT:
  Input: 'skyflat(type of filter)' or 'domeflat(type of filter)'
         -These contain trimmed image files grouped by filter grabbed from the i         mage headers seen at the beginning of this code. 
  Output: For combined flats --> 'cskyflat(type of filter).fits' or 'cdomeflat(t          ype of filter).fits'.
          For normalized flats --> 'nskyflat(type of filter).fits' or 'ndomeflat          (type of filter).fits'.

  REQUIRED MODULES:
  pyraf

  NOTES:
  in junkfile ftr flats still show. We changed the gethead requirements to only   bring in files that start with tr but the ftr files will not go away! =(
  
  WRITTEN BY:
  Rose A. Finn
  EDITED BY:
  Natasha Collova, Tiffany Flood, and Kaitlyn Hoag 5/29/15
  Grant Boughton, Natasha Collova, and Sandy Spicer 6/3/16
  UPDATES:
  Now combines and normalizes the grouped flat files 
    Input ex. = skyflatR
    Output ex. = cskyflatR.fits (combined sky flats in R band) & nskyflatR.fits     (normalized sky flats in R band)
  
'''
import glob
import os
import numpy as np
from pyraf import iraf

iraf.noao()
iraf.imred()
iraf.ccdred()    

os.system('gethead tr*f00.fits CMMTOBS > tempflats')    #tempflats is the name of a "junk file" that contains the gethead information from all the flat images. This file will be deleted after the information is read out in the future. We are left with multiple arrays of information from the headers.
#we assume that the flat images are trimmed and the file name starts with 'tr'
infile=open('tempflats','r')
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
        print 'grouping files for filter type = ',element
        indices=np.where((array_ftype == f)&(array_filter == element))
        if len(indices[0]) > 0:
            outfile = open(ftype_filter,'w')
            for i in indices[0]:
                outfile.write(fnames[i]+'\n')
            outfile.close()
os.remove('tempflats')            
flats = glob.glob('*flat*')
flats=set(flats)-set(glob.glob('c*.fits'))-set(glob.glob('n*.fits'))     # doesn't include flat files that have already been combined or normalized in the following loop
for f in flats:
    iraf.imcombine(input='@'+f, output='c'+f, combine='median', scale='median') # combine flat, create e.g. cdomeflatR
    flat_mean = iraf.imstat(images='c'+f, fields='mean', lower='INDEF', format=0, Stdout=1) # find mean of combined flat image
    iraf.imarith(operand1='c'+f, op='/', operand2=flat_mean[0], result='n'+f) # normalize the combined flat
    

                

