#!/usr/bin/env python

'''
  GOAL:
  make lists of files that contain
  - the same object
  - the same filter
  - are exposed for over a minute


  REQUIREMENTS:
  - only looks at images with the title in the format f*o00.fits
  - must be run in folder containing files to sort
  - will create the output file in the format object_filter

  EXAMPLE:
   In the directory containing all flattened objects with fixed headers to be sorted type in the command line:
      '/home/share/research/pythonCode/uat_HDIsortobjects.py'(or whatever the path is to where this program is stored)

written by Rose A. Finn
edited by Grant Boughton & Kelly Whalen

   
'''
import glob
import os
import numpy as np
import argparse
from astropy.io import fits

parser = argparse.ArgumentParser(description ='group objects by filter and target for combining with swarp')
parser.add_argument('--filestring', dest = 'filestring', default = 'hcftr', help = 'string to use to get input files (default = "hcftr" which grabs all files "hcftr*o00.fits")')

args = parser.parse_args()


# group files using gethead
# CMMTOBS - filter
# OBJECT - title of observed object
# EXPTIME - exposure time

filestring = args.filestring+'*o00.fits'
try:
    os.system('gethead '+filestring+' FILTER, OBJECT, EXPTIME > junkfile2')
    infile=open('junkfile2','r')
    fnames=[]      #creates empty list to contain file name
    ftype=[]       #creates empty list to contain type of filter
    fobject=[]     #creates empty list to contain name of object
    exptime=[]     #creates empty list to contain time(s) exposed
    for line in infile:
        t=line.split()      #splits input at each space and puts each element into                          list t
        fnames.append(t[0]) #adds the element in index 0 of list t to the empty                             list fnames
        fobject.append(t[2].rstrip('\n').replace(' ',''))
        exptime.append(t[3].rstrip('\n').replace(' ',''))
        ftype.append(t[1])
    
    infile.close()

except:
    files = glob.glob(filestring)

    fnames=[]      #creates empty list to contain file name
    ftype=[]       #creates empty list to contain type of filter
    fobject=[]     #creates empty list to contain name of object
    exptime=[]     #creates empty list to contain time(s) exposed

    for f in files:
        data, header = fits.getdata(f,header=True)
        fnames.append(f)
        ftype.append(header['FILTER'])
        fobject.append(header['OBJECT'])
        exptime.append(header['EXPTIME'])

filters=set(ftype)    #takes all inputs from ftype and adds the unique ones to                        the set of filters
objecttypes=set(fobject)
  
ftype=np.array(ftype) # make into character array
fobject=np.array(fobject)
exptime=np.array(exptime,'f')

for f in filters:
    for t in objecttypes:
        objectgroup=str(t)+ '_' + str(f) #creates a variable that is the                                                 combined characteristics of the image
        indices=np.where((f == ftype) & (t == fobject) & (exptime > 60)) #finds                                                                       indices where the filter and name                                                            of object match the values searched                                                          for in the for loop
        if len(indices[0]) > 0:
            if objectgroup.find(' ') > -1: # get rid of long filter names
                t = objectgroup.split()
                filename = t[0]+'-'+t[1]
            else:
                filename = objectgroup
            outfile = open(filename,'w') #writes the name of file to a
                                            #file titled the characteristics
                                            #stated in the variable objectgroup
            for i in indices[0]:
                outfile.write(fnames[i]+'\n')
            outfile.close()
            
                
os.remove('junkfile2')



