#!/usr/bin/env python

'''
  GOAL:
  make lists of files that contain
  - dome flats with same filter
  - sky flats with same filter
  - bias frames
  - dark frames with the same exposure time


  mv c7133t014[3-9]f00.fits junk/.
'''
import glob
import os
import numpy as np

filters=['V','R','ha8','ha12','ha16']
# group files using gethead
# OBSTYPE (object,flat,bias)
# CMMTOBS - filter
# bias frames
os.system('ls *b00.fits > biasframes')

object_files=glob.glob('*o00.fits')
flat_files=glob.glob('*f00.fits')
all_files=glob.glob('*.fits')
os.system('gethead *f00.fits CMMTOBS > junkfile')
infile=open('junkfile','r')
fnames=[]
ftype=[]
for line in infile:
    t=line.split('.fits')
    fnames.append(t[0]+'.fits')
    ftype.append(t[1].rstrip('\n').replace(' ',''))
infile.close()

ftype=np.array(ftype) # make into character array
flattypes=['domeflat','skyflat']
for f in filters:
    for t in flattypes:
        flatgroup=t+f
        print 'flatgroup = ',flatgroup
        indices=np.where(ftype == flatgroup)
        print indices, len(indices[0])
        if len(indices[0]) > 1:
            outfile = open(flatgroup,'w')
            for i in indices[0]:
                print fnames[i],ftype[i]
                outfile.write(fnames[i]+'\n')
            outfile.close()
            
                

