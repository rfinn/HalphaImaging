#!/usr/bin/env python

'''
writing this to change filter from 

e.g. pointing_24 r  

to 

r

This is an issue for Feb 2020 data.  
not sure what we did wrong when taking data.  
As Greg suggested, should just get filter from 
filter wheel info, rather than relying on 
observer to enter OBJECT and COMMENT correctly.

'''
import argparse
import glob
from astropy.io import fits

parser = argparse.ArgumentParser(description ='Edit image headers to correct filter name')
parser.add_argument('--filestring', dest='filestring', default='hdz', help='match string for input files (default =  hdz, which operates on all hdz*.fits files)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)
i=1
for f in files:
    print('FIXING HEADER FOR FILE %i OF %i'%(i,nfiles))
    data, header = fits.getdata(f,header=True)
    t = header['FILTER']
    o = header['OBJECT']
    new_filter = t.split()[-1]
    #print('\t older filter = ',t)
    print(f,'\t',o,'\t',new_filter)
    header['FILTER']=new_filter
    header['OBJECT'] = o.replace('-','_')


    print('WRITING UPDATED FILE')
    fits.writeto(f,data,header,overwrite=True)
    i += 1
    print('\n')

