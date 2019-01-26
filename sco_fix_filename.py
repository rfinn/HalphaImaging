#!/usr/bin/env python

import glob
import os


files = glob.glob('*.fit')
if len(files) > 0:
    for f in files:
        output_name = f.split('.fit')[0]+'.fits'
        os.rename(f,output_name)


subdirs = ['DARKS','FLATS','BIAS']
prefixes = ['dark','flat','bias','zero']
filedict = {'dark':'DARKS','flat':'FLATS','bias':'BIAS','zero':'BIAS'}
for s in subdirs:
    if os.path.exists(s):
        print(s+' directory exists')
    else:
        os.mkdir(s)

files = glob.glob('*.fits')
for f in files:
    for p in prefixes:
        if f.find(p) > -1:
            os.rename(f,filedict[p]+f)
        
