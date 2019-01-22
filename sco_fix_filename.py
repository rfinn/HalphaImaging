#!/usr/bin/env python

import glob
import os


files = glob.glob('*.fit')
for f in files:
    output_name = f.split('.fit')[0]+'.fits'
    os.rename(f,output_name)

