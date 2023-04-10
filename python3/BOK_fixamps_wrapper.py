#!/usr/bin/env python

import os
import subprocess
import sys
import glob
import time
HOME = os.getenv("HOME")


program= f"{HOME}/github/HalphaImaging/python3/BOK_pipeline_fixampoffsets.py"
filelist = glob.glob('ksb*ooi_r_v1.fits')
filelist.sort()
print(f'found {len(filelist)} files to process...')
for i in range(len(filelist)):
    print()
    print(f"processing {i}/{len(filelist)}")
    
    
    image = filelist[i]
    outimage = 'zm'+image
    if os.path.exists(outimage):
        print(f'{outimage} found - moving to next object')
    else:
        print(f"processing {image} -> {outimage}")
        cmd = f"python {program} {image} r"
        os.system(cmd)
        #cmds = ['python3', program,image]
        #print(f"Submitting job to process {input_file}")
        #process = subprocess.Popen(cmds)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #stdout, stderr = process.communicate()
        #print(stdout.decode())
        #print(stderr.decode())
