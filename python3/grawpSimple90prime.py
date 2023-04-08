#!/usr/bin/env python

'''
GOAL:
* This will create a bash script to run multiple serial jobs in parallel using slurms array option

USAGE:
* when logging into grawp
  module load Python3

* Run this in the directory that has all of the ksb files from the NOAO pipeline
* create a list of input images

rfinn@localhost:/mnt/qnap_home/rfinn/Halpha/BokNOAOPipeline$ ls ksb*ooi*v1.fits > all_ksb_images

* Try running this once and check the output of the JOB_{}.sh script.
* If it looks good, run again using the --submit flag to submit the script to slurm.

'''

import os
import subprocess
import sys
import glob
import time
HOME = os.getenv("HOME")


program= f"{HOME}/github/HalphaImaging/python3/BOK_pipeline_fixampoffsets.py"
input_file = "all_ksb_images"
infile = open(input_file,'r')
i = 0
for line in infile:
    image = line.rstrip()
    outimage = 'zm'+image
    if not os.path.exists(outimage):
        cmds = ['python3', program,image]
        print(f"Submitting job to process {input_file}")
        process = subprocess.Popen(cmds)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #stdout, stderr = process.communicate()
        #print(stdout.decode())
        #print(stderr.decode())
        i += 1
        if i > 20:
            # submit 20 at a time
            break
        #if i%20 == 0:
        #    time.sleep(60)

