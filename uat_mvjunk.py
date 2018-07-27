#!/usr/bin/env python

'''
GOAL
- make a subdirectory called junk (if it doesn't already exist)
- move unwanted files to junk

PROCEDURE
- user creates a file called junkfiles with filenumbers to move
- to create file, use 
- example:
0060-0068
0074
0100-0153

USAGE:
from within ipython

%run ~/github/HalphaImaging/uat_mvjunk.py --filename junkfiles

'''


import os
import argparse

parser = argparse.ArgumentParser(description ='Reads in images numbers and moves images to a subdirectory called junk.  File should contain 4 digit numbers that uniquely identify the image.  Ranges can be specified on one line using a dash.  \n For example: \n\n 0001 \n 0034-0042 \n')
parser.add_argument('--filename', dest = 'filename', default = 'junkfiles', help = 'file containing list of 4-digit image numbers to move to junk.  Default filename is junkfiles.')

args = parser.parse_args()

if not(os.path.exists('junk')):
    os.mkdir('junk')

infile = open(args.filename)
for line in infile:
    if line.find('-') > -1:
        t = line.rstrip().split('-')
        #print t
        a = int(t[0])
        b = int(t[1])+1
        for i in range(a,b):
            os.system('mv *%04i*.fits junk/.'%(i))
                       
    else:
        os.system('mv *'+line.rstrip()+'*.fits junk/.')

infile.close()
