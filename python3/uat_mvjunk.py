#!/usr/bin/env python

'''
GOAL
- make a subdirectory called junk (if it doesn't already exist)
- move unwanted files to junk

PROCEDURE
- user creates a file called junkfiles with filenumbers to move
- to create file, use 
- example:
f 0060-0068 
o 0074 
b 0100-0153

- first letter signifies flat, object, bias
- followed by file number

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
    t = line.split()
    prefix = t[0]
    number = t[1]
    if number.find('-') > -1:
        t = number.rstrip().split('-')
        print(t)
        a = int(t[0])
        b = int(t[1])+1
        for i in range(a,b):
            print('mv *%04i*%s00.fits junk/.'%(i,prefix))
            os.system('mv *%04i*%s00.fits junk/.'%(i,prefix))
                       
    else:
        print(f'mv *{number.rstrip():04i}*{prefix}00.fits junk/.')
        os.system(f'mv *{number.rstrip():04i}*{prefix}00.fits junk/.')
        pass

infile.close()
