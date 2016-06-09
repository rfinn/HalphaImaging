#!/usr/bin/env python
import numpy as np
import argparse

parser = argparse.ArgumentParser(description ='convert output from mscgetcatalog to ds9 regions file')
parser.add_argument('--input', dest='input', default=None, help='output from mscgetcatalog')
parser.add_argument('--output', dest='output', default=None, help='ds9 region file (should end with .reg)')
parser.add_argument('--radius', dest='radius', default='10.', help='radius of circle in arcsec')
args = parser.parse_args()

c=np.loadtxt(args.input,dtype={'names':('ra','dec','u','g','r'),'formats':('S12','S12','f','f','f')},unpack=True)

outfile=open(args.output,'w')
outfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile.write('fk5\n')

for i in range(len(c[0])):
    s= 'circle('+c[0][i]+','+c[1][i]+','+str(args.radius)+'") \n'
    outfile.write(s)
outfile.close()
