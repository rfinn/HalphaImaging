#!/usr/bin/env python
'''
GOAL:
- generate cutouts for CO sample
- call plot_cutouts_ha for each galaxy

NOTES:
- run for cutouts folder
- program will make a directory of the prefix name, if it doesn't already exist
- program then moves to prefix directory and gets cutouts
- move back up to cutouts directory and repeat

'''
import os
homedir = os.getenv("HOME")
os.sys.path.append(homedir+'/github/HalphaImaging/python3/')
import plot_cutouts_ha as pc
from astropy.table import Table
vftabledir = homedir+'/research/Virgo/tables-north/v0/'

# read in vf main file
vfmain = Table.read(vftabledir+'vf_north_v0_main.fits')

# select CO galaxies
vfmain_co = vfmain[vfmain['COflag']]

# run plotcutouts_ha for each galaxy
#for i in range(2):    # for testing
for i in range(len(vfmain_co)):
    if not os.path.exists(vfmain_co['prefix'][i]):
        os.mkdir(vfmain_co['prefix'][i])
    os.chdir(vfmain_co['prefix'][i])
    if vfmain_co['radius_flag'][i]:
        rad = 6*vfmain_co['radius'][i]
    else:
        rad = 100
    c = pc.cutouts('',ra=vfmain_co['RA'][i],dec=vfmain_co['DEC'][i],size=rad,galid=vfmain_co['prefix'][i])
    c.runall()
    c.plotallcutouts()
    c.plotsfrcutouts()    
    os.chdir('..')
