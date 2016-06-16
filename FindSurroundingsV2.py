#!/usr/bin/python
# coding: utf-8


import numpy as np
from astropy.io import fits
import time
import argparse

parser = argparse.ArgumentParser(description ='Creates line matched tables for area surrounding input RA and DEC')
parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Use shortened version of NSA catalog (nsa_uat.fits)')
parser.add_argument('--RA', dest = 'r', default = 0, help = 'Target RA')
parser.add_argument('--DEC', dest = 'd', default = 0, help = 'Target DEC')
parser.add_argument('--radius', dest = 'rd', default = 1, help = 'Radius around target to include in tables')
parser.add_argument('--prefix', dest = 'pre', default = 'test', help = 'Prefix of new file names')
args = parser.parse_args()


def findSurroundings2(ra, dec, radius, prefix):
    start_time = time.time()
    if arg.s:
        petro = fits.getdata('YangDR7PetroToNSA_uat.fits')
        model = fits.getdata('YangDR7ModelToNSA_uat.fits')
        nsadat = fits.getdata('nsa_uat.fits')
        simard1 = fits.getdata('Simard1ToNSA_uat.fits')
        simard2 = fits.getdata('Simard2ToNSA_uat.fits')
        simard3 = fits.getdata('Simard3ToNSA_uat.fits')
    else:
        petro = fits.getdata('YangDR7PetroToNSA.fits')
        model = fits.getdata('YangDR7ModelToNSA.fits')
        nsadat = fits.getdata('nsa_v0_1_2.fits')
        simard1 = fits.getdata('Simard1ToNSA.fits')
        simard2 = fits.getdata('Simard2ToNSA.fits')
        simard3 = fits.getdata('Simard3ToNSA.fits')
    d=np.sqrt((ra-nsadat.RA)**2 + (dec-nsadat.DEC)**2)
    index=np.arange(len(d))
    matches=index[d<=radius]
    print '%i entries found' % len(matches)
    if arg.s:
        fits.writeto(prefix+'_NSA_uat.fits',nsadat[matches],clobber=True)
        fits.writeto(prefix+'_Petro_uat.fits',petro[matches],clobber=True)
        fits.writeto(prefix+'_Model_uat.fits',model[matches],clobber=True)
        fits.writeto(prefix+'_Simard1_uat.fits',simard1[matches],clobber=True)
        fits.writeto(prefix+'_Simard2_uat.fits',simard2[matches],clobber=True)
        fits.writeto(prefix+'_Simard3_uat.fits',simard3[matches],clobber=True) 
    else:
        fits.writeto(prefix+'_NSA.fits',nsadat[matches],clobber=True)
        fits.writeto(prefix+'_Petro.fits',petro[matches],clobber=True)
        fits.writeto(prefix+'_Model.fits',model[matches],clobber=True)
        fits.writeto(prefix+'_Simard1.fits',simard1[matches],clobber=True)
        fits.writeto(prefix+'_Simard2.fits',simard2[matches],clobber=True)
        fits.writeto(prefix+'_Simard3.fits',simard3[matches],clobber=True)
    print("--- %s seconds ---" % (time.time() - start_time))
    
if __name__ == '__main__':
    findSurroundings2(arg.r,arg.d,arg.rd,arg.pre)