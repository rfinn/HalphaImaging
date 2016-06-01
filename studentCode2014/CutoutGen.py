# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/obsastro5/.spyder2/.temp.py
"""

#! /usr/bin/env python

import numpy as np
import pylab as pl
from galaxyCutout import cutout
from astropy.io import fits
from astropy.wcs import WCS


filtercenter={'MKW11':6731.,'MKW10':6695.,'MKW8':6731}
filterwidth=80.

class cluster:
    def __init__(self,clustername):
        self.prefix=clustername
        infile='/home/share/research/LocalClusters/NSAmastertables/' + self.prefix + '_NSA.fits'
        hdulist=fits.open(infile)
        self.nsa=hdulist[1].data
        self.filtercenter=filtercenter[self.prefix]
        self.Zmax=(((self.filtercenter + (filterwidth/2))/6563) -1)
        self.Zmin=(((self.filtercenter - (filterwidth/2))/6563) -1)
        
    def plotradec(self):
        pl.figure()
        pl.plot(self.nsa.RA,self.nsa.DEC,'k.')

    def makecuts(self,clustername):
        hdulist= fits.open('/home/share/research/LocalClusters/NSAmastertables/' + self.prefix + '_NSA.fits')
        fdata=hdulist[1].data
        
        # select galaxies that lie in Halpha filter
        self.zFlag = (fdata.Z > self.Zmin)&(fdata.Z < self.Zmax)
        print self.zFlag
        print fdata.Z
        
        #loop through flagged galaxies and get cutouts
        fdulist = fits.open('Hac_final.fits')
        prihdr = fdulist[0].header
        n1 = prihdr['NAXIS1']
        n2 = prihdr['NAXIS2']
        print n1
        print n2


        w = WCS('Hac_final.fits')
        px,py = w.wcs_world2pix(fdata.RA,fdata.DEC,1)
        print px[0:10]
        print fdata.RA[0:10]
        onimageflag=(px < n1) & (px >0) & (py < n2) & (py > 0)
               
        for i in range(len(self.zFlag)):
            if (self.zFlag[i] & onimageflag[i]):
                cutout('test.fits', fdata.RA[i], fdata.DEC[i],xw=100,yw=100, outfile='testout' +str(i) + '.fits')
                
                
        
        #pl.figure()
        #pl.hist(fdata.Z)
        #pl.hist(fdata.Z[zFlag])

mkw11 = cluster('MKW11')
mkw10 = cluster('MKW10')
mkw8 = cluster('MKW8')
