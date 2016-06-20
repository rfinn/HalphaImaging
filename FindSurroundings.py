#!/usr/bin/python
# coding: utf-8

# Goal:  
# Make code to take in RA, DEC, and Radius and return a matched catalog with all information on galaxies in that area
# 


import numpy as np
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(description ='Creates line matched tables for area surrounding input RA and DEC')
parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Use shortened version of NSA catalog (nsa_uat.fits)')
parser.add_argument('--RA', dest = 'r', default = 0, help = 'Target RA')
parser.add_argument('--DEC', dest = 'd', default = 0, help = 'Target DEC')
parser.add_argument('--radius', dest = 'rd', default = 1, help = 'Radius around target to include in tables')
parser.add_argument('--prefix', dest = 'pre', default = 'test', help = 'Prefix of new file names')
args = parser.parse_args()


def findSurroundings(ra, dec, radius, prefix):
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
    matches=(d<=radius) # list of indices that are within the radius
    print '%i entries found' % len(matches)
    nsamatch = []
    petromatch = []
    modelmatch = []
    simard1match = []
    numcolumns = [len(nsadat[0]),len(petro[0]),len(model[0]),len(simard1[0]),len(simard2[0]),len(simard3[0])]
    for i in range(len(matches)):
        for j in range(max(numcolumns)):
            if len(nsamatch)<len(nsadat[0]):
                nsamatch.append([])
            if len(simard1match)<len(simard1[0]):
                simard1match.append([]) 
            if len(simard2match)<len(simard2[0]):
                simard2match.append([])
            if len(simard3match)<len(simard3[0]):
                simard3match.append([])
            if len(petromatch)<len(petro[0]):
                petromatch.append([])
            if len(modelmatch)<len(model[0]):
                modelmatch.append([])
            
            try:
                nsamatch[j].append(nsadat[matches[i]][j])
                simard1match[j].append(simard1[matches[i]][j])
                simard2match[j].append(simard2[matches[i]][j])
                simard3match[j].append(simard3[matches[i]][j])
                petromatch[j].append(petro[matches[i]][j])
                modelmatch[j].append(model[matches[i]][j])
            except IndexError:
                pass
    nsacol = []
    petrocol = []
    modelcol = []
    simard1col = []
    simard3col = []
    simard3col = []
    for i in range(max(numcolumns)):
        try:
            print nsadat.names[i],nsadat.formats[i],nsadat.dtype[i]
            ncol = fits.Column(name=nsadat.names[i],format=nsadat.formats[i],array=nsamatch[i])
            nsacol.append(ncol)
            s1col = fits.Column(name=simard1.names[i],format=simard1.formats[i],array=simard1match[i])
            simard1col.append(s1col)
            s2col = fits.Column(name=simard2.names[i],format=simard2.formats[i],array=simard2match[i])
            simard2col.append(s2col)
            s3col = fits.Column(name=simard3.names[i],format=simard3.formats[i],array=simard3match[i])
            simard3col.append(s3col)
            pcol = fits.Column(name=petro.names[i],format=petro.formats[i],array=petromatch[i])
            petrocol.append(pcol)
            mcol = fits.Column(name=model.names[i],format=model.formats[i],array=modelmatch[i])
            modelcol.append(mcol)   
        except IndexError:
            pass
    nsadefs = fits.ColDefs(nsacol)
    petrodefs = fits.ColDefs(petrocol)
    modeldefs = fits.ColDefs(modelcol)
    simard1defs = fits.ColDefs(simard1col)
    simard2defs = fits.ColDefs(simard2col)
    simard3defs = fits.ColDefs(simard3col)
    nsatable = fits.BinTableHDU.from_columns(nsadefs)
    petrotable = fits.BinTableHDU.from_columns(petrodefs)
    modeltable = fits.BinTableHDU.from_columns(modeldefs)
    simard1table = fits.BinTableHDU.from_columns(simard1defs)
    simard2table = fits.BinTableHDU.from_columns(simard2defs)
    simard3table = fits.BinTableHDU.from_columns(simard3defs)
    if arg.s:
        nsatable.writeto(prefix+'_NSA_uat.fits',clobber=True)
        petrotable.writeto(prefix+'_Petro_uat.fits',clobber=True)
        modeltable.writeto(prefix+'_Model_uat.fits',clobber=True)
        simard1table.writeto(prefix+'_Simard1_uat.fits',clobber=True)
        simard2table.writeto(prefix+'_Simard2_uat.fits',clobber=True)
        simard3table.writeto(prefix+'_Simard3_uat.fits',clobber=True)
    else:
        nsatable.writeto(prefix+'_NSA.fits',clobber=True)
        petrotable.writeto(prefix+'_Petro.fits',clobber=True)
        modeltable.writeto(prefix+'_Model.fits',clobber=True)
        simard1table.writeto(prefix+'_Simard1.fits',clobber=True)
        simard2table.writeto(prefix+'_Simard2.fits',clobber=True)
        simard3table.writeto(prefix+'_Simard3.fits',clobber=True)
        
if __name__ == '__main__':
    findSurroundings(arg.r,arg.d,arg.rd,arg.pre)

