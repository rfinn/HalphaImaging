
# coding: utf-8

# Goal:  
# Make code to take in RA, DEC, and Radius and return a matched catalog with all information on galaxies in that area
# 

# In[2]:

import numpy as np
from astropy.io import fits


# In[67]:

def findSurroundings(ra, dec, radius, prefix):
    petro = fits.getdata('YangDR7PetroToNSA.fits')
    model = fits.getdata('YangDR7ModelToNSA.fits')
    nsadat = fits.getdata('nsa_v0_1_2.fits')
    simard = fits.getdata('Simard1ToNSA.fits')
    d=np.sqrt((ra-nsadat.RA)**2 + (dec-nsadat.DEC)**2)
    index=np.arange(len(d))
    matches=(d<=radius) # list of indices that are within the radius
    print '%i entries found' % len(matches)
    fits.writeto('yang_subcat.fits',petro[matches])
    fits.writeto('simard_subcat.fits',simard[matches])
    nsamatch = []
    petromatch = []
    modelmatch = []
    simardmatch = []
    numcolumns = [len(nsadat[0]),len(petro[0]),len(model[0]),len(simard[0])]
    for i in range(len(matches)):
        for j in range(max(numcolumns)):
            if len(nsamatch)<len(nsadat[0]):
                nsamatch.append([])
            if len(petromatch)<len(petro[0]):
                petromatch.append([])
            if len(modelmatch)<len(model[0]):
                modelmatch.append([])
            if len(simardmatch)<len(simard[0]):
                simardmatch.append([]) 
            try:
                nsamatch[j].append(nsadat[matches[i]][j])
                simardmatch[j].append(simard[matches[i]][j])
                petromatch[j].append(petro[matches[i]][j])
                modelmatch[j].append(model[matches[i]][j])
            except IndexError:
                pass
            
        
    nsacol = []
    petrocol = []
    modelcol = []
    simardcol = []
    for i in range(max(numcolumns)):
        try:
            print nsadat.names[i],nsadat.formats[i],nsadat.dtype[i]
            ncol = fits.Column(name=nsadat.names[i],format=nsadat.formats[i],array=nsamatch[i])
            scol = fits.Column(name=simard.names[i],format=simard.formats[i],array=simardmatch[i])
            pcol = fits.Column(name=petro.names[i],format=petro.formats[i],array=petromatch[i])
            mcol = fits.Column(name=model.names[i],format=model.formats[i],array=modelmatch[i])
            nsacol.append(ncol)
            petrocol.append(pcol)
            modelcol.append(mcol)
            simardcol.append(scol)
        except IndexError:
            pass
    nsadefs = fits.ColDefs(nsacol)
    petrodefs = fits.ColDefs(petrocol)
    modeldefs = fits.ColDefs(modelcol)
    simarddefs = fits.ColDefs(simardcol)
    nsatable = fits.BinTableHDU.from_columns(nsadefs)
    petrotable = fits.BinTableHDU.from_columns(petrodefs)
    modeltable = fits.BinTableHDU.from_columns(modeldefs)
    simardtable = fits.BinTableHDU.from_columns(simarddefs)
    nsatable.writeto(prefix+'_NSA.fits',clobber=True)
    petrotable.writeto(prefix+'_Petro.fits',clobber=True)
    modeltable.writeto(prefix+'_Model.fits',clobber=True)
    simardtable.writeto(prefix+'_Simard.fits',clobber=True)
    


# In[80]:



