
# coding: utf-8

# Goal:  
# To match the NSA catalog with both the Yang Petro catalog and the Yang Model catalog.  
# 
# Required files:  
# * Must have the output files from both CombineGroupYang.ipynb and CompleteYang.ipynb in the current directory
# * Specifically completePetro.fits and completeModel.fits
# * The NASA Sloan Atlas catalog (nsa_v0_1_2.fits) found here http://www.nsatlas.org/data
# 
# Outputs:
# * YangDR7PetroToNSA.fits 
# * YangDR7ModelToNSA.fits

# In[1]:

import csv
import numpy as np
from astropy.io import fits
import glob


# In[2]:

def findnearest(x1,y1,x2,y2,delta):
    matchflag=1
    nmatch=0
    d=np.sqrt((x1-x2)**2 + (y1-y2)**2)
    index=np.arange(len(d))
    t=index[d<delta]
    matches=t
    if len(matches) > 0:
        nmatch=len(matches)
        if nmatch > 1:
            imatch=index[(d == min(d[t]))]
        else:
            imatch=matches[0]
    else:
        imatch = 0
        matchflag = 0

    return imatch, matchflag,nmatch


# In[4]:

petro = fits.getdata('completePetro.fits')
model = fits.getdata('completeModel.fits')
nsadat =fits.getdata('nsa_v0_1_2.fits')


# In[7]:

matchRadius=0.1/3600

imatch=np.zeros(len(nsadat.RA),'i')
matchflag=np.zeros(len(nsadat.RA),'bool')
nmatch=np.zeros(len(nsadat.RA),'i')
for i in range(len(nsadat.RA)):
    imatch[i],matchflag[i],nmatch[i]  = findnearest(nsadat.RA[i],nsadat.DEC[i],petro.field(2),petro.field(3),matchRadius)      
outfile='YangDR7PetroToNSA.fits'
matchedarray=np.zeros(len(nsadat),dtype=petro.dtype)
matchedarray[matchflag] = petro[imatch[matchflag]]
fits.writeto(outfile,matchedarray,clobber=True)


# In[8]:

imatch=np.zeros(len(nsadat.RA),'i')
matchflag=np.zeros(len(nsadat.RA),'bool')
nmatch=np.zeros(len(nsadat.RA),'i')
for i in range(len(nsadat.RA)):
    imatch[i],matchflag[i],nmatch[i]  = findnearest(nsadat.RA[i],nsadat.DEC[i],model.field(2),model.field(3),matchRadius)      
outfile='YangDR7ModelToNSA.fits'
matchedarray=np.zeros(len(nsadat),dtype=model.dtype)
matchedarray[matchflag] = model[imatch[matchflag]]
fits.writeto(outfile,matchedarray,clobber=True)

