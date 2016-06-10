
# coding: utf-8

# In[2]:

import numpy as np
from astropy.io import fits
import time


# In[3]:

def findSurroundings2(ra, dec, radius, prefix):
    start_time = time.time()
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
    fits.writeto(prefix+'_NSA.fits',nsadat[matches],clobber=True)
    fits.writeto(prefix+'_Petro.fits',petro[matches],clobber=True)
    fits.writeto(prefix+'_Model.fits',model[matches],clobber=True)
    fits.writeto(prefix+'_Simard1.fits',simard1[matches],clobber=True)
    fits.writeto(prefix+'_Simard2.fits',simard2[matches],clobber=True)
    fits.writeto(prefix+'_Simard3.fits',simard3[matches],clobber=True)
    print("--- %s seconds ---" % (time.time() - start_time))


# In[ ]:



