
# coding: utf-8

# GOAL:  
#     To complete the Yang tables so they include information on each galaxy and the group it is in

# In[3]:

import csv
import numpy as np
from astropy.io import fits
import glob


# In[4]:

def combineTable(t1,t2):
    t1 = np.array(t1) # SDSS7.csv
    t2 = np.array(t2) # combinePetro.csv / combineModel.csv
    combinedArray = np.empty((639359,28),'f') 
    combinedArray[:,0:13] = t1[0:13,:].T
    combinedArray[:,13:28] = t2[2:17,:].T
    return combinedArray

def makeFits(cT,outname):
    p0 = fits.Column(name='galaxy_ID',format='J',array=cT[:,0])
    p1 = fits.Column(name='galaxy_ID_NYU',format='J',array=cT[:,1])
    p2 = fits.Column(name='RA_gal',format='E',array=cT[:,2])
    p3 = fits.Column(name='DEC_gal',format='E',array=cT[:,3])
    p4 = fits.Column(name='Z_gal',format='E',array=cT[:,4])
    p5 = fits.Column(name='app_mag',format='E',array=cT[:,5])
    p6 = fits.Column(name='mag_lim',format='E',array=cT[:,6])
    p7 = fits.Column(name='completeness',format='E',array=cT[:,7])
    p8 = fits.Column(name='log_h_petro',format='E',array=cT[:,8])
    p9 = fits.Column(name='color_petro',format='E',array=cT[:,9])
    p10 = fits.Column(name='log_h_model',format='E',array=cT[:,10])
    p11 = fits.Column(name='color_model',format='E',array=cT[:,11])
    p12 = fits.Column(name='Z_source',format='J',array=cT[:,12])
    p13 = fits.Column(name='group_ID',format='J',array=cT[:,13])
    p14 = fits.Column(name='brightest_gal',format='J',array=cT[:,14])
    p15 = fits.Column(name='massive_gal',format='J',array=cT[:,15])
    p16 = fits.Column(name='RA_group',format='E',array=cT[:,16])
    p17 = fits.Column(name='DEC_group',format='E',array=cT[:,17])
    p18 = fits.Column(name='Z_group',format='E',array=cT[:,18])
    p19 = fits.Column(name='L_group',format='E',array=cT[:,19])
    p20 = fits.Column(name='mass_group',format='E',array=cT[:,20])
    p21 = fits.Column(name='halo_mass1',format='E',array=cT[:,21])
    p22 = fits.Column(name='halo_mass2',format='E',array=cT[:,22])
    p23 = fits.Column(name='d_L',format='E',array=cT[:,23])
    p24 = fits.Column(name='d_M',format='E',array=cT[:,24])
    p25 = fits.Column(name='f_edge',format='E',array=cT[:,25])
    p26 = fits.Column(name='ID1',format='J',array=cT[:,26])
    p27 = fits.Column(name='ID2',format='J',array=cT[:,27])
    
    
    coldefs = fits.ColDefs([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27])
    hdu = fits.BinTableHDU.from_columns(coldefs)
    outfile = outname+'.fits'
    hdu.writeto(outfile,clobber=True)


# In[5]:

doc='SDSS7' 
with open(doc) as fin, open(doc+'.csv','w') as fout:
    o=csv.writer(fout)
    for line in fin:
        o.writerow(line.split())


# In[6]:

petro = np.loadtxt('combinePetro.csv',delimiter=',',unpack=True,dtype=float)
model = np.loadtxt('combineModel.csv',delimiter=',',unpack=True,dtype=float)
SDSS7 = np.loadtxt('SDSS7.csv',delimiter=',',unpack=True,dtype=float)


# In[ ]:

completePetro = combineTable(SDSS7,petro)
completeModel = combineTable(SDSS7,model)
makeFits(completePetro,'completePetro')
makeFits(completeModel,'completeModel')
np.savetxt("completePetro.csv", completePetro, delimiter=",")
np.savetxt("completeModel.csv", completeModel, delimiter=",")


# In[ ]:



