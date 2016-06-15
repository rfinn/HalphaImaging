
# coding: utf-8

# GOAL:<br/>
#    To combine the Yang group and Yang galaxy catalogs into one large catalog containing all their information
# 
# PROCEDURE:
#     1. Convert space seperated variable tables to csv for easier use
#     2. Read in new csv files
#     3. Combine tables
#     4. Write new csv files containing the combined tables
# 
# REQUIRED MODULES:
#     1. Python's built in csv module
#     2. Numerical Python (numpy)
#     3. Astropy
#     3. Glob
# 
# REQUIRED FILES:  
#     The original SSDS DR7 Yang group and galaxy catalogs found here: http://gax.shao.ac.cn/data/Group.html  
#     Only uses the "C" set of data because that set has the most accurate redshift data  
#     NOTE: currently only works on these files  
#     
# OUTPUTS:  
#     combinePetro.csv (this combines ipetroC_1 and petroC_group)  
#     combinePetro.fits  
#     combineModel.csv (this combines imodelC_1 and modelC_group)  
#     combineModel.fits      
#     (can be changed in the 5th cell in the np.savetxt command)  
#     These files contain the following information by column:
#     1. Galaxy ID
#     2. Galaxy ID in NYU_VAGC
#     3. Group ID (0 if not in any group)
#     4. Brightest Galaxy in group (1 = True, 2 = False)
#     5. Most Massive Galaxy in group (M_stellar) (1 = True, 2 = False)
#     6. RA of group
#     7. DEC of group
#     8. Z of group
#     9. group L_{-19.5} log L_{\odot}/h^2; (characteristic group luminosity)
#     10. group stellar mass of galaxies brighter than M_r<-19.5      
#     11. halo mass1 in log M_halo/ (M_{\odot}/h); (estimated using the ranking of group L_{-19.5})
#     12. halo mass2 in log M_halo/ (M_{\odot}/h); (estimated using the ranking of M_stellar)
#     13. mean separation d=(1/n)^{1/3} Mpc/h of groups brighter than this; (estimated using the ranking of group L_{-19.5})
#     14. mean separation d=(1/n)^{1/3} Mpc/h of groups more massive than this;  (estimated using the ranking of M_stellar)
#     15. f_{edge}; edge effect: See Yang et al. 2007 
#     16. ID1: 1: mass estimated using the ranking of L_{-19.5}, -1: extrapolation
#     17. ID2: 1: mass estimated using the ranking of M_stellar, -1: extrapolation 
#       

# In[97]:

import csv
import numpy as np
from astropy.io import fits
import glob


# In[1]:

# Combines the two inputs into one table, matched by group ID
# Returns a numpy array

def combineTable(t1,t2):
    t1 = np.array(t1)
    t2 = np.array(t2)
    combinedArray = np.empty((639359,17),'f') # creates an empty numpy array that is the same size as the final table
    i = 0
    combinedArray[:,0:5] = t1[:,0:5]
    for row in combinedArray:
        if row[2] !=0: # only matches galaxies that belong in a group
            combinedArray[i,5:18] = t2[1:13,row[2]+2] # uses row[2]+2 because there is an offset of 2 between group ID and index
        i= i+1
    return combinedArray

def makeFits(cT,outname):
    p0 = fits.Column(name='galaxy_ID',format='J',array=cT[:,0])
    p1 = fits.Column(name='galaxy_ID_NYU',format='J',array=cT[:,1])
    p2 = fits.Column(name='group_ID',format='J',array=cT[:,2])
    p3 = fits.Column(name='brightest_gal',format='J',array=cT[:,3])
    p4 = fits.Column(name='massive_gal',format='J',array=cT[:,4])
    p5 = fits.Column(name='RA_group',format='E',array=cT[:,5])
    p6 = fits.Column(name='DEC_group',format='E',array=cT[:,6])
    p7 = fits.Column(name='Z_group',format='E',array=cT[:,7])
    p8 = fits.Column(name='L_group',format='E',array=cT[:,8])
    p9 = fits.Column(name='mass_group',format='E',array=cT[:,9])
    p10 = fits.Column(name='halo_mass1',format='E',array=cT[:,10])
    p11 = fits.Column(name='halo_mass2',format='E',array=cT[:,11])
    p12 = fits.Column(name='d_L',format='E',array=cT[:,12])
    p13 = fits.Column(name='d_M',format='E',array=cT[:,13])
    p14 = fits.Column(name='f_edge',format='E',array=cT[:,14])
    p15 = fits.Column(name='ID1',format='J',array=cT[:,15])
    p16 = fits.Column(name='ID2',format='J',array=cT[:,16])
    coldefs = fits.ColDefs([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16])
    hdu = fits.BinTableHDU.from_columns(coldefs)
    outfile = outname+'.fits'
    hdu.writeto(outfile,clobber=True)


# In[3]:

# Converts space seperated variable tables to .csv files

infiles = glob.glob('*C_*')
infiles = set(infiles)-set(glob.glob('*.csv')) # doesn't convert files that are already .csv
for doc in infiles:
    with open(doc) as fin, open(doc+'.csv','w') as fout:
        o=csv.writer(fout)
        for line in fin:
            o.writerow(line.split())


# In[101]:

# Opens the new .csv files

petro1 = np.loadtxt('ipetroC_1.csv',delimiter=',',unpack=True,dtype=float)
petro3 = np.loadtxt('petroC_group.csv',delimiter=',',unpack=True,dtype=float)
model1 = np.loadtxt('imodelC_1.csv',delimiter=',',unpack=True,dtype=float)
model3 = np.loadtxt('modelC_group.csv',delimiter=',',unpack=True,dtype=float)


# In[142]:

# Writes out the combined tables to a .csv file

combinePetro = combineTable(petro1,petro3)
combineModel = combineTable(model1,model3)
makeFits(combinePetro,'combinePetro')
makeFits(combineModel,'combineModel')
np.savetxt("combinePetro.csv", combinePetro, delimiter=",")
np.savetxt("combineModel.csv", combineModel, delimiter=",")

