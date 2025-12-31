#!/usr/bin/env python

'''
adding INSTRMNT keyword to header to all WFC*PA.fits files for INT data so scamp can fit the geometric distortion for each chip separately.

USAGE:

python ~/github/Virgo/programs/INT_add_inst_header.py 

'''


from astropy.io import fits
from astropy import wcs
import glob
import os

matchstrings = ['WFC*1PA.fits','WFC*2PA.fits','WFC*3PA.fits','WFC*4PA.fits']
instruments = ['INTWFC1','INTWFC2','INTWFC3','INTWFC4']
print(f"\nWorking on directory {os.getcwd()}")
for i in range(len(matchstrings)):
    files = glob.glob(matchstrings[i])
    print('chip ',i+1,' updating ',len(files),' files')
    j=0
    for f in files:
        #print(j)
        # read in image and header for the image that needs to be updated
        
        hdu = fits.open(f)
        #hdu.verify()

        # trying again after implementing two commands that might fix the issue with CD1_1
        hdu[0].header.set('INSTRMNT',instruments[i])
        try:
            hdu.writeto(f,overwrite=True,output_verify='ignore')
            
        except TypeError: # getting buffer too small - maybe files are corrupted?
            print(f"OH NO! astropy won't write {f}, file may be corrupted.")
            if '2019feb' in os.getcwd():
                # look for earlier version
                print(f"\tlooking for an earlier version of file {f}")
                try:
                    if 'WFC.r.' in f:
                        altfile = f"/data-pool/Halpha/processed/2019febINT-all/allimages/rband/{f}"
                    else:
                        altfile = f"/data-pool/Halpha/processed/2019febINT-all/allimages/Halpha/{f}"
                    if os.path.exists(altfile):
                        hdu.close() # close first attempt
                        print(f"found {altfile}!")
                        os.system(f"cp {altfile} .")
                        
                        hdu = fits.open(f)
                        hdu[0].header.set('INSTRMNT',instruments[i])
                        hdu.writeto(f,overwrite=True,output_verify='ignore')
                        print('\talt version seemed to work!')
                    else:
                        print("\tcould not find alt version of file ",altfile)
                except TypeError:
                    print("\tso sorry, but was not able to get alt version of file to work")
                    # now try 
        hdu.close()
        j+= 1
