#!/usr/bin/env python

'''
BASIC INFORMATION ABOUT THIS CODE:
  -This is the FIFTH program you need to run in order to perform the reduction o  f HDI images.
  -This code will fix all the headers to contain basic WCS information so that i  n the future we can feed the fits images into stacking programs like SCamp and  SWarp. 
  -Updates CMMTOBS --> FILTER
        RASTRNG --> CRVAL1
        DECSTRNG --> CRVAL2
  -Adds CRPIX1, CRPIX2, CD1_1, CD2_2, CTYPE1, CTYPE2

BEFORE RUNNING THIS CODE:
  -Be sure to type ur_setup in the terminal everytime you open a new terminal     window so that ureka is activated.
  -Ensure that pyraf is still activated by retyping the commands listed in the c  omments of the FIRST program titled "uat_HDIgroupflatfiles.py".

PROCEDURE:
  -This code uses the python task .rename and .append to the list of header information so that the header contains basic WCS information that can be easily read by astronomy programs. 
  -In addition, we used different arguments using args.parse to call on adding the different header information.
     
GOAL:
    The goal of this code is to successfully add header information that is not included in the HDI raw data frames.

EXAMPLE:
   In the directory containing all flattened objects with incorrect headers type in the command line:
      /home/share/research/pythonCode/uat_HDIfixheader.py

      (or whatever the path is to where this program is stored)
   
INPUT/OUPUT:
    Input --> all cftr*.fits in directory
    Output --> hcftr*.fits

REQUIRED MODULES:
    argparse
    glob
    astropy
    
EXTRA NOTES:
can be used on dome flattened images. To do so type '--filestring "dtr*.fits"' after the command

WRITTEN BY:
Rose Finn

EDITED BY:
Research Team Summer 2015 --> Grant Boughton, Natasha Collova, Tiffany Flood, Kaitlyn Hoag, Kelly Whalen

'''

import argparse
import glob
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

parser = argparse.ArgumentParser(description ='Edit image headers to include basic WCS information to the HDI image headers')
parser.add_argument('--filestring', dest='filestring', default='d', help='match string for input files (default =  d, which operates on all d*.fits files)')
parser.add_argument('--pixscalex', dest='pixelscalex', default='0.00011805', help='pixel scale in x (default = 0.00011808)')
parser.add_argument('--pixscaley', dest='pixelscaley', default='0.00011805', help='pixel scale in y (default = 0.00011808)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)
i=1
for f in files:
    print 'FIXING HEADER FOR FILE %i OF %i'%(i,nfiles)
    data, header = fits.getdata(f,header=True)
    header.rename_keyword('FILTER1','FWHEEL1')
    header.rename_keyword('FILTER2','FWHEEL2')

    FILTER = header['CMMTOBS']

    try:    
        ccdsec = header['CCDSEC']
    except:
        print 'WARNING'
        print '\t no CCDSEC in header'
        print '\t assuming HDI values'
        ccdsec = '[1:4095, 1:4109]'


    RASTRNG = header['RASTRNG']
    DECSTRNG = header['DECSTRNG']
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    RA = coord.Angle(RASTRNG,unit=u.hour)
    DEC = coord.Angle(DECSTRNG,unit=u.degree)
    EQUINOX = header['EQUINOX']
    INSTRUMENT = header['INSTRUME']
    OBJECT = header['OBJECT']
    EXPTIME = header['EXPTIME']
    # process coordinates to J2000 epoch
    c = SkyCoord(ra=RA.deg*u.degree,dec=DEC.degree*u.degree,frame='fk5',equinox='J'+str(EQUINOX))
    #print 'original coords = ',c
    c2000 = c.transform_to(coord.FK5(equinox='J2000.0'))
    #print 'J2000 coords = ',c2000
    # trying to clear header and add the min info required
    header.clear()
    header.append(card=('FILTER',FILTER,'FILTER'))
    header.append(card=('DATASEC',ccdsec,'DATA SECTION'))
    header.append(card=('INSTRUME',INSTRUMENT,'INSTRUMENT'))
    if 'CRVAL1' in header:
        header['CRVAL1']=(c.ra.value,'RA of reference point')
    else:
        header.append(card=('CRVAL1',c.ra.value,'RA of reference point'))

    if 'CRVAL2' in header:
        header['CRVAL2']=(c.dec.value,'DEC of reference point')
    else:
        header.append(card=('CRVAL2',c.dec.value,'DEC of reference point'))

    if 'CRPIX1' in header:
        header['CRPIX1']=(1990.,'X reference pixel')
        header['CRPIX2']=(2048.,'Y reference pixel')
    else:
        header.append(card=('CRPIX1',1990.,'X reference pixel'))
        header.append(card=('CRPIX2',2048.,'Y reference pixel'))

    if 'CD1_1' in header:
        header['CD1_1']=(float(args.pixelscalex),'Pixel scale in X')
        header['CD2_2']=(float(args.pixelscaley),'Pixel scale in Y')
    else:
        header.append(card=('CD1_1',float(args.pixelscalex),'Pixel scale in X'))
        header.append(card=('CD2_2',float(args.pixelscaley),'Pixel scale in Y'))
    if 'CD1_2' in header:
        header['CD1_2']=(0,'')# value from Becky's quick_WCS_instructions
        header['CD2_1']=(0,'')# value from Becky's quick_WCS_instructions
    else:
        header.append(card=('CD1_2',1.226985e-6,''))# value from Becky's quick_WCS_instructions
        header.append(card=('CD2_1',-1.2404496e-6,''))# value from Becky's quick_WCS_instructions
    if 'CTYPE1' in header:
        header['CTYPE1']=('RA---TAN','')# value from Becky's quick_WCS_instructions
        header['CTYPE2']=('DEC--TAN','')# value from Becky's quick_WCS_instructions
    else:
        header.append(card=('CTYPE1','RA---TAN',''))# value from Becky's quick_WCS_instructions
        header.append(card=('CTYPE2','DEC--TAN',''))# value from Becky's quick_WCS_instructions
    #if 'EQUINOX_OBS' in header:
    #    header['EQUINOX_OBS'] = EQUINOX
    #else:
    #    header.append(card=('EQUINOX_OBS',EQUINOX,'Equinox at time of observations'))
    header.append(card=('EQUINOX', EQUINOX,'equinox of RA,Dec'))
    header.append(card=('EPOCH', EQUINOX,'equinox of RA,Dec'))
    header.append(card=('OBJECT', OBJECT,'Object Name'))
    header.append(card=('EXPTIME', EXPTIME,'Exposure Time (sec)'))
    #if 'WCSDIM' in header:
    #    header['WCSDIM']=(2,'')
    #else:
    #    header.append(card=('WCSDIM',2,''))
    #if 'WAT0_001' in header:
    #    header['WAT0_001']=('system=image','')
    #else:
    #    header.append(card=('WAT0_001','system=image',''))
        
    #if 'WAT0_002' in header:
    #    header['WAT0_002']=('system=image','')
    #else:
    #    header.append(card=('WAT0_002','system=image',''))
    
    if 'GAIN' in header:
        header['GAIN']=('1.3','gain (e-/ADU)')
    else:
        header.append(card=('GAIN','1.3','gain (e-/ADU)'))
    
    # add code to transform RA and Dec to epoch 2000
    # c = SkyCoord(ra=123.67*u.degree, dec = 19.5545*u.degree,obstime='J2015.25',frame='fk5',equinox='J2015.25')
    # for mosaic
    #if 'TELRA' in header:
    #    header['TELRA']=(c2000.ra.value,'RA of reference point')
    #    header['TELDEC']=(c2000.dec.value,'DEC of reference point')
    #    header['TELEQUIN']=('J2000.0','Epoch (years)')
    #else:
    #    header.append(card=('TELRA',c2000.ra.value,'RA of reference point'))
    #    header.append(card=('TELDEC',c2000.dec.value,'DEC of reference point'))
    #    header.append(card=('TELEQUIN','J2000.0','Epoch (years)'))    
    print 'WRITING UPDATED FILE'
    fits.writeto('h'+f,data,header,clobber=True)
    i += 1
    print '\n'

    
    


                   
