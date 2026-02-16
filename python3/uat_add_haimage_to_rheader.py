#!/usr/bin/env python

"""
Run in the coadd directory e.g. /data-pool/HalphaGroups/allcoaddsv2

get list of *R.fits images

look for corresponding halpha image [ha4, ha8, ha12, ha16]

"""

import glob
import os
import argparse
from astropy.io import fits

def get_params_from_name_uat(image_name):
    '''
    coadd names should be as follows.
    if float(dec) < 0:
        outfile = f'UAT-{ra:07.3f}-{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'        
    else:
        outfile = f'UAT-{ra:07.3f}+{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'    

    each pointing will have an internal '-', like ABELL1367-h01
    '''
    t = os.path.basename(image_name).split('-')
    #print(t)
    if len(t) == 7:
        ra, dec = t[1].split('+')
        telescope = t[2]
        dateobs = t[3]
        pointing = t[4]+'-'+t[5]
    elif len(t) == 8: # meant to catch negative declinations
        ra = t[1]
        dec = -1*t[2]
        telescope = t[3]
        dateobs = t[4]
        pointing = t[5]+'-'+t[6]
        
    else:
        print("ruh roh - trouble getting info from ",image_name, len(t))
        print(image_name)
        print(t)
        return
    return telescope,dateobs,pointing,ra,dec

def get_params_from_name_vfs(image_name):
    '''
    coadd names should be as follows.
    if float(dec) < 0:
        outfile = f'UAT-{ra:07.3f}-{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'        
    else:
        outfile = f'UAT-{ra:07.3f}+{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'    

    each pointing will have an internal '-', like ABELL1367-h01
    '''
    t = os.path.basename(image_name).split('-')
    #print(t)
    if len(t) == 6:
        ra, dec = t[1].split('+')
        telescope = t[2]
        dateobs = t[3]
        pointing = t[4]
    elif len(t) == 7: # meant to catch negative declinations
        ra = t[1]
        dec = -1*t[2]
        telescope = t[3]
        dateobs = t[4]
        pointing = t[5]
        
    else:
        print("ruh roh - trouble getting info from ",image_name, len(t))
        #print(image_name)
        print("t = ",t, len(t))
        return
    return telescope,dateobs,pointing, ra, dec

def add_field_2_header(ffile, field, value):
    """
    PARAMS:
    * ffile - file to update header with
    * field - new field for header
    * value - value to put in the header card 
    """
    
    hdu = fits.open(ffile)
    hdu[0].header.set(field,value)
    hdu.writeto(ffile,overwrite=True)
    hdu.close()

if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description ="This program adds the Halpha imaging name to the r-band images, searching for all UAT*-R.fits for the r-band images.")
    parser.add_argument('--filestring', dest = 'filestring', default = 'UAT', help = 'filestring to match. default is UAT, which will grab all UAT*R.fits files. Note: weight images will be renamed with their corresponding science coadd.')
    #parser.add_argument('--se', dest = 'se', default = False, action='store_true', help = 'Run source extractor')    
    #parser.add_argument('--scamp', dest = 'scamp', default = False, action='store_true', help = 'Run scamp')
    #parser.add_argument('--submedian', dest = 'submedian', default = False, action='store_true', help = 'Subtract a median sky value from each image.')
    parser.add_argument('--verbose', dest = 'verbose', default = False, action='store_true', help = 'Print more information/diagnostic statements.')
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'Will print matches but will not edit headers.')
    parser.add_argument('--vfs', dest = 'vfs', default = False, action='store_true', help = 'Set this if running of VFS coadds (different naming convention).')    
    
    args = parser.parse_args()

    halpha_options = [f"ha{i}" for i in [4,8,12,16]]
    halpha_options = halpha_options + ['Ha4','Halpha','Ha6657'] # adding more options so we can use with Virgo VFS images
    print("halpha filters: ",halpha_options)

    rfilelist1 = glob.glob(f"{args.filestring}*R.fits")
    rfilelist2 = glob.glob(f"{args.filestring}*r.fits")
    rfilelist = rfilelist1 + rfilelist2

    rfilelist.sort()

    # loop through r-band coadds

    for rfile in rfilelist:
        print()

        # get the target name
        if args.vfs:
            telescope, dateobs, pointing, ra, dec = get_params_from_name_vfs(rfile)
        else:
            telescope, dateobs, pointing, ra, dec = get_params_from_name_uat(rfile)
            print(f"I think: telescope={telescope}, dateobs={dateobs}, pointing={pointing}")
        print(f"working on file {rfile}, target={pointing}")
        # look for target name with halpha image
        for h in halpha_options:

            if args.vfs:
                # don't include dateobs because it could be different for r vs halpha
                searchstring = f"{args.filestring}*{ra}*{dec}*{telescope}*{pointing}*{h}.fits"
                if args.verbose:
                    print(searchstring)
            else:
                searchstring = f'UAT*{ra}*{dec}*{pointing}*-{h}.fits'
                if args.testing:
                    print(f"looking for halpha image in filter {h} using search string = {searchstring}")
            hfiles = glob.glob(searchstring)

            if len(hfiles) == 1:
                # add hfile to rband image header
                print(f"{rfile}: adding HAIMAGE={hfiles[0]} to header")
                if not args.testing:
                    add_field_2_header(rfile,'HAIMAGE',hfiles[0])

                # add the rband image to the halpha image
                print(f"{hfiles[0]}: adding RIMAGE={rfile} to header")
                if not args.testing:
                    add_field_2_header(hfiles[0],'RIMAGE',rfile)            
