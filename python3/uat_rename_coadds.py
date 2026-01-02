#!/usr/bin/env python

'''

User provides:
- filestring (look for *filestring*.fits)
- input directory (defaults to currect directory)
- output directory (defaults to current directory)

output:
- renamed files

USAGE:
in directory: /data-pool/HalphaGroups/rename_test/

python ~/github/HalphaImaging/python3/uat_rename_coadds.py

renaming  ./NRGb177-h02_R.coadd.fits -> ./VF-181.362+20.406-HDI-2015-04-17-NRGb177-h02-R.fits
renaming  ./NRGb177-h02_R.coadd.weight.fits -> ./VF-181.362+20.406-HDI-2015-04-17-NRGb177-h02-R.weight.fits
renaming  ./NRGb177-h02_ha12.coadd.fits -> ./VF-181.362+20.406-HDI-2015-04-17-NRGb177-h02-ha12.fits
renaming  ./NRGb177-h02_ha12.coadd.weight.fits -> ./VF-181.362+20.406-HDI-2015-04-17-NRGb177-h02-ha12.weight.fits



python ~/github/HalphaImaging/python3/uat_rename_coadds.py --outdir /data-pool/HalphaGroups/allcoadsv0


UPDATES:
* 2026 Jan 02:
 mosaic images have underscore in object names, and the renaming script expects a dash.  adding mosaic flag


'''

import os
import numpy as np
import glob
from astropy.io import fits
import argparse

def get_updated_uat_coadd_name(imname,mosflag=False):
    h = fits.getheader(imname)
    ra = float(h['CRVAL1'])
    dec = float(h['CRVAL2'])

    try:
        dateobs_time = h['OBSDATE']
    except KeyError:
        dateobs_time = h['DATE-OBS']        
    dateobs = dateobs_time.split('T')[0]
    dateobs = dateobs.replace('-','')

    temptel = h['INSTRUME']
    if (temptel.find('hdi') > -1) | (temptel.find('HDI') > -1):
        telescope = "HDI"
    elif (temptel.find('mos') > -1) | (temptel.find('MOS') > -1) |  (temptel.find('Mos') > -1) :
        telescope = "MOS"
    else:
        print("WARNING: could not parse header field INSTRUME")
        print("\t setting instrument to HDI")
        telescope = "HDI"

    # let's not switch everything to underscores b/c it's wreaking havoc with halphagui
    try:
        pointing,filterwithsuffix = imname.split('_')
    except ValueError:
        cpointing, hnumber,filterwithsuffix = imname.split('_')
        pointing = f"{cpointing}-{hnumber}"
        
    #pointing = pointing.replace('-','_')

    if mosflag & pointing.startswith('f'):
        # remove preceding 'f' from pointing
        pointing = pointing[1:]
        
    # remove coadd from name
    filterwithsuffix = filterwithsuffix.replace('.coadd','')
    # create string for output name

    if float(dec) < 0:
        #outfile = f'UAT_{ra:07.3f}-{dec:06.3f}_{telescope}_{dateobs}_{pointing}_{filterwithsuffix}'
        outfile = f'UAT-{ra:07.3f}-{np.abs(dec):06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'        
    else:
        #outfile = f'UAT_{ra:07.3f}+{dec:06.3f}_{telescope}_{dateobs}_{pointing}_{filterwithsuffix}'
        outfile = f'UAT-{ra:07.3f}+{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'        
    return outfile


if __name__ == '__main__':
    telescope = 'BOK'
    
    parser = argparse.ArgumentParser(description ="rename UAT coadded images to f'UAT-{ra:07.3f}+{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix} ")
    parser.add_argument('--filestring', dest = 'filestring', default = 'coadd', help = 'filestring to match. default is coadd, which will grab all *coadd*.fits.  Note: weight images will be renamed with their corresponding science coadd.')
        
    parser.add_argument('--indir', dest = 'indir', default = '.', help = 'directory of input images.  Default is current directory')
    parser.add_argument('--outdir', dest = 'outdir', default = '.', help = 'directory to write output images to.  Default is current directory')
    parser.add_argument('--mosaic', dest = 'mosaic', default = False, action='store_true', help = 'set this is the object has an underscore instead of a dash.  e.g. NGC5846_h01 instead of NGC5846-h01')    
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'new filenames are printed when this is set, but the files are not renamed.  This is useful when testing changes to the code.')    

    
    args = parser.parse_args()


    # grab all the coadds in a particular directory
    if args.mosaic: # look for coadds that have been flattened by getzp and have therefore start with an f
        filelist = glob.glob(os.path.join(args.indir,"f*"+args.filestring+"*.fits"))
    else:
        filelist = glob.glob(os.path.join(args.indir,"*"+args.filestring+"*.fits"))
    filelist.sort()


    for fname in filelist:
        if 'weight' in fname:
            continue
        if 'CS-ZP' in fname:
            continue
        new_name = get_updated_uat_coadd_name(os.path.basename(fname), mosflag=args.mosaic)
        # rename as UAT-{RA}-{DEC}-{TEL}-{OBSDATE}-{POINTING}-{FILTER}
        new_output_image = os.path.join(args.outdir, new_name)
        print('renaming ',fname,'->',new_output_image)
        if not args.testing:
            #os.rename(fname,new_output_image)
            os.system(f"cp {fname} {new_output_image}")

        # rename the weight file
        weightfile = fname.replace('.fits','.weight.fits')
        # remove preceding 'f' if in mosaic mode
        if args.mosaic:
            weightfile = weightfile[1:]
        if os.path.exists(weightfile):
            print('renaming ',weightfile,'->',new_output_image.replace('.fits','.weight.fits'))
            if not args.testing:
                #os.rename(weightfile,new_output_image.replace('.fits','.weight.fits'))
                os.system(f"cp {weightfile} {new_output_image.replace('.fits','.weight.fits')}")
        else:
            print("WARNING: weight image does not exist for ",fname, weightfile)



