#!/usr/bin/env python

'''

User provides:
- input directory (defaults to currect directory)
- output directory (defaults to current directory)
- instrument

'''

import os
import glob
from astropy.io import fits
import argparse

def get_updated_uat_coadd_name(imname):
    h = fits.getheader(imname)
    ra = float(h['CRVAL1'])
    dec = float(h['CRVAL2'])
    
    dateobs_time = h['OBSDATE']
    dateobs = dateobs_time.split('T')[0]

    temptel = h['INSTRUME']
    if (temptel.find('hdi') > -1) | (temptel.find('HDI') > -1):
        telescope = "HDI"
    elif (temptel.find('mos') > -1) | (temptel.find('MOS') > -1):
        telescope = "MOS"
    else:
        print("WARNING: could not parse header field INSTRUME")
        print("\t setting instrument to HDI")
        
    pointing,filterwithsuffix = imname.split('_')

    # remove coadd from name
    filterwithsuffix = filterwithsuffix.replace('.coadd','')
    # create string for output name

    if float(dec) < 0:
        outfile = f'VF-{ra:07.3f}-{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'
    else:
        outfile = f'VF-{ra:07.3f}+{dec:06.3f}-{telescope}-{dateobs}-{pointing}-{filterwithsuffix}'
    return outfile


if __name__ == '__main__':
    telescope = 'BOK'
    
    parser = argparse.ArgumentParser(description ='stack the 90prime images after noao pipeline')
    parser.add_argument('--filestring', dest = 'filestring', default = 'coadd', help = 'filestring to match. default is coadd, which will grab all *coadd*.fits.  Note: weight images will be renamed with their corresponding science coadd.')
        
    parser.add_argument('--indir', dest = 'indir', default = '.', help = 'directory of input images.  Default is current directory')
    parser.add_argument('--outdir', dest = 'outdir', default = '.', help = 'directory to write output images to.  Default is current directory')
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'new filenames are printed when this is set, but the files are not renamed.  This is useful when testing changes to the code.')    

    
    args = parser.parse_args()


    # grab all the coadds in a particular directory
    filelist = glob.glob(os.path.join(args.indir,"*"+args.filestring+"*.fits"))
    filelist.sort()


    for fname in filelist:
        if 'weight' in fname:
            continue
        if 'CS-ZP' in fname:
            continue
        new_name = get_updated_uat_coadd_name(os.path.basename(fname))
        # rename as UAT-{RA}-{DEC}-{TEL}-{OBSDATE}-{POINTING}-{FILTER}
        new_output_image = os.path.join(args.outdir, new_name)
        print('renaming ',fname,'->',new_output_image)
        if not args.testing:
            os.rename(fname,new_output_image)

        # rename the weight file
        weightfile = fname.replace('.fits','.weight.fits')
        if os.path.exists(weightfile):
            print('renaming ',weightfile,'->',new_output_image.replace('.fits','.weight.fits'))
            if not args.testing:
                os.rename(weightfile,new_output_image.replace('.fits','.weight.fits'))
        else:
            print("WARNING: weight image does not exist for ",fname)



