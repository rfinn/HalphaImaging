#!/usr/bin/env python

"""
trying to help scamp solve for the solutions by putting images back into MEF format

"""

from astropy.io import fits
import os
import sys
import glob


def build_mef_from_ccds(base_name, ccd_ids=(1, 2, 3, 4), suffix="PA"):
    """
    from chatgpt

    Reconstruct a multi-extension FITS file from individual CCD images.

    Parameters
    ----------
    base_name : str
        Base filename before _XPA.fits
    ccd_ids : iterable
        CCD identifiers (default: 1–4)
    suffix : str
        Suffix string (default: 'PA')

        
    """

    hdulist = fits.HDUList()

    # Create a minimal primary HDU (no data)
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header['EXTEND'] = True
    hdulist.append(primary_hdu)

    for ccd in ccd_ids:
        fname = f"{base_name}_{ccd}{suffix}.fits"
        if not os.path.exists(fname):
            raise FileNotFoundError(f"Missing file: {fname}")

        with fits.open(fname) as h:
            data = h[0].data
            header = h[0].header

            # Clean up header keywords that shouldn't propagate
            header.pop('SIMPLE', None)
            header.pop('BITPIX', None)
            header.pop('NAXIS', None)
            for k in list(header.keys()):
                if k.startswith('NAXIS'):
                    header.pop(k, None)

            # Add CCD identification (SCAMP-friendly)
            header['CCDID'] = ccd
            header['EXTNAME'] = f'CCD{ccd}'

            img_hdu = fits.ImageHDU(data=data, header=header)
            hdulist.append(img_hdu)

    outname = f"{base_name}_MEF.fits"
    hdulist.writeto(outname, overwrite=True)
    print(f"Wrote {outname}")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ="This program reconstructs the MEF files from the split images.  Hopefully this improves the performace of scamp for some challenging cases...")
    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC, which will grab all WFC*PA.fits.')
    parser.add_argument('--verbose', dest = 'verbose', default = False, action='store_true', help = 'Print more information')        
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'Will run on one directory only.')
    
    args = parser.parse_args()    

    # glob all WFC*PA.fits (or mWFC*PA.fits)
    flist = glob.glob(f'{args.filestring}*PA.fits')
    flist.sort()
    if args.verbose:
        print("flist = ",flist)
    # look for unique set of all image prefixes among all the *1PA.fits *2PA.fits, etc files
    allfilerootlist = [t.split('_')[0] for t in flist]
    filerootlist = list(set(allfilerootlist))

    for froot in filerootlist:
        build_mef_from_ccds(base_name, ccd_ids=(1, 2, 3, 4), suffix="PA")

        
        if args.testing:
            sys.exit()

