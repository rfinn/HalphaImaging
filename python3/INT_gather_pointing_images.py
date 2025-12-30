#!/usr/bin/env python

'''

PROCEDURE:
* get list of current folders in the directory
* combine all filter images from the same pointing into one directory
* I had already split by pointing and filter, but I want to run scamp and swarp in the directory that contains all images from a pointing


'''

import os
import shutil
from astropy.io import fits
import glob
import argparse

if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description ="This program will combine images from e.g. pointing001-r and pointing001-Halpha into pointing001 and remove the filter-based directories.")
    #parser.add_argument('--filestring', dest = 'filestring', default = 'coadd', help = 'filestring to match. default is coadd, which will grab all *coadd*.fits.  Note: weight images will be renamed with their corresponding science coadd.')
        
    #parser.add_argument('--indir', dest = 'indir', default = '.', help = 'directory of input images.  Default is current directory')
    #parser.add_argument('--outdir', dest = 'outdir', default = '.', help = 'directory to write output images to.  Default is current directory')
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'Commands are printed when this is set, but the files are not moved.  This is useful when testing changes to the code.')    

    
    args = parser.parse_args()
    # get list of current directory
    flist1 = os.listdir()
    flist1.sort()

    # get list of pointings - cut off filter from dir name
    pointings = []
    for f in flist1:
        if( (f.find('pointing') > -1) | (f.find('target') > -1)) & (f.find('-') > -1):
            pointings.append(f.split('-')[0])

    # pointings will now have two listing for each filter b/c of r and Halpha
    upointings = set(pointings) # unique pointings
    print(f"found {len(upointings)} unique pointings")
    
    subdirs = ['plots','shortexposure']
    hafilters = ['Halpha','Ha6657']
    for p in upointings:
        print(f"working on {p}")
        if not(os.path.exists(p)):
            if args.testing:
                print(f"os.mkdir({p})")
            else:
                os.mkdir(p)
        if os.path.exists(p+'-r'):
            if args.testing:
                print('mv '+p+'-r/* '+p+'/.')
                print('rmdir '+p+'-r')
            else:
                os.system('mv '+p+'-r/* '+p+'/.')
                os.system('rmdir '+p+'-r')
        for h in hafilters:
            hadir = p+'-'+h
            if os.path.exists(hadir):
                if args.testing:
                    print('mv {}/* {}/.'.format(hadir,p))
                else:
                    os.system('mv {}/* {}/.'.format(hadir,p))
                for s in subdirs:
                    subdirpath = '{}-{}/{}'.format(p,h,s)
                    if os.path.exists(subdirpath):
                        command = 'mv {}/* {}/{}/.'.format(subdirpath,p,s)
                        if args.testing:
                            print(command)
                        else:
                            os.system(command)
                    if args.testing:
                        print('rmdir '+subdirpath)
                    else:
                        os.system('rmdir '+subdirpath)
                if args.testing:
                    print('rmdir {}'.format(hadir))
                else:
                    os.system('rmdir {}'.format(hadir))




