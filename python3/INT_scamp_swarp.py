#!/usr/bin/env python

'''
PROCEDURE:
get list of directories
move into each directory
run scamp and swarp

NOTES:
This assumes that the files are already sorted by pointing, but not by filter.

USAGE:
Run this from, e.g. /home/rfinn/data/reduced/scratch-int-feb2019
- this directory has a subdirectory for each pointing
- the subdirectory contains both the r and Halpha image

'''

import os
import sys
import glob
import shutil
from astropy.io import fits
import argparse
import multiprocessing as mp

def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i+1
    except UnboundLocalError:
        return 0

def run_one_se(filename):
    #print(('RUNNING SEXTRACTOR ON FILE %i OF %i'%(i,nfiles)))
    t = filename.split('.fits')
    froot = t[0]
    # DONE:TODO - check what needs to be updated in default.sex.INT - checked this an it's all ok
    os.system('sex ' + filename + ' -c default.sex.INT -CATALOG_NAME ' + froot + '.cat')
    
def get_updated_coadd_name(imname):
    f = imname
    h = fits.getheader(f)
    ra = float(h['CRVAL1'])
    dec = float(h['CRVAL2'])

    t,dateobs,telescope,pointing,filterwsuffix = f.split('-')
    # create string for output name

    if float(dec) < 0:
        outfile = 'VF-{:07.3f}-{:06.3f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,vfid,filterwsuffix)
    else:
        outfile = 'VF-{:07.3f}+{:06.3f}-{:s}-{:s}-{:s}-{:s}'.format(ra,abs(dec),telescope,dateobs,pointing,filterwsuffix)
    return outfile

    
if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description ="This program does most of the post-theli processing for the INT data.  Fingers crossed!")
    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC, which will grab all WFC*PA.fits.')
    parser.add_argument('--se', dest = 'se', default = False, action='store_true', help = 'Run source extractor')    
    parser.add_argument('--scamp', dest = 'scamp', default = False, action='store_true', help = 'Run scamp')
    parser.add_argument('--submedian', dest = 'submedian', default = False, action='store_true', help = 'Subtract a median sky value from each image.')
    parser.add_argument('--swarp', dest = 'swarp', default = False, action='store_true', help = 'Run swarp')    
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'Will run on one directory only.')
    
    args = parser.parse_args()

    
    homedir = os.getenv("HOME")
    telescope = 'INT'
    # get list of current directory
    flist1 = os.listdir()
    working_dir = os.getcwd()
    # overwrite output files if they exist
    overwrite = True
    flist1.sort()

    #flist1 = ['pointing022','pointing026']
    # setting this to false just b/c I am redoing mosaics for these two pointings
    submedian=False
    for subdir in flist1: # loop through list
        if os.path.isdir(subdir) & ((subdir.find('pointing') > -1) | ('target_' in subdir)):# & (subdir.find('-') > -1):
            print('##########################################')
            print('##########################################')        
            print('WORKING ON DIRECTORY: ',subdir)
            print('##########################################')
            print('##########################################')

            # move to subdirectory
            os.chdir(subdir)

            if args.se: # run source extractor
                os.system('cp ~/github/HalphaImaging/astromatic/default.* .')        
                filelist = glob.glob(args.filestring+'*PA.fits')
                filelist.sort()
                # link the astromatic files
        
                print(f"found {len(filelist)} files to run source extractor on")
                se_pool = mp.Pool(mp.cpu_count())
                seresults = [se_pool.apply_async(run_one_se,args=(filename,),callback=se_collect_results) for filename in filelist]
        
                se_pool.close()
                se_pool.join()
                se_results = [r.get() for r in seresults]
                
            # run scamp
            if args.scamp:
                scampflag=False # keeps track if scamp finishes successfully
                try:
                    os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --scamp --int --filestring WFC')
                    scampflag=True
                    os.system('ls WFC.r*PA.fits > '+subdir+'_r')
                    os.system('ls WFC.Halpha*PA.fits > '+subdir+'_Halpha')
                    os.system('ls WFC.Ha6657*PA.fits > '+subdir+'_Ha6657')            
                except:
                    print('##########################################')
                    print('WARNING: problem running scamp for ',subdir)
                    print('##########################################')

            if args.moveshort:
                # move short exposure times
                os.system('python ~/github/HalphaImaging/python3/move_short_exposures.py --filestring WFC')        

            if args.submedian:
                # subtract median
                os.system('python ~/github/HalphaImaging/python3/subtract_median.py --filestring WFC')
                
            # run swarp
            if args.swarp:

                # gather files - why am I not getting the mWFC files?
                os.system('ls WFC.r*PA.fits > '+subdir+'_r')
                os.system('ls WFC.Halpha*PA.fits > '+subdir+'_Halpha')
                os.system('ls WFC.Ha6657*PA.fits > '+subdir+'_Ha6657')            

                # count lines r band file, run if more than 2 lines
                suffix = ['_r','_Halpha','_Ha6657']
                filelists = [subdir+i for i in suffix]
                #################################################3
                # make r-band mosaic
                #################################################3            
                nlines = count_lines(filelists[0])

                if nlines > 2:
                    os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[0])
                    refimage = filelists[0]+'.coadd.fits'
                    for i in [1,2]:
                        nlines = count_lines(filelists[i])
                        if nlines > 3:
                            #os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[i]+' --refimage '+refimage)
                            os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[i]+' --refimage '+refimage+' --noback')
                    ## don't need to rebuilt the r-band mosaic using the r-band mosaic as reference now that I am using the swarp flags correctly
                    pass
                    # run swarp again on the rband data, using the same refimage
                    #os.system('cp '+refimage+' refimage.fits')
                    #os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[0]+' --refimage refimage.fits')
                    #os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+filelists[0]+' --refimage refimage.fits --noback')
                else:
                    print('WARNING: No r-band mosaic, making remaining mosaics without alignment')
                    for f in filelists:
                        nlines = count_lines(f)
                        if nlines > 2:
                            #os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+f)
                            os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --swarp --int --l '+f+' --noback')
                        else:
                            print('WARNING: not enough images to make mosaic in ',f)
                            
                ## not removing median subtracted images - these take so long to make!
                # remove median subtracted images
                #os.system('rm mWFC*.fits')

                ## RENAME COADDS
                new_output_image = get_updated_BOK_coadd_name(ha_coadd)
                print('renaming ',ha_coadd,'->',new_output_image)
                os.rename(ha_coadd,new_output_image)
                # rename the weight file
                os.rename(ha_coadd.replace('.fits','.weight.fits'),new_output_image.replace('.fits','.weight.fits'))


                newname = get_updated_BOK_coadd_name(rband_coadd)
                print('renaming ',rband_coadd,'->',newname)
                os.rename(rband_coadd,newname)
                # rename the weight file    
                os.rename(rband_coadd.replace('.fits','.weight.fits'),newname.replace('.fits','.weight.fits'))

    
            os.chdir(working_dir)
            if args.testing:
                # just running on one directory for testing purposes
                sys.exit()
