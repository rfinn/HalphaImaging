#!/usr/bin/env python

'''

ORGANIZING DATA
* working in /mnt/qnap_home/rfinn/Halpha/Bok

* moving short exposure time images to subdirectory junk

%run ~/github/HalphaImaging/python3/move_short_exposures.py --filestring ksb --exptime 31

###########

checking object names

gethead object exptime RA DEC *ooi*.fits > header_info

fixed a few issues

############

data from 04/15 is pretty crappy
- photometric zps are sometimes 2-3 mags lower
- VFID2911 has data from 04/16, so moving the 04/15 data to junk
- going to keep the others, but the depth is def not good

'''

import os
from astropy.io import fits
import argparse



def combine_masks(weight_image,dq_image):
    '''
    combine the weight and data quality image
    combined weight = weight_image + 1000*data_quality_image
    then use the following flag in swarp
    WEIGHT_THRESH 1000

    INPUT: 
    * weight_image : this is the weight image
    * dq_image : data quality image, with nonzero indicating bad pixels

    OUTPUT:
    * combined_weight.fits : this is the combined weight image
    
    '''
    weight_hdu = fits.open(weight_image)
    dq_hdu = fits.open(dq_image)
    
    # loop over the 4 images, extensions 1-4
    for i in range(1,5):
        print('combining image number ',i)
        weight_hdu[i].data = weight_hdu[i].data + 1000*dq_hdu[i].data
    weight_hdu.writeto('combined_weight.fits',overwrite=True)

def combine_all_masks(filelist):
    for f in filelist:

        combined_mask = f.replace('.fits','.combweight.fits')
        if os.path.exists(combined_mask):
            continue
        else:
            weight_image = f.replace('ooi','oow')
            dq_image = f.replace('ooi','ood')
            combine_masks(weight_image,dq_image)
            os.rename('combined_weight.fits',combined_mask)
def run_swarp(image_list, weight_list,refimage=None):
    sstring = 'swarp @{} --WEIGHT_IMAGE {} --WEIGHT_SUFFIX .combweight.fits --COMBINE_TYPE WEIGHTED '.format(image_list,weight_list)
    if refimage is not None:
        sstring += '--refimage {} '.format(refimage)
    else:
        
    os.system(sstring)
    


def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i+1
    except UnboundLocalError:
        return 0

def write_filelists(targets,header_table):
    for t in targets:
        outfile = open(t,'w')
        filenames = header_table['OBJECT'] == t
        for f in filenames:
            outfile.write('{} \n'.format(f))
        outfile.close

homedir = os.getenv("HOME")
telescope = 'INT'
# get list of current directory
flist1 = os.listdir()
working_dir = os.getcwd()
# overwrite output files if they exist
overwrite = True
flist1.sort()
runscamp=False
runswarp=True


rawdir = os.getcwd()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='stack the 90prime images after noao pipeline')
        
    parser.add_argument('--filestring', dest = 'filestring', default = 'ksb', help = 'filestring to match. default is ksb')
    args = parser.parse_args()

    os.system('gethead object exptime FILTER RA DEC '+args.filstring+'*ooi*.fits > header_info')
    t = Table.read('header_info',data_start=0,delimiter=' ',format='ascii',guess=False,fast_reader=False,names=['FILENAME','OBJECT','EXPTIME','FILTER','RA','DEC'])

    

    # sort images by location and filter
    # alternatively could use object name,
    # but not all are correct, so need to fix names
    targets = sort(list(set(t['OBJECT'])))
    
    
    write_filelists(targets,t)

    # combine masks
    combine_all_masks(t['FILENAME'])
    # get list of r-band objects only
    primary_targets = []
    for t in targets:
        if t.ends_with('_r'):
            primary_targets.append(t)

    # run swarp

    # run swarp to mosaic r-band

    # run swarp to mosaic halpha
    
    # alight r and halpha imaging
    

