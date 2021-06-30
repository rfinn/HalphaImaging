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


#########

problem running median subtract on  ksb_210315_104538_ooi_r_v1.fits - moving this file to temp, and continuing.

updated subtract_median_sky to calculate median using astropy.stats.sigma_clipped_stats when the first attempt
returns a median that is == nan.

added file and median-subtracted file back into main directory

#####################

'''

import os

import argparse

from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time

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
        #combined_mask = 'm'+combined_mask
        if os.path.exists(combined_mask):
            print('output image already exists: ',combined_mask)
            print('proceeding to next image')
            print()
            continue
        else:
            weight_image = f.replace('ooi','oow')
            dq_image = f.replace('ooi','ood')
            combine_masks(weight_image,dq_image)
            # prepend the m to match the name of the median subtracted image
            os.rename('combined_weight.fits',combined_mask)
def run_swarp(image_list,refimage=None):
    '''

    RETURNS:
    * name of output image from swarp
    '''
    print(image_list)
    vfid,filter = image_list.split('_')
    weight_list = image_list+'_weights'

    # get date of observation from the first image in the list
    images = open(image_list,'r').readline()
    dateobs = images.split('_')[1]
    dateobs = '20'+dateobs
    
    output_image = 'VF-{}-BOK-{}-{}.fits'.format(dateobs,vfid,filter)
    # start building swarp command
    commandstring = 'swarp @{} -WEIGHT_IMAGE @{} -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME {}'.format(image_list,weight_list,output_image)
    if refimage is not None:
        # copying this from uat_astr_mosaic.py
        # still need to fix this.
        data,header = fits.getdata(refimage,header=True)
        w = WCS(header)
        image_size = data.shape

        ra,dec = w.wcs_pix2world(int(image_size[0]/2.),int(image_size[1]/2.),1)
        center = str(ra)+','+str(dec)
        mosaic_image_size = str(image_size[1])+','+str(image_size[0])
        
        commandstring = commandstring + ' -CENTER_TYPE MANUAL -CENTER {} -PIXEL_SCALE {} -IMAGE_SIZE {} '.format(center,pixel_scale,mosaic_image_size)

        
    os.system(commandstring)
    return output_image
    
def run_swarp_all_filters(target):
    '''
    INPUT:
    * image_list : list containing r-band images, like VFID0422_r

    PROCEDURE:
    * run swarp on r-band image, then run on Halpha using r as reference, then rerun on r using r as reference
    '''
    # run swarp on r-band mosaic
    rfilelist = target
    rband_coadd = run_swarp(rfilelist)
    
    # run swarp on Halpha, using r-band mosaic as ref image
    hafilelist = target.replace('_r','_Ha4')
    ha_coadd = run_swarp(hafilelist,refimage=rband_coadd)
    # run swarp on r-band, using r-band mosaic as ref image
    temp = run_swarp(rfilelist,refimage=rband_coadd)    
    pass

def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i+1
    except UnboundLocalError:
        return 0


def write_filelists(targets,header_table,medsub=False):
    for t in targets:
        # open file to store the list of science images
        outfile = open(t,'w')
        # open file to store the list of weight images
        weightfile = open(t+'_weights','w')
        # keep the filenames that match the current target name
        filenames = header_table['FILENAME'][header_table['OBJECT'] == t]
        # loop over the filenames and write each science and
        # weight image to the corresponding list
        for f in filenames:
            if medsub:
                outfile.write('m{} \n'.format(f))
            else:
                outfile.write('{} \n'.format(f))
            combined_mask = f.replace('.fits','.combweight.fits')
            weightfile.write('{} \n'.format(combined_mask))
        outfile.close()
        weightfile.close()



if __name__ == '__main__':
    telescope = 'BOK'
    
    parser = argparse.ArgumentParser(description ='stack the 90prime images after noao pipeline')
    parser.add_argument('--filestring', dest = 'filestring', default = 'ksb', help = 'filestring to match. default is ksb')
        
    parser.add_argument('--submedian', dest = 'submedian', default = False, action='store_true',help = 'set this to subtract the median from images.')
    parser.add_argument('--combinemasks', dest = 'combinemasks', default = False, action='store_true',help = 'set this to combine weight image and bad pixel mask.')
    parser.add_argument('--sortfiles', dest = 'sortfiles', default = False, action='store_true',help = 'write image and weights to files')
    parser.add_argument('--swarp', dest = 'swarp', default = False, action='store_true',help = 'run swarp to create coadded images')                
    args = parser.parse_args()

        
    os.system('gethead -a object exptime FILTER RA DEC '+args.filestring+'*ooi*v1.fits > header_info')
    filetable = Table.read('header_info',data_start=0,delimiter=' ',format='ascii',guess=False,fast_reader=False,names=['FILENAME','OBJECT','EXPTIME','FILTER','RA','DEC'])

    

    # sort images by location and filter
    # alternatively could use object name,
    # but not all are correct, so need to fix names

    # get a list of all the unique objects
    # this is like, e.g. VFID2911_r or VFID2911_Ha4
    targets = list(set(filetable['OBJECT']))
    targets.sort()

    # get list of r-band objects only
    primary_targets = []
    for t in targets:
        if t.endswith('_r'):
            primary_targets.append(t)
    print('{} primary targets'.format(len(primary_targets)))
    print(primary_targets)

    
    # subtract median from sky
    if args.submedian:
        # subtract median
        os.system('python ~/github/HalphaImaging/python3/subtract_median.py --filestring {} --filestring2 {} --mef '.format(args.filestring,'ooi_r_v1.fits'))
        #os.system('python ~/github/HalphaImaging/python3/subtract_median.py --filestring {} --filestring2 {} --mef '.format(args.filestring,'ooi_Ha+4nm_v1.fits'))        

    
    if args.combinemasks:
        # combine masks
        # this combines weight image and bad pixel masks
        combine_all_masks(filetable['FILENAME'])
    
    #print(targets)
    # need to update to write median-subtracted images to filelist instead of ksb files
    if args.sortfiles:
        write_filelists(targets,filetable,medsub=True)


    if args.swarp:
        for target in primary_targets:
            run_swarp_all_filters(target)

    

