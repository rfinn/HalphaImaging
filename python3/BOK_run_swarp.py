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
            # prepend the m to match the name of the median subtracted image
            if os.path.exists('m'+f):
                os.rename('combined_weight.fits','m'+combined_mask)
            else:
                os.rename('combined_weight.fits',combined_mask)
def run_swarp(image_list,refimage=None):
    sstring = 'swarp @{} --WEIGHT_IMAGE {} --WEIGHT_SUFFIX .combweight.fits --COMBINE_TYPE WEIGHTED '.format(image_list,weight_list)
    if refimage is not None:
        # copying this from uat_astr_mosaic.py
        # still need to fix this.
        data,header = fits.getdata(args.refimage,header=True)
        w = WCS(header)
        image_size = data.shape

        ra,dec = w.wcs_pix2world(image_size[0]/2.,image_size[1]/2.,1)
        center = str(ra)+','+str(dec)
        mosaic_image_size = str(image_size[1])+','+str(image_size[0])
        
        sstring += ' --refimage {} '.format(refimage)
        
        commandstring = 'swarp @' + args.l + ' -c '+defaultswarp+' -IMAGEOUT_NAME ' + args.l + outimage+' -WEIGHTOUT_NAME ' + args.l + weightimage+' -CENTER_TYPE MANUAL -CENTER '+center+' -PIXEL_SCALE '+str(pixel_scale)+' -IMAGE_SIZE '+mosaic_image_size         

        
    os.system(sstring)
    
def run_swarp_all(image_list):
    # run swarp on r-band mosaic

    # run swarp on Halpha, using r-band mosaic as ref image

    # run swarp on r-band, using r-band mosaic as ref image
    pass

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



if __name__ == '__main__':
    telescope = 'BOK'
    
    parser = argparse.ArgumentParser(description ='stack the 90prime images after noao pipeline')
    parser.add_argument('--filestring', dest = 'filestring', default = 'ksb', help = 'filestring to match. default is ksb')
        
    parser.add_argument('--submedian', dest = 'submedian', default = False, action='store_true',help = 'set this to subtract the median from images.')
    parser.add_argument('--combinemasks', dest = 'combinemasks', default = False, action='store_true',help = 'set this to combine weight image and bad pixel mask.')        
    args = parser.parse_args()

        
    os.system('gethead object exptime FILTER RA DEC '+args.filstring+'*ooi*.fits > header_info')
    t = Table.read('header_info',data_start=0,delimiter=' ',format='ascii',guess=False,fast_reader=False,names=['FILENAME','OBJECT','EXPTIME','FILTER','RA','DEC'])

    

    # sort images by location and filter
    # alternatively could use object name,
    # but not all are correct, so need to fix names

    # get a list of all the unique objects
    # this is like, e.g. VFID2911_r or VFID2911_Ha4
    targets = sort(list(set(t['OBJECT'])))
    write_filelists(targets,t)

    # get list of r-band objects only
    primary_targets = []
    for t in targets:
        if t.ends_with('_r'):
            primary_targets.append(t)
    print('{} primary targets'.format(len(primary_targets)))
    print(primary_targets)

    if args.submedian:
        # subtract median
        os.system('python ~/github/HalphaImaging/python3/subtract_median.py --filestring {} --filestring2 {} --mef '.format(args.filestring,'ooi'))

    
    if args.combinemasks:
        # combine masks
        # this combines weight image and bad pixel masks
        combine_all_masks(t['FILENAME'])
    

    # subtract median from sky
    
    # run swarp

    # run swarp to mosaic r-band

    # run swarp to mosaic halpha
    
    # alight r and halpha imaging
    

