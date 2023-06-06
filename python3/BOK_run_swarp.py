#!/usr/bin/env python

'''
USAGE:

# add subtract median, although at this point all files have had the median subtracted...

python ~/github/HalphaImaging/python3/BOK_run_swarp.py --filestring mksb --combinemasks

python ~/github/HalphaImaging/python3/BOK_run_swarp.py --filestring mksb --se

python ~/github/HalphaImaging/python3/BOK_run_swarp.py --filestring mksb --scamp

python ~/github/HalphaImaging/python3/BOK_run_swarp.py --filestring mksb --sortfiles

python ~/github/HalphaImaging/python3/BOK_run_swarp.py --filestring mksb --swarp



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

2023-06-02

* Reducing 2022 data on draco (new server)

* updating code to run in multiprocessing mode to take advantage of server

2023-06-05

* need to add options to run source extractor and scamp b/c some of pipeline images have bad wcs


'''

import os
import glob
import argparse
import multiprocessing as mp

from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time


mcombine_results = []
def mcombine_collect_results(result):

    global results
    mcombine_results.append(result)

swarp_results = []
def swarp_collect_results(result):

    global results
    swarp_results.append(result)

se_results = []
def se_collect_results(result):

    global results
    se_results.append(result)
    
getzp_results = []
def getzp_collect_results(result):

    global results
    getzp_results.append(result)

getzp_results2 = []
def getzp_collect_results2(result):

    global results
    getzp_results2.append(result)

def combine_masks(imname):
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
    

    2023-06-02: updating to use multiprocessing

    '''

    combined_mask = imname.replace('.fits','.combweight.fits')
    #combined_mask = 'm'+combined_mask
    if os.path.exists(combined_mask):
        print('output image already exists: ',combined_mask)
        print('proceeding to next image')
        print()
        return
    else:
        print("combining weights for ",imname)
        weight_image = imname.replace('ooi','oow')
        dq_image = imname.replace('ooi','ood')
        weight_hdu = fits.open(weight_image)
        dq_hdu = fits.open(dq_image)
    
        # loop over the 4 images, extensions 1-4
        for i in range(1,5):
            #print('combining image number ',i)
            weight_hdu[i].data = weight_hdu[i].data + 1000*dq_hdu[i].data
        weight_hdu.writeto(combined_mask,overwrite=True)
        
        # prepend the m to match the name of the median subtracted image
        #os.rename('combined_weight.fits',combined_mask)
    

def combine_all_masks(filelist):

    ##
    # implement multiprocessing
    ##

    mcombine_pool = mp.Pool(mp.cpu_count())
    myresults = [mcombine_pool.apply_async(combine_masks,args=(im,),callback=mcombine_collect_results) for im in filelist]
    
    mcombine_pool.close()
    mcombine_pool.join()
    mcombine_results = [r.get() for r in myresults]
    return mcombine_results

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

    # remove the nm from Ha4nm
    filter = filter.replace('nm','')
    
    os.system('cp ~/github/HalphaImaging/astromatic/default.swarp.BOK .')

    ##
    # add RA and DEC - 2023-05-10
    # didn't actually do this for the 2021 data, but changing now so that this is
    # implemented for 2022 data
    #
    # 2023-06-02
    # however, the RA and DEC of the first image is not necessarily going to be
    # the RA and DEC of the mosaiced image -
    # maybe this is why I didn't fully implement
    ##
    output_image = 'VF-{}-BOK-{}-{}.fits'.format(dateobs,vfid,filter)
    output_weight_image = 'VF-{}-BOK-{}-{}.weight.fits'.format(dateobs,vfid,filter)    
    # start building swarp command
    commandstring = 'swarp @{} -WEIGHT_IMAGE @{} -c default.swarp.BOK -IMAGEOUT_NAME {} -WEIGHTOUT_NAME {} '.format(image_list,weight_list,output_image,output_weight_image)
    
    if refimage is not None:
        # copying this from uat_astr_mosaic.py
        # still need to fix this.
        data,header = fits.getdata(refimage,header=True)
        w = WCS(header)
        image_size = data.shape
        pixel_scale = 0.453 # pixel scale for 90prime
        ra,dec = w.wcs_pix2world(int(image_size[0]/2.),int(image_size[1]/2.),1)
        center = str(ra)+','+str(dec)
        mosaic_image_size = str(image_size[1])+','+str(image_size[0])
        
        commandstring = commandstring + ' -CENTER_TYPE MANUAL -CENTER {} -PIXEL_SCALE {} -IMAGE_SIZE {} '.format(center,pixel_scale,mosaic_image_size)


    print('')
    print('Running swarp with the following command:\n',commandstring)
    os.system(commandstring)





    
    return output_image

def get_updated_BOK_coadd_name(imname):
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
def update_header(image,refimage):
    '''
    GOAL:
    scamp is not passing the header keywords through to the coadded image,
    so the goal of this routine is to add the keywords manually...

    INPUT:
    * image - coadded image to add header keywords too
    * refimage - image to get header keywords from
    '''
    header_fields = ['OBJECT','EXPTIME','FILTER','TELESCOP','INSTRUME','GAIN1','DATE-OBS','AIRMASS','MAGZERO','MAGZSIG','SEEING']# List of FITS keywords to propagate

    idata,iheader = fits.getdata(image,header=True)

    rheader = fits.getheader(refimage)

    for f in header_fields:
        try:
            iheader.set(f,rheader[f])
        except KeyError:
            print(f"WARNING: Keyword {f} not found")
    fits.writeto(image,idata,header=iheader,overwrite=True)

    
def run_swarp_all_filters(target):
    '''
    INPUT:
    * image_list : list containing r-band images, like VFID0422_r

    PROCEDURE:
    * run swarp on r-band image, 
    * then run on Halpha using r as reference, 
    * then rerun on r using r as reference

    '''
    # run swarp on r-band mosaic
    rfilelist = target
    rband_coadd = run_swarp(rfilelist)

    #try:
    # run swarp on Halpha, using r-band mosaic as ref image
    hafilelist = target.replace('_r','_Ha4')
    if not os.path.exists(hafilelist):
        hafilelist = target.replace('_r','_Ha4nm')
        if not os.path.exists(hafilelist):
            print("Warning - couldn't find Halpha images")
            return
    ha_coadd = run_swarp(hafilelist,refimage=rband_coadd)
    
    # run swarp on r-band, using r-band mosaic as ref image
    rband_coadd = run_swarp(rfilelist,refimage=rband_coadd)

    # update headers
    # but how does ref image have a 
    
    refimage = open(rfilelist).readline().rstrip()
    update_header(rband_coadd,refimage)
    
    refimage = open(hafilelist).readline().rstrip()
    update_header(ha_coadd,refimage)

    ##
    # rename r-band coadd using the updated naming convention
    # that includes RA and DEC
    ##
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

def count_lines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    try:
        return i+1
    except UnboundLocalError:
        return 0


def write_filelists(targets,header_table,medsub=True):
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
            #if medsub:
            #    outfile.write('m{} \n'.format(f))
            #else:
            outfile.write('{} \n'.format(f))
            combined_mask = f.replace('.fits','.combweight.fits').replace('mksb','ksb')
            weightfile.write('{} \n'.format(combined_mask))
        outfile.close()
        weightfile.close()


def getonezp(imname,filter):
    getzpstring = 'python ~/github/HalphaImaging/python3/getzp.py --image {} --instrument b --filter {} --normbyexptime'.format(imname,filter)
    os.system(getzpstring)


def run_one_se(filename):
    #print(('RUNNING SEXTRACTOR ON FILE %i OF %i'%(i,nfiles)))
    t = filename.split('.fits')
    froot = t[0]
    # DONE:TODO - check what needs to be updated in default.sex.INT - checked this an it's all ok
    os.system('sex ' + filename + ' -c default.sex.INT -CATALOG_NAME ' + froot + '.cat')
    





if __name__ == '__main__':
    telescope = 'BOK'
    
    parser = argparse.ArgumentParser(description ='stack the 90prime images after noao pipeline')
    parser.add_argument('--filestring', dest = 'filestring', default = 'ksb', help = 'filestring to match. default is ksb')
        
    parser.add_argument('--submedian', dest = 'submedian', default = False, action='store_true',help = 'set this to subtract the median from images.')
    parser.add_argument('--se', dest = 'se', default = False, action='store_true',help = 'run source extractor to create catalogs for scamp')
    parser.add_argument('--scamp', dest = 'scamp', default = False, action='store_true',help = 'run scamp to solve for WCS')
    parser.add_argument('--fixamps', dest = 'fixamps', default = False, action='store_true',help = 'fix ZP offsets between amplifiers in individual MEF images')    
    parser.add_argument('--combinemasks', dest = 'combinemasks', default = False, action='store_true',help = 'set this to combine weight image and bad pixel mask.')
    parser.add_argument('--sortfiles', dest = 'sortfiles', default = False, action='store_true',help = 'write image and weights to files')
    parser.add_argument('--swarp', dest = 'swarp', default = False, action='store_true',help = 'run swarp to create coadded images')
    parser.add_argument('--getzp', dest = 'getzp', default = False, action='store_true',help = 'run getzp to determine photometric zp of r and Halpha images')                    
    args = parser.parse_args()


    #print('command: gethead -a object exptime FILTER RA DEC '+args.filestring+'*ooi*v1.fits > header_info')
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
        os.system('python ~/github/HalphaImaging/python3/subtract_median.py --filestring {} --filestring2 {} --mef '.format(args.filestring,'ooi_Ha4nm_v1.fits'))
        
 
    if args.se:
        os.system('ln -s ~/github/HalphaImaging/astromatic/default.* .')        
        filelist = glob.glob('mksb*v1.fits')
        filelist.sort()
        # link the astromatic files
        
        print(f"found {len(filelist)} files to run source extractor on")
        se_pool = mp.Pool(mp.cpu_count())
        seresults = [se_pool.apply_async(run_one_se,args=(filename,),callback=se_collect_results) for filename in filelist]
    
        se_pool.close()
        se_pool.join()
        se_results = [r.get() for r in seresults]
       

    if args.scamp:
        os.system('ln -s ~/github/HalphaImaging/astromatic/default.* .')        
        os.system('ls '+args.filestring+'*.cat > scamp_input_cats')
        print('RUNNING SCAMP')
        # TODO - check to see what needs to be updated in default.scamp.INT
        os.system('scamp @scamp_input_cats -c default.scamp.BOK')
        pass

    if args.combinemasks:
        # combine masks
        # this combines weight image and bad pixel masks
        mcombine_results = combine_all_masks(filetable['FILENAME'])

    # TODO - add a function to fix ZP offsets in individual images.
    # like BOK_pipeline_fixampoffsets.py - but no median subtraction

    if args.fixamps:
        # call BOK_pipeline_fixampoffsets.py
        pass

    #print(targets)
    # need to update to write median-subtracted images to filelist instead of ksb files
    if args.sortfiles:
        write_filelists(targets,filetable,medsub=args.submedian)

        
    if args.swarp:
        #for target in primary_targets:
        #    run_swarp_all_filters(target)
            # break below is for debugging purposes
            #break
        
        swarp_pool = mp.Pool(mp.cpu_count())
        swresults = [swarp_pool.apply_async(run_swarp_all_filters,args=(target,),callback=swarp_collect_results) for target in primary_targets]
    
        swarp_pool.close()
        swarp_pool.join()
        swarp_results = [r.get() for r in swresults]

        
    if args.getzp:
        rfiles = glob.glob('VF*r.fits')
        getzp_pool = mp.Pool(mp.cpu_count())
        zpresults = [getzp_pool.apply_async(getonezp,args=(imname,'r'),callback=getzp_collect_results) for imname in rfiles]
    
        getzp_pool.close()
        getzp_pool.join()
        getzp_results = [r.get() for r in zpresults]
        
        #hfiles = glob.glob('VF*Ha4.fits')
        ##
        # 2022 data have Ha4nm
        hfiles = glob.glob('VF*Ha4nm.fits')
        getzp_pool2 = mp.Pool(mp.cpu_count())
        zpresults2 = [getzp_pool2.apply_async(getonezp,args=(imname,'ha'),callback=getzp_collect_results2) for imname in hfiles]
    
        getzp_pool2.close()
        getzp_pool2.join()
        getzp_results = [r.get() for r in zpresults2]


            
        
