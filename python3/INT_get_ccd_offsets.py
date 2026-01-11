#!/usr/bin/env python
"""
OVERVIEW:
This program compares each ccd from an image with panstarrs to determine a relative photometric offset.  Without this step, we get clear systematics between ccds when running getzp.py.

This can be run after median subtraction - just need to specify the file prefix.

PROCEDURE:
* get list of images
* get the matching SE catalogs

"""

import numpy as np
import os
import sys
from astropy.io import fits


def panstarrs_query(ra_deg, dec_deg, rad_deg, maxmag=19,
                    maxsources=10000):
    """
    FOUND THIS ON THE WEB 
    https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/

    
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object

    https://vizier.cds.unistra.fr/viz-bin/VizieR-2
    """
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    from astropy import units as u
    
    pan_columns =['objID', 'RAJ2000', 'DEJ2000','e_RAJ2000', 'e_DEJ2000', 'f_objID', 'Qual','gmag', 'e_gmag','rmag', 'e_rmag','imag', 'e_imag','zmag', 'e_zmag','ymag', 'e_ymag']
    #print(pan_columns)
    
    vquery = Vizier(columns=pan_columns,column_filters={"gmag":("<%f" % maxmag)},row_limit=maxsources)
    print('HEY!!! in panstarrs_query, ra,dec,rad = ',ra_deg,dec_deg,rad_deg)
    field = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg),frame='icrs')
    t = vquery.query_region(field, width=("%fd" % rad_deg), catalog="II/349/ps1")
    #print("in panstarrs_query...")
    #print(t)
    return t[0]


def get_panstarrs(rootname, centerRA, centerDEC, radius_deg):
    """
    PARAMS:
    * rootname - to use in constructing the panstarrs catalog name
    * centerRA - RA in deg
    * centerDEC - DEC in deg
    * radius_deg - radius of search area in deg

    RETURN:
    * pan - panstarrs table

    """
    from astropy.table import Table
    import os

    ###################################
    # get Pan-STARRS catalog over the same region
    ###################################
    ptab_name = f"{rootname}_pan_tab.csv" # don't see how this works - there is no fits in the basename

    #print(f"image={image}, panstarrs table = {ptab_name}")
    if os.path.exists(ptab_name):
        print('panstarrs table already downloaded')
        ptab = Table.read(ptab_name)

    else:
        print()
        print("no previous panstarrs catalog found :(")
        pan = panstarrs_query(centerRA, centerDEC, radius_deg)
        ptab = Table(pan)
        ptab.write(ptab_name,format='csv',overwrite=True)
    #color_correct_panstarrs()
    return ptab
        
def read_se_cat(se_cat):
    '''
    PARAMS:
    se_cat is an hdu of a se catalog (FITS_LDAC) with data in hdu 2

    keep only flag=0 sources and stars
    
    RETURN
    data table with good sources
    '''
    from astropy.table import Table
    htab = fits.getdata(se_cat,2)
    keepflag =  (htab['FLAGS'] <  1) #& (htab['CLASS_STAR'] > 0.95) #& (htab['MAG_AUTO'] > -11.)
    #print("Hello!")
    print(f"Number of sources in SE catalog = {np.sum(keepflag)} {se_cat}")
    return htab[keepflag]

def get_ra_dec_range(table_list):
    """
    PARAMS:
    * table_list: list of SE tables (hdu 2 in SE FITS_LDAC catalog)

    RETURNS:
    * racenter - central ra position of all the tables (deg)
    * deccenter - central dec position of all the tables (deg)
    * radius_deg - radius (deg) needed to enclose all the data points #, accounts for cos(dec) effect on delta_ra
    """
    ramin=10000
    ramax=0
    decmin=10000
    decmax=0

    for t in table_list:
        #print(t)
        if np.min(t['ALPHA_J2000']) < ramin:
            ramin = np.min(t['ALPHA_J2000'])
        if np.max(t['ALPHA_J2000']) > ramax:
            ramax = np.max(t['ALPHA_J2000'])
        if np.min(t['DELTA_J2000']) < decmin:
            decmin = np.min(t['DELTA_J2000'])
        if np.max(t['DELTA_J2000']) > decmax:
            decmax = np.max(t['DELTA_J2000'])
    racenter = 0.5*(ramax+ramin)
    deccenter = 0.5*(decmax+decmin)    
    raradius = 0.5*(ramax-ramin)#*np.cos(np.radians(deccenter)) # I don't think we need this?
    decradius = 0.5*(decmax-decmin)
    radius_deg = np.sqrt(decradius**2 + raradius**2)
    return racenter, deccenter, radius_deg

def cross_match_tables(table_list,rakey0='RAJ2000',deckey0='DEJ2000', rakey1='ALPHA_J2000',deckey1='DELTA_J2000'):
    '''
    DESCRIPTION:
    * table_list - list of tables

    PARAMS:
    * table_list: 

    RETURN:
    * matched_tables
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    # use ra and dec of first table (panstarrs) as the ref
    cref = SkyCoord(ra = table_list[0][rakey0],dec=table_list[0][deckey0], \
                        unit=(u.degree,u.degree),frame='icrs')

    allmatchflag = np.ones(len(table_list[0]),'bool')
    allindices = [np.arange(len(table_list[0]))]
    matched_tables = [table_list[0]]
    for i in range(1,len(table_list)):
        cmatch = SkyCoord(ra=table_list[i][rakey1],\
                              dec=table_list[i][deckey1],\
                              unit=(u.degree,u.degree),frame='icrs')

        index,dist2d,dist3d = cref.match_to_catalog_sky(cmatch)

        # only keep matches with matched RA and Dec w/in 3 arcsec
        matchflag = dist2d.degree < 5./3600
        print(f"the number of matches = {np.sum(matchflag)}")
        #print(len(cref), len(cmatch), len(index), len(dist2d), len(matchflag))

        matchedarray1=np.zeros(len(cref),dtype=table_list[i].dtype)
        matchedarray1[matchflag] = cmatch[index[matchflag]]
        print(f"median FLUX_AUTO of matched table = {np.median(matchedarray1['FLUX_AUTO'])}")
        matched_tables.append(matchedarray1)
        
        # keep track of sources that are matched in all catalog
        #allmatchflag = allmatchflag & matchflag

        # keep track of indices that point to the cmatch table
        allindices.append(index)

        # index of closest match in secoords cat
        #index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)
        #matchflag = dist2d.degree < 5./3600


        # REFERENCE CODE FROM GETZP.PY
        # index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)
        ## only keep matches with matched RA and Dec w/in 5 arcsec
        # matchflag = dist2d.degree < 5./3600
        # matchedarray1=np.zeros(len(pancoords),dtype=secat.dtype)
        # matchedarray1[matchflag] = secat[index[matchflag]]
    
    matched_tables = [tab[allmatchflag] for tab in matched_tables]
    if args.verbose:
        for t in matched_tables: print(len(t), np.sum(allmatchflag))
    return matched_tables

def match_tables_to_panstarrs(panstarrs_table,table_list,rakey0='RAJ2000',deckey0='DEJ2000', rakey1='ALPHA_J2000',deckey1='DELTA_J2000'):
    '''
    DESCRIPTION:
    * table_list - list of tables

    PARAMS:
    * table_list: 

    RETURN:
    * matched_tables
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.table import Table, Column
    
    # use ra and dec of first table (panstarrs) as the ref
    cref = SkyCoord(ra = panstarrs_table[rakey0],dec = panstarrs_table[deckey0], \
                        unit=(u.degree,u.degree),frame='icrs')

    allmatchflag = np.ones(len(table_list[0]),'bool')
    allindices = [np.arange(len(table_list[0]))]
    matched_tables = []
    pan_columns = ['gmag', 'e_gmag','rmag', 'e_rmag']
    for i in range(len(table_list)):
        secoords = SkyCoord(ra=table_list[i][rakey1],\
                              dec=table_list[i][deckey1],\
                              unit=(u.degree,u.degree),frame='icrs')

        index,dist2d,dist3d = secoords.match_to_catalog_sky(cref) # output is row matched to cmatch/se table

        # only keep matches with matched RA and Dec w/in 3 arcsec
        matchflag = dist2d.degree < 6./3600

        if np.sum(matchflag) < 1:
            print("WARNING: NO MATCHES BETWEEN SE AND PANSTARRS CATALOG!")
            print(f"\tdist2d mean={np.mean(dist2d.to('arcsec')):.2f}, median={np.median(dist2d.to('arcsec')):.2f}")
            print(f"\tdist2d min={np.min(dist2d.to('arcsec')):.2f}, max ={np.min(dist2d.to('arcsec')):.2f}")
        else:
            print(f"NUMBER OF MATCHES BETWEEN PANSTARRS AND SE = {np.sum(matchflag)}")
        outtab = Table(table_list[i])
        for c in pan_columns:
            newcol = Column(np.zeros(len(table_list[i]), 'd'), c)
            # add panstarrs g and r mags to tables
            newcol[matchflag] = panstarrs_table[c][index[matchflag]]
            outtab.add_column(newcol)
                                                   
                       
        # cut table to keep matches
        matched_tables.append(outtab[matchflag])
        #print(outtab.colnames)
        



        # REFERENCE CODE FROM GETZP.PY
        # index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)
        ## only keep matches with matched RA and Dec w/in 5 arcsec
        # matchflag = dist2d.degree < 5./3600
        # matchedarray1=np.zeros(len(pancoords),dtype=secat.dtype)
        # matchedarray1[matchflag] = secat[index[matchflag]]

    
    return matched_tables

def get_median_offsets(panstarrs_table, se_tabs, exptimes):
    '''
    match catalogs by ra and dec
    keep only flag=0 sources and stars
    
    '''
    
    # cross match arrays
    matched_tables = match_tables_to_panstarrs(panstarrs_table, se_tabs)
    medians = np.array([10**(np.nanmedian(tab['MAG_AUTO'] + 2.5*np.log10(exptime) - tab['rmag'])/2.5) for tab,exptime in zip(matched_tables,exptimes)])
 

    scale_factors = medians/np.mean(medians)
    return scale_factors             

def add_scale_to_header(image_list, scale_list, testing=False):
    from astropy.io import fits
    for i,im in enumerate(image_list):
        if not testing:
            hdu = fits.open(im)
            hdu[0].header.set('PANFZP',f"{scale_list[i]:.4f}",'Pan-STARRS photometric zeropoint scale (counts/s from FLUX_AUTO/EXPTIME -> PS1 flux)')
            hdu.writeto(im, overwrite=True)
        else:
            print(f'PANFZP, {scale_list[i]:.4f}, ','relative scale from PANSTARRS')

if __name__ == '__main__':
    import glob
    import argparse
    import os

    parser = argparse.ArgumentParser(description ="This program measures the offset of each ccd relative to PanSTARRS, calculates the scaling required to bring the ccds onto a uniform scale, and adds a keyword PANSCALE to the image headers.  PANSCALE can be used as the flux scaling keyword with SWarp.")
    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC, which will grab all WFC*PA.cat.')
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

    # sort list of unique image names
    filerootlist.sort()

    print("\n################################################")
    print(f"Working on directory {os.getcwd()}")
    print("################################################")        
    
    # why not get all images in a particular filter???
    filerootlist = [f"{args.filestring}.{f}." for f in ['r','Halpha','Ha6657']]
    # loop through list of image rootnames
    for i in range(len(filerootlist)):

        # get all files with the same root
        image_ccd_list = [m for m in flist if filerootlist[i] in m]
        
        if args.verbose:
            print(f'FYI found {len(image_ccd_list)} matching images for {filerootlist[i]}')

        if len(image_ccd_list) < 2:
            print(f"\tno matches to {filerootlist[i]}, moving to next filter.")
            continue

        # get exposure times from the image list
        exptimes = []
        for cimage in range(len(image_ccd_list)):
            header = fits.getheader(cimage)
            try:
                exptimes.append(float(header['EXPTIME']))
            except KeyError:
                print(f"WARNING: could not get exptime in header for {cimage}")
                sys.exit()

        # get se cats and measure offsets
        # might also need to remove the preceding m if running on mWFC images
        se_cat_list = [s.replace('.fits','.cat') for s in image_ccd_list]        
        #if 'mWFC' in args.filestring:
        #    se_cat_list = [s.replace('.fits','.cat') for s in image_ccd_list]
        #else:
        #    se_cat_list = [s.replace('.fits','.cat') for s in image_ccd_list]

        se_cat_list.sort()

        # read in the SE fits tables
        se_tabs = [read_se_cat(s) for s in se_cat_list]
        
        # get range of coordinates
        racenter, deccenter, radius_deg = get_ra_dec_range(se_tabs)
        if args.verbose:
            print(f"return values from get_ra_dec_range: ra={racenter:.5f}, dec={deccenter:.5f}, rad={radius_deg:.5f}")
        # download panstarrs catalog
        panstarrs_table = get_panstarrs(filerootlist[i], racenter, deccenter, 2*radius_deg)
        
        # calculate relative scale from median of FLUX_AUTO
        scale_factor_list = get_median_offsets(panstarrs_table, se_tabs, exptimes)
        #if args.verbose:
        print(f"\tscale factors = {scale_factor_list}")
        
        # add the scale factors to the header for use with SWARP
        # then swarp can use this to scale the images when making a coadd.
        add_scale_to_header(image_ccd_list, scale_factor_list,testing=args.testing)

        if args.testing:
            sys.exit()


