trying to get swarp working correctly on the BOK images.

output from BOK_run_swarp.py:
(virgo) rfinn@draco:/data-pool/Halpha/processed/BOK2021-22$ less temp_output
66 primary targets

Calling run_swarp for rfilelist


In run_swarp, image_list =  VFID4796_r  refimage =  None

Running swarp with the following command:
 swarp @VFID4796_r -WEIGHT_IMAGE @VFID4796_r_weights -c default.swarp.BOK -IMAGEOUT_NAME VF-20210418-BOK-VFID4796-r.fits -WEIGHTOUT_NAME VF-20210418-BOK-VFID4796-r.weight.fits -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.4533 

rband_coadd =  VF-20210418-BOK-VFID4796-r.fits

Calling run_swarp for hafilelist


In run_swarp, image_list =  VFID4796_Ha4  refimage =  VF-20210418-BOK-VFID4796-r.fits
made it into refimage block
at end of ref image block, commandstring =  swarp @VFID4796_Ha4 -WEIGHT_IMAGE @VFID4796_Ha4_weights -c default.swarp.BOK -IMAGEOUT_NAME VF-20210418-BOK-VFID4796-Ha4.fits -WEIGHTOUT_NAME VF-20210418-BOK-VFID4796-Ha4.weight.fits -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.4533  -CENTER_TYPE MANUAL -CENTER 235.01158309529717,12.538089978185097 -IMAGE_SIZE 9186,8819 

Running swarp with the following command:
 swarp @VFID4796_Ha4 -WEIGHT_IMAGE @VFID4796_Ha4_weights -c default.swarp.BOK -IMAGEOUT_NAME VF-20210418-BOK-VFID4796-Ha4.fits -WEIGHTOUT_NAME VF-20210418-BOK-VFID4796-Ha4.weight.fits -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE 0.4533  -CENTER_TYPE MANUAL -CENTER 235.01158309529717,12.538089978185097 -IMAGE_SIZE 9186,8819 
renaming  VF-20210418-BOK-VFID4796-Ha4.fits -> VF-235.012+12.538-BOK-20210418-VFID4796-Ha4.fits
renaming  VF-20210418-BOK-VFID4796-r.fits -> VF-234.988+12.515-BOK-20210418-VFID4796-r.fits

Then checked header info of coadded images, and the r and halpha don't have the same centers.

(virgo) rfinn@draco:/data-pool/Halpha/processed/BOK2021-22$ gethead NAXIS1 NAXIS2 RA DEC CRVAL1 CRPIX1 CRVAL2 CRPIX2 CD1_1 VF-23*.fits
VF-234.988+12.515-BOK-20210418-VFID4796-r.fits          9186 8819   234.9877839081 4593.500000000 12.51504827686 4410.000000000 -1.259166666667E-04
VF-234.988+12.515-BOK-20210418-VFID4796-r.weight.fits   9186 8819   234.9877839081 4593.500000000 12.51504827686 4410.000000000 -1.259166666667E-04
VF-235.012+12.538-BOK-20210418-VFID4796-Ha4.fits        9186 8819   235.0115830953 4593.500000000 12.53808997819 4410.000000000 -1.259166666667E-04
VF-235.012+12.538-BOK-20210418-VFID4796-Ha4.weight.fits 9186 8819   235.0115830953 4593.500000000 12.53808997819 4410.000000000 -1.259166666667E-04


So swarp is setting the RA/DEC as expected, but it does not match the RA/DEC in the reference image, even though I am getting the ra/dec from the refimage...


GOING TO PAUSE AND TRY TO SOLVE THIS WITH BECKY.

POST NOTE JAN 2026:
AND WE DID SOLVE IT!  I WAS NOT USING ALL OF THE REQUIRED PARAMETERS THAT YOU NEED WHEN SPECIFYING THE CENTER AND IMAGE SIZE.
