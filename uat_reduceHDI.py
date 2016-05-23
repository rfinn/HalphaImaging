#!/usr/bin/env python
'''
  GOAL:
  - to reduce HDI data through the flatfielding step

  PROCEDURE:
  - user should move junk files to a junk subdirectory before starting
    - junk files include initial bias frames, pointing, and other garbage frames
  - use gethead to pull relevant header data
  - overscan subtract and trim images
  - combine flats according to flat type (dome vs sky) and filter
  - normalize flat
  - divide object frames by normalized flats
    - initially I am going to make one set of images using dome flats, and one set using sky flats
  
  
  USEAGE:

  REQUIRED MODULES:
  - pyraf

  NOTES:
  - as of this writing (http://stupendous.rit.edu/richmond/wiyn/hdi_oct2013/tech_7/tech_7.html)
    - bias frames should not be used
    - dark current does not need to subtracted (need to verify this)
  - developed to reduce HDI data from Apr 2015
  
  WRITTEN BY:
  Rose A. Finn, 5/27/15

  UPDATES:

'''
from pyraf import iraf

# subtract overscan bias level
ccdproc(images=inlist,output=outlist,ccdtype='',overscan=True,trim=True,fixpix=False,zerocor=False,darkcor=False,flatcor=False,biassec='[4100:4140,:]',trimsec='[1:4112,1:4096]')
# combine flats

# normalize flats

# divide images by flats

# combine sky flats

# normalize sky flats

# divide images by sky flats
