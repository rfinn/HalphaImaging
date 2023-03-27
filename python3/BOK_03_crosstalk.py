#!/usr/bin/env python

import numpy as np
from astropy.io import fits

### from Zou+2017, https://iopscience.iop.org/article/10.3847/1538-3881/aa72d9
#Table 5 
#Crosstalk Coefficients for All Pairs of Amplifiers
#Note. The scale is 10^-5.

# "The crosstalk correction in a particular amplifier is computed by summing the pixel values in all other amplifiers multiplied by their crosstalk coefficients and then subtracted from this amplifier."


#HDU	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	
crosstalk_dict = {
    1: [0,	-23,	-33,	-30,	1,	2,	2,	2,	1,	-1,	1,	1,	2,	2,	1,	3],
    2: [-17,	0,	-24,	-23,	1,	1,	2,	3,	1,	0,	1,	1,	2,	3,	3,	3],
    3: [-10,	-8,	0,	-11,	1,	1,	2,	2,	1,	0,	1,	1,	2,	2,	2,	4],
    4: [-11,	-9,	-11,	0,	2,	-4,	2,	2,	1,	0,	0,	1,	1,	2,	2,	2],
    5: [3,	3,	3,	3,	0,	-33,	-37,	-27,	2,	-1,	2,	1,	0,	0,	1,	2],
    6: [3,	3,	4,	3,	-10,	0,	-17,	-14,	3,	2,	4,	3,	1,	1,	1,	2],
    7: [3,	3,	3,	3,	-11,	-14,	0,	-6,	2,	1,	2,	2,	1,	2,	2,	2],
    8: [3,	3,	4,	4,	-9,	-6,	-5,	0,	1,	3,	4,	2,	0,	1,	1,	1],
    9: [3,	2,	3,	2,	1,	1,	2,	2,	0,	-16,	-23,	-19,	3,	3,	2,	2],
    10: [2,	1,	2,	1,	0,	1,	1,	1,	-31,	0,	-40,	-32,	2,	2,	2,	2],
    11: [3,	2,	2,	1,	0,	-1,	1,	1,	-25,	-16,	0,	-17,	2,	2,	1,	2],
    12: [3,	2,	3,	2,	1,	1,	2,	2,	-15,	-11,	-8,	0,	-2,	1,	2,	1],
    13: [4,	3,	4,	3,	-1,	0,	1,	1,	-9,	0,	1,	3,	0,	-24,	-26,	-22],
    14: [4,	3,	4,	3,	0,	1,	1,	2,	-9,	1,	1,	-2,	-24,	0,	-26,	-22],
    15: [4,	4,	5,	4,	1,	1,	1,	2,	-11,	0,	1,	3,	-13,	-10,	0,	-5 ],
    16: [5,	5,	5,	5,	-1,	1,	1,	2,	-7,	-3,	-1,	0,	-15,	-14,	-12,	0  ]
    }


crosstalk_coeffs = np.zeros(16,16)
for i in range(16):
    coeffs = crosstalk_dict[i+1]
    for j,c in enumerate(coeffs):
        crosstalk_coeffs[i,j] = c
    
crosstalk_coeffs = crosstalk_coeffs*1.e-5

def correct_image(imname):
    hdu = fits.open(imname)
    crosstalk_correction_list = []
    for i in range(1,len(hdu)):
        crosstalk_correction = np.zeros(t[i].data.shape)
        crosstalk_correction
        for j in range(1,len(hdu)):
            crosstalk_correction += crosstalk_coeffs[j]*hdu[j].data
        crosstalk_correction_list.append(crosstalk_correction)

    # now apply the correction to the data

    
