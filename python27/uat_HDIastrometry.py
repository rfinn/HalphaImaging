#!/usr/bin/env python

import glob
import os
import subprocess
#from astropy.io import fits


files = sorted(glob.glob('ftr*o00.fits'))
for f in files:
    print 'RUNNING ASTROMETRY FOR ',f
  #  os.system('gethead '+f+' EXPTIME > timecheckfile')
  #  infile=open('timecheckfile','r')
  #  exposuretime = []
  #  for line in infile:
  #      t=line.split()
  #      exposuretime.append(t[0])
  #  infile.close()
    read_exptime = 'gethead ' + f + ' EXPTIME'
    exptime = subprocess.check_output(read_exptime,shell=True)
    exptime = exptime.rstrip()
    print '   EXPOSURE TIME = ',exposuretime[0]
  #  os.remove('timecheckfile')
    if float(exptime) > 60.:
        t = f.split('.fits')
        froot = t[0]
        os.system('/usr/local/astrometry/bin/solve-field --cpulimit 30 ' + f)
        os.rename( froot + '.new' , 'w' + froot + '.fits')
    else:
        pass
    
if not(os.path.exists('astrometry')):
    os.mkdir('astrometry')

try:
    os.system('mv *.axy astrometry/.')
    os.system('mv *.corr astrometry/.')
    os.system('mv *.xyls astrometry/.')
    os.system('mv *.match astrometry/.')
    os.system('mv *.png astrometry/.')
    os.system('mv *.rdls astrometry/.')
    os.system('mv *.solved astrometry/.')
    os.system('mv *.wcs astrometry/.')
except:
    print 'WARNING:  No astrometry.net files to move.'
