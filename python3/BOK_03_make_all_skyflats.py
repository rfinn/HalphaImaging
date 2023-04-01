#!/usr/bin/env python

"""
wrapper function for BOK_03_sky_flat.py

going to try to use subprocess
"""
from subprocess import Popen
import os
homedir = os.getenv("HOME")
filters =['r']
print(homedir)
codedir = os.path.join(homedir,'github/HalphaImaging/python3/')
print(codedir)
for f in filters:
    for i in range(16):
        
        cmd_str = 'python '+codedir+f'BOK_03_sky_flat.py {f} {i+1}'
        #proc = Popen([cmd_str], shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)
        os.system(cmd_str)

