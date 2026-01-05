#!/usr/env/bin python

"""
run getzp on all of the files
"""

if __name__ == '__main__':
    import argparse
    import glob
    import os
    import sys
    parser = argparse.ArgumentParser(description ="Run getzp on files in the current directory")
    parser.add_argument('--filestring', dest = 'filestring', default = 'UAT', help = 'filestring for coadds, default is UAT')
    parser.add_argument('--filter', dest = 'filter', default = None, help = 'add a filter name to run on one set of filters at a time.  could be useful if you want to run e.g. in two different windows.')        
    parser.add_argument('--testing', dest = 'testing', default = False, action='store_true', help = 'Will run on one directory only.')
    
    args = parser.parse_args()

    if args.filter is None:
        matchstring = f"{args.filestring}*.fits"
    else:
        matchstring = f"{args.filestring}*{args.filter}.fits"

    filelist = glob.glob(matchstring)
    filelist.sort()

    

    for f in filelist:
        if 'weight' in f:
            continue
        
        if 'MOS' in f:
            instrument = 'm'
        else:
            instrument = 'h'

        if 'R.fits' in f:
            filtername =  'R'
        elif 'r.fits' in f:
            filtername = 'r'
        else:
            filtername = 'h'
        # TODO - NEED TO ADD CASES FOR HA8, 12, AND 16 ONCE WE HAVE THOSE TRANSFORMATIONS
        print("running ",f"python ~/github/HalphaImaging/python3/getzp.py --image {f} --instrument {instrument} --filter {filtername}")
        os.system(f"python ~/github/HalphaImaging/python3/getzp.py --image {f} --instrument {instrument} --filter {filtername}")
        # move back to parent directory
        if args.testing:
            sys.exit()
