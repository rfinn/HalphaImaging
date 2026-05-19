#!/usr/bin/env python

from pathlib import Path
import re

COADD_DIR = Path("/data-pool/Halpha/coadds-v20260518")
PSF_DIR = Path("/data-pool/Halpha/psf-images-v20260518")

DRY_RUN = True   # set False after checking output


# def fixed_name(name):
#     return re.sub(
#         r"(VF-\d{3}\.\d{3})([+-])(\d\.\d{3})(-INT-)",
#         r"\1\20\3\4",
#         name,
#     )

def fixed_name(name):
    return re.sub(
        r"(VF-\d{3}\.\d{3})([+-])(\d\.\d{3})(-INT-)",
        r"\g<1>\g<2>0\g<3>\g<4>",
        name,
    )
def rename_bad_dec_files(directory):
    #t1 = glob("*INT*r.fits")
    #t2 = glob("*INT*r-psf.fits")
    #t3 = glob("*INT*.png")
    
    for p in sorted(directory.glob("*INT*")):

        newname = fixed_name(p.name)

        if newname == p.name:
            continue

        #print(p.name,newname)        
        newpath = p.with_name(newname)

        if newpath.exists():
            print(f"SKIP exists: {newpath}")
            continue

        print(f"{p.name}  ->  {newname}")

        if not DRY_RUN:
            p.rename(newpath)


print("COADDS")
rename_bad_dec_files(COADD_DIR)

print("\nPSFS")
rename_bad_dec_files(PSF_DIR)

print("\nPSF plots")
rename_bad_dec_files(PSF_DIR / "plots")
