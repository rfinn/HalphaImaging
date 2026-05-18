#!/usr/bin/env python

import argparse
import subprocess
from pathlib import Path
from astropy.io import fits


def run(cmd):
    print(" ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True)


# def get_ref_params(reffile):
#     hdr = fits.getheader(reffile)

#     nx = hdr["NAXIS1"]
#     ny = hdr["NAXIS2"]
#     ra = hdr["CRVAL1"]
#     dec = hdr["CRVAL2"]

#     # assumes square-ish pixels; fine for SWarp target scale
#     pixscale = abs(hdr["CD1_1"]) * 3600.0 if "CD1_1" in hdr else abs(hdr["CDELT1"]) * 3600.0

#     return ra, dec, nx, ny, pixscale
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


def get_ref_params(reffile):
    data, hdr = fits.getdata(reffile, header=True)
    w = WCS(hdr)

    ny, nx = data.shape

    # center pixel in FITS/image coordinates
    xcen = nx / 2.0
    ycen = ny / 2.0

    # origin=0 because these are numpy/python pixel coordinates
    ra, dec = w.wcs_pix2world(xcen, ycen, 0)

    # pixel scale in arcsec/pixel
    pscale = proj_plane_pixel_scales(w)  # deg/pix
    pixscale = float(abs(pscale[0]) * 3600.0)

    return float(ra), float(dec), int(nx), int(ny), pixscale

def swarp_resample(infile, reffile, outname, weightfile=None, outweight=None, overwrite=False):
    infile = Path(infile)
    reffile = Path(reffile)
    outname = Path(outname)

    if outname.exists() and not overwrite:
        print(f"Output exists, skipping: {outname}")
        return

    outname.parent.mkdir(parents=True, exist_ok=True)

    # ra, dec, nx, ny, pixscale = get_ref_params(reffile)

    # cmd = [
    #     "swarp", str(infile),
    #     "-RESAMPLE", "Y",
    #     "-COMBINE", "N",
    #     "-CENTER_TYPE", "MANUAL",
    #     "-CENTER", f"{ra},{dec}",
    #     "-IMAGE_SIZE", f"{nx},{ny}",
    #     "-PIXELSCALE_TYPE", "MANUAL",
    #     "-PIXEL_SCALE", f"{pixscale:.8f}",
    #     "-IMAGEOUT_NAME", str(outname),
    #     "-WEIGHTOUT_NAME", str(outname).replace(".fits", ".swarp.weight.fits"),
    #     "-SUBTRACT_BACK", "N",
    #     "-FSCALE_DEFAULT", "1.0",
    #     "-VERBOSE_TYPE", "NORMAL",
    # ]
    ra, dec, nx, ny, pixscale = get_ref_params(reffile)

    cmd = [
        "swarp", str(infile),
        "-RESAMPLE", "Y",
        "-COMBINE", "Y",
        "-COMBINE_TYPE", "WEIGHTED",
        "-CENTER_TYPE", "MANUAL",
        "-CENTER", f"{ra},{dec}",
        "-IMAGE_SIZE", f"{nx},{ny}",
        "-PIXELSCALE_TYPE", "MANUAL",
        "-PIXEL_SCALE", f"{pixscale:.8f}",
        "-IMAGEOUT_NAME", str(outname),
        "-WEIGHTOUT_NAME", str(outweight),
        "-SUBTRACT_BACK", "N",
        "-FSCALE_DEFAULT", "1.0",
        ]

    if weightfile is not None:
        cmd += [
            "-WEIGHT_TYPE", "MAP_WEIGHT",
            "-WEIGHT_IMAGE", str(weightfile),
        ]
    else:
        cmd += ["-WEIGHT_TYPE", "NONE"]

    print("swarp command:\n",cmd)
    run(cmd)

    if outweight is not None:
        tmp_weight = Path(str(outname).replace(".fits", ".swarp.weight.fits"))
        outweight = Path(outweight)
        if tmp_weight.exists():
            tmp_weight.replace(outweight)


def main():
    parser = argparse.ArgumentParser(
        description="Use SWarp to resample new INT r-band image onto old INT Halpha grid."
    )
    parser.add_argument("infile", help="New r-band science image")
    parser.add_argument("weightfile", help="New r-band weight image")
    parser.add_argument("reffile", help="Old Halpha reference image")
    parser.add_argument("outname", help="Output reprojected r-band image")
    parser.add_argument("outweight", help="Output reprojected r-band weight image")
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    swarp_resample(
        infile=args.infile,
        weightfile=args.weightfile,
        reffile=args.reffile,
        outname=args.outname,
        outweight=args.outweight,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
