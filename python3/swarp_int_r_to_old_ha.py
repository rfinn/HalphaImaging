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

from astropy.io import fits
from pathlib import Path
import numpy as np


SCIENCE_KEYS_FROM_SOURCE = [
    "FILTER", "FILTER1", "FILTER2",
    "EXPTIME",
    "PHOTZP", "MAGZP", "PAN_FZP",
    "SKYMED", "SKYSTD",
    "BUNIT", "GAIN", "RDNOISE",
    "FSCALE",
    "FLTRATIO",
    "FILTER_RATIO",
]

COORD_KEYS_FROM_REF = [
    "RA", "DEC", "OBJRA", "OBJDEC",
]


def repair_swarp_header(out_r, source_r, ref_ha, out_weight=None):
    """
    Repair SWarp output headers for hybrid INT r-band products.

    out_r    : reprojected r image
    source_r : original new r image before reprojection
    ref_ha   : old Halpha image defining final WCS/grid
    """
    out_r = Path(out_r)
    source_r = Path(source_r)
    ref_ha = Path(ref_ha)

    with fits.open(out_r, mode="update") as hout, \
         fits.open(source_r) as hsrc, \
         fits.open(ref_ha) as href:

        hdr = hout[0].header
        src = hsrc[0].header
        ref = href[0].header

        # science/calibration metadata should follow source r image
        for key in SCIENCE_KEYS_FROM_SOURCE:
            if key in src:
                hdr[key] = src[key]

        # center/object coordinate keywords should follow reference Halpha grid
        for key in COORD_KEYS_FROM_REF:
            if key in ref:
                hdr[key] = ref[key]

        # matched image bookkeeping
        hdr["RIMAGE"] = (out_r.name, "Hybrid reprojected r image")
        hdr["HAIMAGE"] = (ref_ha.name, "Matched old Halpha image")
        hdr["SRC_R"] = (source_r.name, "Source r image before reprojection")
        hdr["REF_HA"] = (ref_ha.name, "Reference Halpha image")
        hdr["HYBRID"] = (True, "Hybrid INT old-Halpha/new-r product")
        hdr["RREDUCT"] = ("v20260330", "r-band reduction version")
        hdr["HARED"] = ("pre2025", "Halpha reduction version")
        hdr["REPROJ"] = (True, "r image reprojected to Halpha grid")
        hdr["REFTYPE"] = ("Halpha", "Reference grid for reprojection")

        hout.flush()

    if out_weight is not None:
        out_weight = Path(out_weight)
        if out_weight.exists():
            with fits.open(out_weight, mode="update") as hw, \
                 fits.open(source_r) as hsrc, \
                 fits.open(ref_ha) as href:

                hdr = hw[0].header
                src = hsrc[0].header
                ref = href[0].header

                for key in COORD_KEYS_FROM_REF:
                    if key in ref:
                        hdr[key] = ref[key]

                hdr["FILTER"] = src.get("FILTER", "r")
                hdr["BUNIT"] = ("weight", "Weight map")
                hdr["RIMAGE"] = (out_r.name, "Associated reprojected r image")
                hdr["HAIMAGE"] = (ref_ha.name, "Matched old Halpha image")
                hdr["SRC_R"] = (source_r.name, "Source r image before reprojection")
                hdr["REF_HA"] = (ref_ha.name, "Reference Halpha image")
                hdr["HYBRID"] = (True, "Hybrid INT old-Halpha/new-r weight")
                hdr["REPROJ"] = (True, "Weight map reprojected to Halpha grid")
                hw.flush()

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
        "-COPY_KEYWORDS", "OBJECT,FILTER,INSTRMNT,TELESCOP,GAIN,EPOCH,DATE-OBS,MJD-OBS,AIRMASS,EXPTIME,PHOTZP,MAGZP,PAN_FZP,FLTRATIO,FILTER_RATIO,SKYMED,SKYSTD,BUNIT,RIMAGE,HAIMAGE",
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

    # repair headers after SWarp output exists
    repair_swarp_header(
        out_r=outname,
        source_r=infile,
        ref_ha=reffile,
        out_weight=outweight,
        )

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
