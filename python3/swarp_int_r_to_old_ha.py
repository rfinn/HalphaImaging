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




def repair_swarp_header(out_image, source_image, ref_image, out_weight=None):
    """
    Add hybrid-product provenance after SWarp.

    SWarp already copies science/calibration metadata through COPY_KEYWORDS.
    This function only adds bookkeeping keywords.
    """
    out_image = Path(out_image)
    source_image = Path(source_image)
    ref_image = Path(ref_image)

    with fits.open(out_image, mode="update") as hdul:
        hdr = hdul[0].header

        hdr["SRCIMG"] = (source_image.name, "Source image before SWarp")
        hdr["REFIMG"] = (ref_image.name, "Reference image/grid for SWarp")
        hdr["HYBRID"] = (True, "Hybrid INT product")
        hdr["REPROJ"] = (True, "Image resampled with SWarp")

        if "-r.fits" in out_image.name:
            hdr["RIMAGE"] = (out_image.name, "Hybrid reprojected r image")
            hdr["HAIMAGE"] = (ref_image.name, "Reference Halpha image")
            hdr["SRC_R"] = (source_image.name, "Source r image before reprojection")
            hdr["REF_HA"] = (ref_image.name, "Reference Halpha image")
            hdr["RREDUCT"] = ("v20260330", "r-band reduction version")
            hdr["HARED"] = ("pre2025", "Halpha reduction version")

        elif "-Halpha.fits" in out_image.name or "-Ha6657.fits" in out_image.name:
            hdr["HAIMAGE"] = (out_image.name, "Hybrid reprojected Halpha image")
            hdr["SRC_HA"] = (source_image.name, "Source Halpha image before reprojection")
            hdr["REF_HA"] = (ref_image.name, "Reference Halpha image")
            hdr["HARED"] = ("pre2025", "Halpha reduction version")

        hdul.flush()

    if out_weight is not None:
        out_weight = Path(out_weight)

        if out_weight.exists():
            with fits.open(out_weight, mode="update") as hdul:
                hdr = hdul[0].header

                hdr["BUNIT"] = ("weight", "Weight map")
                hdr["SRCIMG"] = (source_image.name, "Source image before SWarp")
                hdr["REFIMG"] = (ref_image.name, "Reference image/grid for SWarp")
                hdr["HYBRID"] = (True, "Hybrid INT weight product")
                hdr["REPROJ"] = (True, "Weight map resampled with SWarp")

                if "-r.weight.fits" in out_weight.name:
                    hdr["RIMAGE"] = (out_image.name, "Associated reprojected r image")
                    hdr["SRC_R"] = (source_image.name, "Source r image before reprojection")
                    hdr["REF_HA"] = (ref_image.name, "Reference Halpha image")

                elif "-Halpha.weight.fits" in out_weight.name or "-Ha6657.weight.fits" in out_weight.name:
                    hdr["HAIMAGE"] = (out_image.name, "Associated reprojected Halpha image")
                    hdr["SRC_HA"] = (source_image.name, "Source Halpha image before reprojection")
                    hdr["REF_HA"] = (ref_image.name, "Reference Halpha image")

                hdul.flush()
                
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
def get_ref_params_v2(reffile):
    data, hdr = fits.getdata(reffile, header=True)
    w = WCS(hdr)

    ny, nx = data.shape

    # Use the reference image WCS reference coordinate, not the
    # sky coordinate of the array midpoint.
    ra = float(hdr.get("CRVAL1", w.wcs.crval[0]))
    dec = float(hdr.get("CRVAL2", w.wcs.crval[1]))

    # pixel scale in arcsec/pixel
    pscale = proj_plane_pixel_scales(w)  # deg/pix
    pixscale = float(abs(pscale[0]) * 3600.0)

    return ra, dec, int(nx), int(ny), pixscale
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
        "-COPY_KEYWORDS", "OBJECT,FILTER,INSTRMNT,TELESCOP,GAIN,EPOCH,DATE-OBS,MJD-OBS,AIRMASS,EXPTIME,PHOTZP,MAGZP,PAN_FZP,FLTRATIO,FILTER_RATIO,SKYMED,SKYSTD,BUNIT,RIMAGE,HAIMAGE,SEFWHM,FWHM",
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
