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
