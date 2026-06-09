#!/usr/bin/env python

from pathlib import Path
import shutil

OLD_DIR = Path("/data-pool/Halpha/coadds-pre2025-hapy/")
NEW_DIR = Path("/data-pool/Halpha/coadds-v20260330/")
#OUT_DIR = Path("/data-pool/Halpha/coadds-v20260518/")
# after realizing that both the halpha and r-band INT coadds need to be reprojected
OUT_DIR = Path("/data-pool/Halpha/coadds-v20260518/")

GAIA_SUBDIR = "gaia_catalogs"   # change to "gaia_catalog" if needed


def copy_file(src, dst, overwrite=False):
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() and not overwrite:
        print(f"exists, skipping: {dst}")
        return

    print(f"copying: {src} -> {dst}")
    shutil.copy2(src, dst)

def copy_associated_files(src, src_dir, out_dir, overwrite=False, copy_weight=True):
    """
    Copy weight, pan_tab CSVs, and Gaia catalogs associated with src.
    """
    stem = src.name.replace(".fits", "")

    # weight file
    if copy_weight:
        weight = src.with_name(stem + ".weight.fits")
        if weight.exists():
            copy_file(weight, out_dir / weight.name, overwrite=overwrite)

    # matching pan_tab files in coadd dir
    for p in src_dir.glob(f"{stem}*pan_tab.csv"):
        copy_file(p, out_dir / p.name, overwrite=overwrite)

    # matching Gaia catalogs
    for gaia_name in [GAIA_SUBDIR, "gaia_catalog"]:
        gaia_dir = src_dir / gaia_name
        if gaia_dir.exists():
            out_gaia_dir = out_dir / gaia_name
            for g in gaia_dir.glob(f"{stem}*"):
                copy_file(g, out_gaia_dir / g.name, overwrite=overwrite)



def copy_science_with_aux(src, src_dir, out_dir, overwrite=False):
    copy_file(src, out_dir / src.name, overwrite=overwrite)
    copy_associated_files(src, src_dir, out_dir, overwrite=overwrite)


def is_weight_file(p):
    return p.name.endswith(".weight.fits")


def main(overwrite=False):
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # 1. INT Halpha images from OLD_DIR
    # ------------------------------------------------------------
    int_ha = sorted([
        f for f in OLD_DIR.glob("*INT*.fits")
        if not is_weight_file(f)
        and (
            f.name.endswith("-Halpha.fits")
            or f.name.endswith("-Ha6657.fits")
        )
    ])

    print(f"\nCopying {len(int_ha)} old INT Halpha coadds")
    for f in int_ha:
        #copy_science_with_aux(f, OLD_DIR, OUT_DIR, overwrite=overwrite)
        copy_associated_files(f, OLD_DIR, OUT_DIR, overwrite=overwrite, copy_weight=False)

    # ------------------------------------------------------------
    # 2. All BOK coadds from NEW_DIR
    # ------------------------------------------------------------
    bok = sorted([
        f for f in NEW_DIR.glob("*BOK*.fits")
        if not is_weight_file(f)
    ])

    print(f"\nCopying {len(bok)} BOK coadds")
    for f in bok:
        copy_science_with_aux(f, NEW_DIR, OUT_DIR, overwrite=overwrite)

    # ------------------------------------------------------------
    # 3. All HDI and MOS coadds from NEW_DIR
    # ------------------------------------------------------------
    hdi_mos = sorted([
        f for f in NEW_DIR.glob("*.fits")
        if not is_weight_file(f)
        and ("HDI" in f.name or "MOS" in f.name)
    ])

    print(f"\nCopying {len(hdi_mos)} HDI/MOS coadds")
    for f in hdi_mos:
        copy_science_with_aux(f, NEW_DIR, OUT_DIR, overwrite=overwrite)


if __name__ == "__main__":
    main(overwrite=False)
