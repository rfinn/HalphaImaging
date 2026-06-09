#!/usr/bin/env python

"""
match_int_oldha_newr.py

Match:
    old INT Halpha coadds
to:
    new INT r-band coadds

using:
    telescope + DATEOBS + pointing

Then write:
    1. a CSV summary table
    2. a parallel job file for reprojection

Output reprojection target directory:
    /data-pool/coadds-v20260518/

Example output reprojection command:
    python reproject_int_r_to_old_ha.py NEW_R OLD_HA OUT_R
"""

from pathlib import Path
import pandas as pd
import re
import shutil


OLD_DIR = Path("/data-pool/Halpha/coadds-pre2025-hapy/")
NEW_DIR = Path("/data-pool/Halpha/coadds-v20260330/")
#OUT_DIR = Path("/data-pool/Halpha/coadds-v20260518/")

# redoing after I realized that both halpha and r-band coadds need to be swarped
OUT_DIR = Path("/data-pool/Halpha/coadds-v20260609/")

OUT_CSV = "int_hybrid_matches.csv"
OUT_JOBS = "reproject_jobs.txt"

import re
import numpy as np
from datetime import datetime

def vf_coord_string(ra, dec):
    return f"VF-{float(ra):07.3f}{float(dec):+07.3f}"

def parse_filename(path):
    """
    Parse filenames like:

    VF-266.040+57.998-INT-20190531-p012-r.fits
    VF-266.040+57.998-INT-20190531-p012-Halpha.fits
    VF-266.040+57.998-INT-20190531-p012-Ha6657.fits

    Returns dict with:
        ra, dec, tel, date, date_dt, pointing, filter, path
    """

    name = path.name

    m = re.match(
        r"VF-(?P<ra>[0-9.]+)(?P<dec>[+-][0-9.]+)-"
        r"(?P<tel>[A-Z]+)-(?P<date>\d{8})-"
        r"(?P<pointing>[^-]+)-(?P<filter>Halpha|Ha6657|r)\.fits$",
        name,
    )

    if m is None:
        return None

    d = m.groupdict()
    d["ra"] = float(d["ra"])
    d["dec"] = float(d["dec"])
    d["date_dt"] = datetime.strptime(d["date"], "%Y%m%d")
    d["path"] = str(path)

    return d


def angular_sep_deg(ra1, dec1, ra2, dec2):
    """
    Small-angle angular separation in degrees.
    Good enough for matching coadd centers.
    """
    dra = (ra2 - ra1) * np.cos(np.deg2rad(0.5 * (dec1 + dec2)))
    ddec = dec2 - dec1
    return np.sqrt(dra**2 + ddec**2)


def find_best_r_match(
    ha_info,
    r_infos,
    max_date_diff_days=7,
    max_sep_deg=0.25,
):
    """
    Find best r-band match for one Halpha image.

    Hard requirements:
      - same telescope
      - same pointing
      - angular separation < max_sep_deg
      - date difference <= max_date_diff_days

    Ranking:
      - smallest angular separation
      - then smallest date difference
    """

    candidates = []

    for r in r_infos:
        if r["tel"] != ha_info["tel"]:
            continue
        if r["pointing"] != ha_info["pointing"]:
            continue

        dt_days = abs((r["date_dt"] - ha_info["date_dt"]).days)
        if dt_days > max_date_diff_days:
            continue

        sep = angular_sep_deg(
            ha_info["ra"], ha_info["dec"],
            r["ra"], r["dec"],
        )

        if sep > max_sep_deg:
            continue

        candidates.append((sep, dt_days, r))

    if len(candidates) == 0:
        return None

    candidates.sort(key=lambda x: (x[0], x[1]))
    return candidates[0][2]




# ------------------------------------------------------------
# gather files
# ------------------------------------------------------------
def is_weight_file(path):
    return path.name.endswith(".weight.fits")

EXCLUDE_HA_BASENAMES = {
    "VF-234.434+59.357-INT-20190530-p040-Halpha.fits",
    "VF-255.453+23.075-INT-20190601-p008-Ha6657.fits",
    "VF-265.866+57.999-INT-20190601-p012-Ha6657.fits",
}

old_ha = sorted([
    f for f in OLD_DIR.glob("*INT*.fits")
    if (not is_weight_file(f))
    and (f.name.endswith("-Halpha.fits") or f.name.endswith("-Ha6657.fits"))
    and (f.name not in EXCLUDE_HA_BASENAMES)
    ])

print(f"number of files in old_ha = {len(old_ha)}")

print(f"Excluding {len(EXCLUDE_HA_BASENAMES)} Halpha images by explicit exclude list")
new_r = sorted([
    f for f in NEW_DIR.glob("*INT*-r.fits")
    if not is_weight_file(f)
])

new_r_weight = sorted([
    f for f in NEW_DIR.glob("*INT*-r.weight.fits")
])


# ------------------------------------------------------------
# build lookup for new r images
# ------------------------------------------------------------

r_lookup = {}

for f in new_r:
    p = parse_filename(f)

    if p is None:
        continue

    key = (p["tel"], p["date"], p["pointing"])

    r_lookup[key] = f


# ------------------------------------------------------------
# match old Ha -> new r
# ------------------------------------------------------------

ha_infos = []
failed_ha_parse = []

for f in old_ha:
    p = parse_filename(f)
    if p is None:
        failed_ha_parse.append(f)
    else:
        ha_infos.append(p)

print(f"number of old_ha files = {len(old_ha)}")
print(f"number parsed into ha_infos = {len(ha_infos)}")
print(f"number failed parse = {len(failed_ha_parse)}")

for f in failed_ha_parse[:20]:
    print("FAILED PARSE:", f.name)
    
#ha_infos = [parse_filename(f) for f in old_ha]
#ha_infos = [x for x in ha_infos if x is not None]

r_infos = [parse_filename(f) for f in new_r]
r_infos = [x for x in r_infos if x is not None]

rows = []

#print(f"length of files in ha_infos = {len(ha_infos)}")
for ha in ha_infos:
    # skip bad images
    
    r = find_best_r_match(
        ha,
        r_infos,
        max_date_diff_days=7,
        max_sep_deg=0.25,   # stricter than 0.5 deg
    )

    if r is None:
        print(f"NO MATCH: {Path(ha['path']).name}")
        continue

    sep = angular_sep_deg(ha["ra"], ha["dec"], r["ra"], r["dec"])
    dt_days = abs((r["date_dt"] - ha["date_dt"]).days)


    rfile = Path(r["path"])
    ha_file = Path(ha["path"])

    rweight = rfile.with_name(rfile.name.replace("-r.fits", "-r.weight.fits"))
    ha_weight = ha_file.with_name(ha_file.name.replace(".fits", ".weight.fits"))

    vfroot = vf_coord_string(ha["ra"], ha["dec"])

    # Preserve Halpha vs Ha6657 suffix
    if ha_file.name.endswith("-Halpha.fits"):
        ha_suffix = "Halpha"
    elif ha_file.name.endswith("-Ha6657.fits"):
        ha_suffix = "Ha6657"
    else:
        ha_suffix = "Ha"

    out_ha_name = (
        f"{vfroot}-"
        f"{ha['tel']}-"
        f"{ha['date']}-"
        f"{ha['pointing']}-{ha_suffix}.fits"
    )

    out_r_name = (
        f"{vfroot}-"
        f"{ha['tel']}-"
        f"{ha['date']}-"
        f"{ha['pointing']}-r.fits"
    )

    out_ha_weight_name = out_ha_name.replace(".fits", ".weight.fits")
    out_r_weight_name = out_r_name.replace("-r.fits", "-r.weight.fits")

    out_ha = OUT_DIR / out_ha_name
    out_ha_weight = OUT_DIR / out_ha_weight_name

    out_r = OUT_DIR / out_r_name
    out_r_weight = OUT_DIR / out_r_weight_name

    rows.append({
        "tel": ha["tel"],
        "ha_dateobs": ha["date"],
        "r_dateobs": r["date"],
        "date_diff_days": dt_days,
        "pointing": ha["pointing"],
        "ha_ra": ha["ra"],
        "ha_dec": ha["dec"],
        "r_ra": r["ra"],
        "r_dec": r["dec"],
        "sep_deg": sep,
        "ha_filter": ha["filter"],
        "r_filter": r["filter"],
        "old_ha": str(ha_file),
        "old_ha_weight": str(ha_weight),
        "new_r": str(rfile),
        "new_r_weight": str(rweight),
        "out_ha": str(out_ha),
        "out_ha_weight": str(out_ha_weight),
        "out_r": str(out_r),
        "out_r_weight": str(out_r_weight),
    })


for row in rows:
    if row["sep_deg"] > 0.05 or row["date_diff_days"] > 2:
        print(
            "CHECK:",
            row["pointing"],
            row["ha_dateobs"],
            row["r_dateobs"],
            f"sep={row['sep_deg']:.3f} deg",
            Path(row["old_ha"]).name,
            Path(row["new_r"]).name,
        )



df = pd.DataFrame(rows)

print()
print(f"Matched {len(df)} INT Ha/r pairs")

df.to_csv(OUT_CSV, index=False)
print(f"Wrote {OUT_CSV}")


# ------------------------------------------------------------
# write GNU parallel jobs file
# ------------------------------------------------------------
with open(OUT_JOBS, "w") as outfile:
    for _, row in df.iterrows():

        # First: resample old Halpha through SWarp onto the chosen reference grid.
        outfile.write(
            f"{row['old_ha']} "
            f"{row['old_ha_weight']} "
            f"{row['old_ha']} "
            f"{row['out_ha']} "
            f"{row['out_ha_weight']}\n"
        )

        # Second: resample matched r-band onto the same old-Halpha reference grid.
        outfile.write(
            f"{row['new_r']} "
            f"{row['new_r_weight']} "
            f"{row['old_ha']} "
            f"{row['out_r']} "
            f"{row['out_r_weight']}\n"
        )


print(f"Wrote {OUT_JOBS}")
