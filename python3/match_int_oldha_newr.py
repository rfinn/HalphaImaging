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


OLD_DIR = Path("/data-pool/Halpha/coadds-pre2025-hapy/")
NEW_DIR = Path("/data-pool/Halpha/coadds-v20260330/")
OUT_DIR = Path("/data-pool/Halpha/coadds-v20260518/")

OUT_CSV = "int_hybrid_matches.csv"
OUT_JOBS = "reproject_jobs.txt"


def parse_filename(path):
    """
    Parse filenames like:

    VF-266.040+57.998-INT-20190531-p012-r.fits
    """

    name = path.name

    m = re.match(
        r"VF-[^/]+-(?P<tel>[A-Z]+)-(?P<date>\d{8})-(?P<pointing>p\d+)-(?P<filter>.+)\.fits",
        name,
    )

    if m is None:
        return None

    d = m.groupdict()
    d["path"] = str(path)

    return d


# ------------------------------------------------------------
# gather files
# ------------------------------------------------------------

old_ha = sorted(OLD_DIR.glob("*INT*-Ha*.fits"))
new_r = sorted(NEW_DIR.glob("*INT*-r.fits"))

print(f"Found {len(old_ha)} old INT Halpha images")
print(f"Found {len(new_r)} new INT r-band images")


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

rows = []

for ha in old_ha:

    p = parse_filename(ha)

    if p is None:
        continue

    key = (p["tel"], p["date"], p["pointing"])

    if key not in r_lookup:
        print(f"NO MATCH: {ha.name}")
        continue

    rfile = r_lookup[key]

    # output reprojected r filename
    out_r = OUT_DIR / rfile.name

    rows.append({
        "tel": p["tel"],
        "dateobs": p["date"],
        "pointing": p["pointing"],
        "old_ha": str(ha),
        "new_r": str(rfile),
        "out_r": str(out_r),
    })


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

        outfile.write(
            f"{row['new_r']} "
            f"{row['old_ha']} "
            f"{row['out_r']}\n"
        )

print(f"Wrote {OUT_JOBS}")


# ------------------------------------------------------------
# example parallel command
# ------------------------------------------------------------

print()
print("Example:")
print(
    "parallel --colsep ' ' -j 6 "
    "python reproject_int_r_to_old_ha.py {1} {2} {3} "
    f":::: {OUT_JOBS}"
)
