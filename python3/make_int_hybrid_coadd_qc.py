#!/usr/bin/env python

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.visualization import simple_norm

COADD_DIR = Path("/data-pool/Halpha/coadds-v20260609")
OUTPDF = "int_halpha_r_coadd_qc.pdf"


def get_r_for_ha(hfile):
    rfile = hfile.name.replace("-Halpha.fits", "-r.fits")
    rfile = rfile.replace("-Ha6657.fits", "-r.fits")
    return hfile.with_name(rfile)


def read_downsample(path, step=8):
    data = fits.getdata(path)
    return data[::step, ::step]


def show_image(ax, data, title,percent=99.5):
    norm = simple_norm(data, stretch="asinh", percent=percent)
    ax.imshow(data, origin="lower", cmap="gray", norm=norm)
    ax.set_title(title, fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])


ha_files = sorted(
    list(COADD_DIR.glob("*INT*-Halpha.fits")) +
    list(COADD_DIR.glob("*INT*-Ha6657.fits"))
)

print(f"Found {len(ha_files)} halpha files")

pairs = []
missing_r = []

for h in ha_files:
    rfile = get_r_for_ha(h)

    if rfile.exists():
        pairs.append((h, rfile))
    else:
        missing_r.append(h)

print(f"Found {len(pairs)} INT Halpha/r pairs")

if missing_r:
    print(f"\nMissing r-band match for {len(missing_r)} Halpha files:")
    for h in missing_r:
        print(f"  {h.name}")
    
# ha_files = sorted(
#     list(COADD_DIR.glob("*INT*-Halpha.fits")) +
#     list(COADD_DIR.glob("*INT*-Ha6657.fits"))
# )

# print(f"Found {len(ha_files)} halpha files")

# pairs = [(h, get_r_for_ha(h)) for h in ha_files if get_r_for_ha(h).exists()]

# print(f"Found {len(pairs)} INT Halpha/r pairs")

pairs_per_page = 20
nrows = 10
ncols = 4   # Ha, r, Ha, r
page = 1 
with PdfPages(OUTPDF) as pdf:
    for start in range(0, len(pairs), pairs_per_page):
        subset = pairs[start:start + pairs_per_page]

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(11, 18),
            constrained_layout=True,
        )

        axes = axes.ravel()

        for ax in axes:
            ax.axis("off")

        for j, (hfile, rfile) in enumerate(subset):
            row = j % nrows
            block = j // nrows
            ax_ha = axes[row * ncols + block * 2]
            ax_r = axes[row * ncols + block * 2 + 1]

            try:
                hdata = read_downsample(hfile)
                rdata = read_downsample(rfile)

                root = hfile.name.replace("-Halpha.fits", "").replace("-Ha6657.fits", "")

                show_image(ax_ha, hdata, f"Ha\n{root}",percent=99.9)
                show_image(ax_r, rdata, "r")

            except Exception as e:
                ax_ha.text(0.5, 0.5, f"ERROR\n{hfile.name}\n{e}",
                           ha="center", va="center", fontsize=6)
                ax_r.axis("off")

        fig.suptitle(
            f"INT hybrid coadd QC: pairs {start+1}-{start+len(subset)}",
            fontsize=14,
        )
        if (start % nrows*ncols//2) == 0:
            print(f"done with page {page}")
            page += 1
        pdf.savefig(fig)
        plt.close(fig)

print(f"Wrote {OUTPDF}")
