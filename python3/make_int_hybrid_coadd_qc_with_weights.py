#!/usr/bin/env python

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.visualization import simple_norm, ImageNormalize, LinearStretch

COADD_DIR = Path("/data-pool/Halpha/coadds-v20260518")
OUTPDF = "int_halpha_r_coadd_weight_qc.pdf"


def get_r_for_ha(hfile):
    rfile = hfile.name.replace("-Halpha.fits", "-r.fits")
    rfile = rfile.replace("-Ha6657.fits", "-r.fits")
    return hfile.with_name(rfile)


def get_weight_for_image(imfile):
    candidates = [
        imfile.with_name(imfile.name.replace(".fits", ".weight.fits")),
        imfile.with_name(imfile.name.replace(".fits", "-weight.fits")),
        imfile.with_name(imfile.name.replace(".fits", "_weight.fits")),
        imfile.with_name(imfile.name.replace(".fits", ".wt.fits")),
    ]

    for c in candidates:
        if c.exists():
            return c

    return None


def read_downsample(path, step=8):
    data = fits.getdata(path)
    return data[::step, ::step]


def show_image(ax, data, title):
    norm = simple_norm(data, stretch="asinh", percent=99.5)
    ax.imshow(data, origin="lower", cmap="gray", norm=norm)
    ax.set_title(title, fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])


def show_weight(ax, data, title, norm):
    ax.imshow(data, origin="lower", cmap="gray", norm=norm)
    ax.set_title(title, fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])


def finite_minmax(data):
    d = np.asarray(data, dtype=float)
    good = np.isfinite(d)
    if not np.any(good):
        return None
    return np.nanmin(d[good]), np.nanmax(d[good])
def collect_weight_values(image_files):
    values = []
    for imfile in image_files:
        wfile = get_weight_for_image(imfile)
        if wfile is None:
            continue
        try:
            wdata = read_downsample(wfile)
            good = np.isfinite(wdata)
            if np.any(good):
                values.append(wdata[good].ravel())
        except Exception:
            pass
    return values


def make_weight_norm(values, label):
    if not values:
        print(f"WARNING: no {label} weight images found")
        return None

    allw = np.concatenate(values)

    # ignore exact zeros for stretch if nonzero pixels exist
    nonzero = allw[np.isfinite(allw) & (allw != 0)]
    if len(nonzero) > 0:
        use = nonzero
    else:
        use = allw[np.isfinite(allw)]

    wmin = np.nanpercentile(use, 1)
    wmax = np.nanpercentile(use, 99)

    if not np.isfinite(wmin) or not np.isfinite(wmax) or wmin == wmax:
        wmin = np.nanmin(use)
        wmax = np.nanmax(use)

    print(f"{label} weight stretch: vmin={wmin:.4g}, vmax={wmax:.4g}")

    return ImageNormalize(vmin=wmin, vmax=wmax, stretch=LinearStretch())

ha_files = sorted(
    list(COADD_DIR.glob("*INT*-Halpha.fits")) +
    list(COADD_DIR.glob("*INT*-Ha6657.fits"))
)

pairs = [(h, get_r_for_ha(h)) for h in ha_files if get_r_for_ha(h).exists()]

print(f"Found {len(pairs)} INT Halpha/r pairs")

# --------------------------------------------------
# Find global weight stretch from all downsampled weights
# --------------------------------------------------
weight_values = []

for hfile, rfile in pairs:
    for imfile in [hfile, rfile]:
        wfile = get_weight_for_image(imfile)
        if wfile is None:
            continue

        try:
            wdata = read_downsample(wfile)
            good = np.isfinite(wdata)
            if np.any(good):
                weight_values.append(wdata[good].ravel())
        except Exception:
            pass

if weight_values:
    allw = np.concatenate(weight_values)
    wmin = np.nanpercentile(allw, 1)
    wmax = np.nanpercentile(allw, 99)

    if not np.isfinite(wmin) or not np.isfinite(wmax) or wmin == wmax:
        wmin = np.nanmin(allw)
        wmax = np.nanmax(allw)

    #weight_norm = ImageNormalize(vmin=wmin, vmax=wmax, stretch=LinearStretch())

    ha_weight_norm = make_weight_norm(
        collect_weight_values([h for h, r in pairs]),
        "Ha",
        )

    r_weight_norm = make_weight_norm(
        collect_weight_values([r for h, r in pairs]),
        "r",
        )

    print(f"Global weight stretch: vmin={wmin:.4g}, vmax={wmax:.4g}")
else:
    weight_norm = None
    print("WARNING: no weight images found")


pairs_per_page = 10
nrows = 10
ncols = 4   # Ha, Ha weight, r, r weight

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
            row = j
            ax_ha = axes[row * ncols + 0]
            ax_hw = axes[row * ncols + 1]
            ax_r = axes[row * ncols + 2]
            ax_rw = axes[row * ncols + 3]

            try:
                hdata = read_downsample(hfile)
                rdata = read_downsample(rfile)

                hweight = get_weight_for_image(hfile)
                rweight = get_weight_for_image(rfile)

                root = hfile.name.replace("-Halpha.fits", "").replace("-Ha6657.fits", "")

                show_image(ax_ha, hdata, f"Ha\n{root}")
                show_image(ax_r, rdata, "r")

                if hweight is not None and ha_weight_norm is not None:
                    hwdata = read_downsample(hweight)
                    show_weight(ax_hw, hwdata, "Ha weight", weight_norm)
                else:
                    ax_hw.text(0.5, 0.5, "No Ha weight", ha="center", va="center", fontsize=6)

                if rweight is not None and r_weight_norm is not None:
                    rwdata = read_downsample(rweight)
                    show_weight(ax_rw, rwdata, "r weight", weight_norm)
                else:
                    ax_rw.text(0.5, 0.5, "No r weight", ha="center", va="center", fontsize=6)

            except Exception as e:
                ax_ha.text(
                    0.5, 0.5,
                    f"ERROR\n{hfile.name}\n{e}",
                    ha="center", va="center", fontsize=6,
                )

        fig.suptitle(
            f"INT hybrid coadd + weight QC: pairs {start+1}-{start+len(subset)}",
            fontsize=14,
        )
        pdf.savefig(fig)
        plt.close(fig)

print(f"Wrote {OUTPDF}")
