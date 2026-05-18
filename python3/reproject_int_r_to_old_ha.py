#!/usr/bin/env python

import argparse
from pathlib import Path
from hapy.imagetools.imutils import reproject_image


def main():
    parser = argparse.ArgumentParser(
        description="Reproject one FITS image onto the WCS/grid of a reference FITS image."
    )
    parser.add_argument("infile", help="Input FITS image to reproject")
    parser.add_argument("reffile", help="Reference FITS image defining output WCS/grid")
    parser.add_argument("outname", help="Output reprojected FITS image")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output")

    args = parser.parse_args()

    infile = Path(args.infile)
    reffile = Path(args.reffile)
    outname = Path(args.outname)

    if not infile.exists():
        raise FileNotFoundError(infile)
    if not reffile.exists():
        raise FileNotFoundError(reffile)

    if outname.exists() and not args.overwrite:
        print(f"Output exists, skipping: {outname}")
        return

    if outname.exists() and args.overwrite:
        outname.unlink()

    outname.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reprojecting: {infile}")
    print(f"Reference:    {reffile}")
    print(f"Output:       {outname}")

    reproject_image(str(infile), str(reffile), str(outname))


if __name__ == "__main__":
    main()
