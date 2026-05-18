#!/usr/bin/env python

import argparse
from pathlib import Path
from hapy.imagetools.imutils import reproject_image

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reproject an r-band FITS image and optional weight image "
            "onto the WCS/grid of a reference Halpha FITS image."
        )
    )
    parser.add_argument("infile", help="Input science FITS image to reproject")
    parser.add_argument("weightfile", help="Input weight FITS image to reproject")
    parser.add_argument("reffile", help="Reference FITS image defining output WCS/grid")
    parser.add_argument("outname", help="Output reprojected science FITS image")
    parser.add_argument("outweight", help="Output reprojected weight FITS image")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")

    args = parser.parse_args()

    infile = Path(args.infile)
    weightfile = Path(args.weightfile)
    reffile = Path(args.reffile)
    outname = Path(args.outname)
    outweight = Path(args.outweight)

    for p in [infile, weightfile, reffile]:
        if not p.exists():
            raise FileNotFoundError(p)

    outname.parent.mkdir(parents=True, exist_ok=True)
    outweight.parent.mkdir(parents=True, exist_ok=True)

    if outname.exists() and not args.overwrite:
        print(f"Science output exists, skipping: {outname}")
    else:
        if outname.exists() and args.overwrite:
            outname.unlink()

        print(f"Reprojecting science: {infile}")
        print(f"Reference:            {reffile}")
        print(f"Output:               {outname}")

        reproject_image(
            str(infile),
            str(reffile),
            str(outname),
            overwrite=args.overwrite,
        )

    if outweight.exists() and not args.overwrite:
        print(f"Weight output exists, skipping: {outweight}")
    else:
        if outweight.exists() and args.overwrite:
            outweight.unlink()

        print(f"Reprojecting weight:  {weightfile}")
        print(f"Reference:            {reffile}")
        print(f"Output weight:        {outweight}")

        reproject_image(
            str(weightfile),
            str(reffile),
            str(outweight),
            overwrite=args.overwrite,
        )




if __name__ == "__main__":
    main()
