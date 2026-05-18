# Overview
The duplicate validation and old/new reduction comparisons indicate that the new INT r-band reductions are robust and significantly improved relative to the older reductions, with very small spatially uniform residuals (≲0.02 dex) and improved duplicate reproducibility. In contrast, the new INT Hα reductions exhibit substantially larger duplicate scatter and clear spatially coherent residual structure across the detector. Rebuilding the CS-ZP images locally from the cutouts reduced some of the spatial systematics, but did not fully recover the duplicate consistency seen in the older Hα reductions. Additional tests using polynomial flattening with two rounds of illumination correction also failed to reproduce the quality of the older reductions, suggesting that the issue is more fundamental to the newer INT Hα processing. Based on these results, the current plan is to adopt a hybrid INT reduction for the catalog paper: use the older INT Hα coadds together with the newer INT r-band coadds, reprojecting the r-band images onto the native Hα grid to preserve the Hα image quality and PSF.

## Current decision on workflow:

```text
old INT Hα  +  new INT r
```

with:

```text
reproject new r → old Hα grid
```

and keep standard downstream naming:

```text
*-R.fits
*-Ha.fits
```

## Process



### 1. Build matching script 

* match:

  ```text
  old INT Hα
  ↔
  new INT r
  ```

using:

* telescope
* DATEOBS
* pointing


Script:
```
python ~/github/HalphaImaging/python3/match_int_oldha_newr.py 
```

## 2. Reproject new INT r-band images and weights images onto prior INT Hα coadds
For each hybrid INT field:

- source science image:
  - new INT `*-r.fits`
- source weight image:
  - new INT `*-r.weight.fits`
- reference image:
  - old INT `*-Halpha.fits`

Create:

- reprojected science image:
  - hybrid `*-r.fits`
- reprojected weight image:
  - hybrid `*-r.weight.fits`
  
```
parallel --bar -j 8 --joblog swarp_int_hybrid.joblog --results swarp_int_hybrid_logs --colsep ' ' python ~/github/HalphaImaging/python3/swarp_int_r_to_old_ha.py {1} {2} {3} {4} {5} :::: reproject_jobs.txt
```
  
## 3. Copy old INT Hα coadds, new BOK, HDI, and MOS coadds into new coadd directory
so the hybrid dataset becomes self-contained.


The new coadd directory:
```text
/data-pool/coadds-v20260518/
```

Script to copy the files:
```
python  ~/github/HalphaImaging/python3/build_hybrid_int_coadds.py 
```
NOTE: this also copies the gaia catalogs and panstarrs csv files.

## 4. Make INT r had `HAIMAGE` in header and INT Halpha has `RIMAGE` in header

Update:

```text
RIMAGE
HAIMAGE
```

```
python ~/github/HalphaImaging/python3/uat_add_haimage_to_rheader.py --filestring VF --filestring2 INT --vfs
```

ls 
## 5. Rebuild INT r-band PSF images

### Approach

create `/data-pool/psf-images-v20260518/`

Populate it with:
- All non-INT / non-hybrid PSFs copied from:
`/data-pool/psf-images/`
- Old INT Hα PSFs copied from:
`/data-pool/psf-images-pre2025/`
- New INT r-band PSFs generated for the reprojected INT r coadds.

Also copy matching diagnostic plots into:

`/data-pool/psf-images-v20260518/plots/`

This keeps provenance clean and avoids contaminating the existing PSF directory.

Rationale:
- /data-pool/psf-images/ remains the current/new-coadd PSF set.
- /data-pool/psf-images-pre2025/ remains the old PSF set.
- /data-pool/psf-images-v20260518/ becomes the hybrid production PSF set.
- run_analysis --psf-dir /data-pool/psf-images-v20260518/ is unambiguous.



### Building new psf images 


```
/data-pool/Halpha/psf-images-v20260518
```

Create an input file list:
```
find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-INT-*-r.fits" | sort > int_hybrid_r_images.txt
```
```
cp /data-pool/github/hapy/hapy/astromatic/default.param .
cp /data-pool/github/hapy/hapy/astromatic/default.conv .
cp /data-pool/github/hapy/hapy/astromatic/default.nnw .
```

Test one:
```
python -m hapy.imagetools.buildpsf --image /data-pool/Halpha/coadds-v20260518/VF-256.898+60.794-INT-20190603-p010-r.fits --int --overwrite
```


Then build psfs in parallel:
```
parallel --bar -j 16 --joblog buildpsf_int_hybrid_r.joblog --results buildpsf_int_hybrid_r_logs python -m hapy.imagetools.buildpsf --image {} --int --overwrite :::: int_hybrid_r_images.txt
```

Rerunning on old INT halpha, just to be on the safe side
```
find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-INT-*-Halpha.fits" | sort > int_hybrid_Halpha_images.txt
```

```
find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-INT-*-Ha6657.fits" | sort >> int_hybrid_Halpha_images.txt
```

```
parallel --bar -j 16 --joblog buildpsf_int_hybrid_Halpha.joblog --results buildpsf_int_hybrid_Halpha_logs python -m hapy.imagetools.buildpsf --image {} --int --overwrite :::: int_hybrid_Halpha_images.txt
```


and for all coadds while I'm at it...

```
 find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-BOK-*-r.fits" | sort > non_int_hybrid_images.txt
 find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-BOK-*-Ha4.fits" | sort >> non_int_hybrid_images.txt
 find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-HDI-*-r.fits" | sort >> non_int_hybrid_images.txt
 find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-HDI-*-R.fits" | sort >> non_int_hybrid_images.txt
 find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-MOS-*-R.fits" | sort >> non_int_hybrid_images.txt
 find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-MOS-*-ha4.fits" | sort >> non_int_hybrid_images.txt
 find /data-pool/Halpha/coadds-v20260518 -maxdepth 1 -type f -name "*-HDI-*-ha4.fits" | sort >> non_int_hybrid_images.txt
```


```
parallel --bar -j 16 --joblog buildpsf_non_int_hybrid.joblog --results buildpsf_non_int_hybrid_logs python -m hapy.imagetools.buildpsf --image {} --int --overwrite :::: non_int_hybrid_images.txt
```

