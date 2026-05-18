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



### 1. Added generic reprojection utility to `hapy`

```python
reproject_image(infile, reffile, outname)
```

using:

```python
reproject_adaptive(..., conserve_flux=True)
```

---

### 2. Build matching script concept

Need script to:

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

---

```
python ~/github/HalphaImaging/python3/match_int_oldha_newr.py 
```

## Reproject r-band weight images too

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
parallel --bar -j 4 --memfree 50G python ~/github/HalphaImaging/python3/reproject_int_r_to_old_ha.py {1} {2} {3} {4} {5} :::: reproject_jobs.txt
```
  
### 3. Write reprojection jobs file

Output:

```text
reproject_jobs.txt
```

containing:

```text
NEW_R   OLD_HA   OUTPUT_R
```

for GNU parallel.

---

### 4. Reproject outputs go to

```text
/data-pool/coadds-v20260518/
```

---

### 5. Copy old INT Hα coadds into

```text
/data-pool/coadds-v20260518/
```

so the hybrid dataset becomes self-contained.

---

### 6. Need header bookkeeping

Update:

```text
RIMAGE
HAIMAGE
```

plus provenance keywords:

```text
HYBRID
RREDUCT
HARED
REPROJ
```

---

## Immediate next step

I think the next concrete action is:

### Run the matching script

to generate:

```text
int_hybrid_matches.csv
reproject_jobs.txt
```

Then:

```bash
parallel --colsep ' ' -j 6 \
python reproject_int_r_to_old_ha.py {1} {2} {3} \
:::: reproject_jobs.txt
```

After that:

* copy old Hα files into the new directory
* update headers
* rebuild cutouts
* rerun `run_analysis` for INT only.
