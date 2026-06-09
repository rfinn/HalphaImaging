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



In new coadd directory (e.g. `/data-pool/Halpha/coadds-v20260609`):
Script:
```
python ~/github/HalphaImaging/python3/match_int_oldha_newr.py 
```

Output:
```
(hapy) rfinn@draco:/data-pool/Halpha/coadds-v20260609$ python ~/github/HalphaImaging/python3/match_int_oldha_newr.py
number of files in old_ha = 119
Excluding 3 Halpha images by explicit exclude list
number of old_ha files = 119
number parsed into ha_infos = 119
number failed parse = 0
CHECK: VFID1168 20220501 20220501 sep=0.054 deg VF-179.400+53.345-INT-20220501-VFID1168-Halpha.fits VF-179.321+53.372-INT-20220501-VFID1168-r.fits
CHECK: VFID1304 20220501 20220501 sep=0.110 deg VF-181.400+50.500-INT-20220501-VFID1304-Halpha.fits VF-181.237+50.538-INT-20220501-VFID1304-r.fits
CHECK: VFID3025 20220502 20220505 sep=0.141 deg VF-187.100+28.700-INT-20220502-VFID3025-Halpha.fits VF-186.957+28.764-INT-20220505-VFID3025-r.fits
CHECK: VFID6397 20220502 20220502 sep=0.140 deg VF-215.900+1.700-INT-20220502-VFID6397-Halpha.fits VF-215.763+01.729-INT-20220502-VFID6397-r.fits
CHECK: VFID6183 20220505 20220505 sep=0.079 deg VF-224.518+02.955-INT-20220505-VFID6183-Halpha.fits VF-224.593+02.979-INT-20220505-VFID6183-r.fits
CHECK: VFID4086 20220504 20220504 sep=0.100 deg VF-233.701+16.331-INT-20220504-VFID4086-Halpha.fits VF-233.609+16.377-INT-20220504-VFID4086-r.fits
CHECK: VFID4037 20220502 20220502 sep=0.106 deg VF-234.076+16.541-INT-20220502-VFID4037-Halpha.fits VF-233.977+16.588-INT-20220502-VFID4037-r.fits

Matched 119 INT Ha/r pairs
Wrote int_hybrid_matches.csv
Wrote reproject_jobs.txt

```

### Check that all INT images have a pair

Halpha without a matching r:
```
for h in *INT*-Halpha.fits *INT*-Ha6657.fits; do [[ -e "$h" ]] || continue; r="${h/-Halpha.fits/-r.fits}"; r="${r/-Ha6657.fits/-r.fits}"; [[ ! -e "$r" ]] && echo "$h -> missing $r"; done
```

r without a matching Halpha:
```
for r in *INT*-r.fits; do [[ -e "$r" ]] || continue; h1="${r/-r.fits/-Halpha.fits}"; h2="${r/-r.fits/-Ha6657.fits}"; [[ ! -e "$h1" && ! -e "$h2" ]] && echo "$r -> no Halpha/Ha6657"; done
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
  
  
go into `virgo` conda environment to have access to swarp:

```
conda activate virgo
```

run one:

```bash
python ~/github/HalphaImaging/python3/swarp_int_r_to_old_ha.py --overwrite /data-pool/Halpha/coadds-pre2025-hapy/VF-177.200+56.055-INT-20220502-VFID0957-Halpha.fits /data-pool/Halpha/coadds-pre2025-hapy/VF-177.200+56.055-INT-20220502-VFID0957-Halpha.weight.fits /data-pool/Halpha/coadds-pre2025-hapy/VF-177.200+56.055-INT-20220502-VFID0957-Halpha.fits /data-pool/Halpha/INT-test-hybrid-alignment/VF-177.200+56.055-INT-20220502-VFID0957-Halpha.fits /data-pool/Halpha/INT-test-hybrid-alignment/VF-177.200+56.055-INT-20220502-VFID0957-Halpha.weight.fits
```

```bash
python ~/github/HalphaImaging/python3/swarp_int_r_to_old_ha.py --overwrite  /data-pool/Halpha/coadds-v20260330/VF-177.143+56.083-INT-20220502-VFID0957-r.fits /data-pool/Halpha/coadds-v20260330/VF-177.143+56.083-INT-20220502-VFID0957-r.weight.fits /data-pool/Halpha/coadds-pre2025-hapy/VF-177.200+56.055-INT-20220502-VFID0957-Halpha.fits /data-pool/Halpha/INT-test-hybrid-alignment/VF-177.200+56.055-INT-20220502-VFID0957-r.fits /data-pool/Halpha/INT-test-hybrid-alignment/VF-177.200+56.055-INT-20220502-VFID0957-r.weight.fits
```

```
python ~/github/HalphaImaging/python3/swarp_int_r_to_old_ha.py  --overwrite /data-pool/Halpha/coadds-v20260330/VF-139.119+41.961-INT-20190205-p026-r.fits /data-pool/Halpha/coadds-v20260330/VF-139.119+41.961-INT-20190205-p026-r.weight.fits /data-pool/Halpha/coadds-pre2025-hapy/VF-139.115+41.958-INT-20190205-p026-Halpha.fits /data-pool/Halpha/coadds-v20260518/VF-139.115+41.958-INT-20190205-p026-r.fits /data-pool/Halpha/coadds-v20260518/VF-139.115+41.958-INT-20190205-p026-r.weight.fits
```

```
export OMP_NUM_THREADS=1
parallel --bar -j 1 --joblog swarp_int_hybrid.joblog --results swarp_int_hybrid_logs --colsep ' ' python ~/github/HalphaImaging/python3/swarp_int_r_to_old_ha.py {1} {2} {3} {4} {5} --overwrite :::: reproject_jobs.txt
```
This takes about 10 min.  swarp is already multithreaded, and I've
    found some contamination/scrambling of images when running it in
    parallel, so I just do one at a time.
  
  
```
conda deactivate
```


### Check images
Make a pdf of the halpha and r-band coadds to check that nothing got corrupted:
```
python ~/github/HalphaImaging/python3/make_int_hybrid_coadd_qc.py
```

Output:
```
(virgo) rfinn@draco:/data-pool/Halpha/coadds-v20260609$ python ~/github/HalphaImaging/python3/make_int_hybrid_coadd_qc.py
Found 119 halpha files
Found 119 INT Halpha/r pairs

```

This runs in serial, so this is a good time to refresh your
beverage...

To review the results:
```
(virgo) rfinn@s64247 tables % mkdir coadds-v20260609
(virgo) rfinn@s64247 tables % cd coadds-v20260609 
(virgo) rfinn@s64247 coadds-v20260609 % scp draco:/data-pool/Halpha/coadds-v20260609/int_halpha_r_coadd_qc.pdf .
int_halpha_r_coadd_qc.pdf                                                                                  100% 2498KB   4.9MB/s   00:00    
(virgo) rfinn@s64247 coadds-v20260609 % open int_halpha_r_coadd_qc.pdf 
```
### [OUTDATED - SKIP THIS] 

On the first few passes, the output names for the low dec fields did not contain a leading 0, e.g the dec was `+2.123` instead of `+02.123`.
```
python ~/github/HalphaImaging/python3/rename_low_dec_reprojected_INT_rband.py
```

## 3. Copy new BOK, HDI, and MOS coadds into new coadd directory
so the hybrid dataset becomes self-contained.


The new coadd directory:
```text
/data-pool/coadds-v20260609/
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

## 5. Get filter ratio



## 6. Rebuild INT r-band PSF images

### If the psf images already exist, you can copy the psf info from the psf images:

```
python ~/github/hapy/scripts/copy_psf_header_fields.py --newdir /data-pool/Halpha/coadds-v20260518/ --psfdir /data-pool/Halpha/psf-images-v20260518/ --pattern "*INT*r.fits"
```

### Strategy for assembling new set of psf images

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



# 7. Check Number of coadds

|Instrument | N(R) |N(r) | N(Halpha) |
|---|---|---|
| BOK  |0    | 66 | 66  |
| HDI  | 19 | 12 | 31|
| MOS | 7 |  0 | 7 |
| INT | 0 | 118 | 110 + 8 | 

Total : r 66 + 31 + 12+7 +123 = 222

Total : halpha : = 222

Phew!


# 8. Move on to running-pipeline.md

hapy/docs/running-pipeline.md

