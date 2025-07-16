#!/bin/bash
#script to compare future tiles to calculated pivot points
#call via: bash /home/1024/ma/sherman/NPS_SIP/pivot_pts/src/bash_scripts/future_xint_comp.sh

WBS="aet deficit pet accumswe soilwater"
MODS="BNU-ESM HadGEM2-CC365 NorESM1-M  IPSL-CM5A-LR GFDL-ESM2G MIROC-ESM-CHEM CanESM2 MIROC5 CSIRO-Mk3-6-0 inmcmc4 CNRM-CM5 CCSM4 MRI-CGCM3"
RCPS="rcp85"
#WBS="accumswe aet deficit pet runoff rain soilwater tmin tmax tavg precip"

#don't need to do for both VI's because the pivot points are the same ( = average wb )
for rcp in $RCPS
do
    for mod in $MODS
    do
        for wb in $WBS
        do
            gdal_calc.py --calc="A>B"  -A /home/1024/ma/sherman/NPS_SIP/pivot_pts/results/CONUS/pivot_pts/window_14days/cutoff_2.5/bilinear/EVI_${wb}_xint.tif -B /home/1024/ma/sherman/NPS_SIP/pivot_pts/results/CONUS/wb_proj/futures/V_1_5_annual_${mod}_${rcp}_${wb}_resampled.tif --allBands B --outfile=/home/1024/ma/sherman/NPS_SIP/pivot_pts/results/CONUS/future_proj/comp_xint_${mod}_${rcp}_${wb}.tif 
        done
    done
done