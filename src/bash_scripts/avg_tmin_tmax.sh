#!/bin/bash
#script to calculate average yearly temperature as the mean of tmin and tmax
#call via: bash /home/ubuntu/data/src/bash_scripts/avg_tmin_tmax.sh


for yr in {2000..2023}
do
	gdal_calc.py -A /home/ubuntu/data/results/CONUS/gm_proj/bilinear/tmin/tmmn_${yr}_resampled.tif -B /home/ubuntu/data/results/CONUS/gm_proj/bilinear/tmax/tmmx_${yr}_resampled.tif --outfile=/home/ubuntu/data/results/CONUS/gm_proj/bilinear/tavg/tavg_${yr}_resampled.tif --calc="(A+B)/2"
done
