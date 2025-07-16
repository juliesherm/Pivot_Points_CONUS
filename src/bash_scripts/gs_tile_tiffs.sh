#!/bin/bash
#script to create tiles from growing season dates
#call via: source /home/ubuntu/data/src/bash_scripts/avg_tmin_tmax.sh

xtilesize=10110
ytilesize=2611
height=13054
width=30329

echo $(seq 0 $height $ytilesize)

for yr in {2000..2023}
do
	for yoff in $(seq 0 $ytilesize $height)
	do
		for xoff in $(seq 0 $xtilesize $width)
		do
			gdal_translate -srcwin $xoff $yoff $xtilesize $ytilesize /home/ubuntu/data/results/CONUS/gs_dates/window_5days/projected/gs_${yr}_resampled.tif /home/ubuntu/data/results/CONUS/gs_dates/window_5days/projected/tiled/gs_${yr}_${xoff}_$yoff.tif
		done
	done
done
