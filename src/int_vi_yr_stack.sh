#!/bin/bash
#script to combine years from growing season integrated vi
#call via: source /home/ubuntu/data/src/bash_scripts/int_vi_yr_stack.sh

VIS="SAVI EVI"

xtilesize=10110
ytilesize=2611
height=13054
width=30329


for vi in $VIS
do
	cd /home/ubuntu/data/results/CONUS/integrated_vi_gs/window_14days/cutoff_2.5/bilinear/${vi}/tiled
	for yoff in $(seq 0 $ytilesize $height)
	do
		for xoff in $(seq 0 $xtilesize $width)
		do
			find "$PWD" |grep "_${xoff}_${yoff}" |sort  > ${vi}_int_${xoff}_${yoff}_files.txt
			gdalbuildvrt -separate ${vi}_${xoff}_${yoff}.vrt -input_file_list /home/ubuntu/data/results/CONUS/integrated_vi_gs/window_14days/cutoff_2.5/bilinear/${vi}/tiled/${vi}_int_${xoff}_${yoff}_files.txt
			gdal_translate /home/ubuntu/data/results/CONUS/integrated_vi_gs/window_14days/cutoff_2.5/bilinear/${vi}/tiled/${vi}_${xoff}_${yoff}.vrt /home/ubuntu/data/results/CONUS/integrated_vi_gs/window_14days/cutoff_2.5/bilinear/${vi}/tiled/${vi}_int_${xoff}_${yoff}.tif -co TILED=YES -co COMPRESS=LZW -co COPY_SRC_OVERVIEWS=YES -co BIGTIFF=YES
		done
	done
done