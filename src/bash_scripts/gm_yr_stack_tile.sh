#!/bin/bash
#script to combine yearly gridmet data into single multibanded tif, then tile into 15 separate files
#call via: source /home/ubuntu/data/src/bash_scripts/gm_yr_stack_tile.sh

GMS="precip tmin tmax tavg"

xtilesize=10110
ytilesize=2611
height=13054
width=30329


for vi in $GMS
do
	cd /home/ubuntu/data/results/CONUS/gm_proj/bilinear/${vi}
	find "$PWD" |grep ".tif$" |sort  > ${vi}_files.txt
	gdalbuildvrt -separate ${vi}.vrt -input_file_list /home/ubuntu/data/results/CONUS/gm_proj/bilinear/${vi}/${vi}_files.txt
	for yoff in $(seq 0 $ytilesize $height)
	do
		for xoff in $(seq 0 $xtilesize $width)
		do
			gdal_translate -srcwin $xoff $yoff $xtilesize $ytilesize /home/ubuntu/data/results/CONUS/gm_proj/bilinear/${vi}/${vi}.vrt /home/ubuntu/data/results/CONUS/gm_proj/bilinear/tiled/${vi}_${xoff}_${yoff}.tif -co TILED=YES -co COMPRESS=LZW -co COPY_SRC_OVERVIEWS=YES -co BIGTIFF=YES
		done
	done
done