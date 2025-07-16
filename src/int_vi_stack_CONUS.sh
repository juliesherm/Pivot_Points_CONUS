#!/bin/bash
#script to combine tiles from growing season integrated vi
#call via: source /home/ubuntu/data/src/bash_scripts/int_vi_stack_CONUS.sh

VIS="SAVI EVI"

for vi in $VIS
do
	cd /home/ubuntu/data/results/CONUS/integrated_vi_gs/window_14days/cutoff_2.5/${vi}/tiled
	find "$PWD" |grep ".tif" |sort  > ${vi}_int_files.txt
	gdalbuildvrt ${vi}_int_CONUS.vrt -input_file_list ${vi}_int_files.txt
	gdal_translate ${vi}_int_CONUS.vrt /home/ubuntu/data/results/CONUS/integrated_vi_gs/window_14days/cutoff_2.5/${vi}/${vi}_int_CONUS.tif -co TILED=YES -co COMPRESS=LZW -co COPY_SRC_OVERVIEWS=YES -co BIGTIFF=YES
done