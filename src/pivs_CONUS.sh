#!/bin/bash

# band 4 = rsquared
# band 5 = p-value
#['y-intercept', 'slope', 'x-intercept','rsqaured','pvalue','stderr','n']

VIS="EVI SAVI"
#WBS="accumswe aet deficit pet runoff rain soilwater tmin tmax tavg precip"
WBS="precip"

cd /home/ubuntu/data/results/CONUS/pivot_pts/window_14days/cutoff_2.5/bilinear/tiled

for vi in $VIS
do
	for wb in $WBS
	do
		find "$PWD" |grep "${vi}_${wb}" |sort  > ${vi}_${wb}_piv_files.txt
        	echo $vi
        	gdalbuildvrt ${vi}_${wb}_pval_CONUS.vrt -input_file_list ${vi}_${wb}_piv_files.txt -b 4
        	gdal_translate ${vi}_${wb}_pval_CONUS.vrt ${vi}_${wb}_rsq.tif -co TILED=YES -co COMPRESS=LZW -co COPY_SRC_OVERVIEWS=YES -co BIGTIFF=YES
		rm ${vi}_${wb}_pval_CONUS.vrt
	done
done