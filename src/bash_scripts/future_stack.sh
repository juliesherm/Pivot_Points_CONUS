#!/bin/bash

#call via: bash /home/1024/ma/sherman/NPS_SIP/pivot_pts/src/bash_scripts/future_stack.sh

# band 1 = total year below pivot
# band 2 = intensity below pivot
# band 3 = longest  consecutive below pivot



WBS="aet deficit accumswe pet runoff soilwater"
#WBS="deficit"
RCPS="rcp45 rcp85"
#RCPS="rcp85"

#cd /home/1024/ma/sherman/NPS_SIP/pivot_pts/results/CONUS/future_comparison
xoff=0
yoff=5222
for rcp in $RCPS
do
	for wb in $WBS
	do
		find "$PWD" |grep "comp_EVI_midcent_sd_.*_${rcp}_${wb}_${xoff}_${yoff}.tif$" |sort  > ${rcp}_${wb}_midfiles.txt
       		gdalbuildvrt -input_file_list ${rcp}_${wb}_midfiles.txt -b 1 -separate ${rcp}_${wb}_sd_midcent_CONUS.vrt
       		gdal_translate ${rcp}_${wb}_sd_midcent_CONUS.vrt ${rcp}_${wb}_EVI_midcent_sd_${xoff}_${yoff}.tif -co TILED=YES -co COMPRESS=LZW -co COPY_SRC_OVERVIEWS=YES -co BIGTIFF=YES
		rm ${rcp}_${wb}_sd_midcent_CONUS.vrt
		rm ${rcp}_${wb}_midfiles.txt
		find "$PWD" |grep "comp_EVI_endcent_sd_.*_${rcp}_${wb}_${xoff}_${yoff}.tif$" |sort  > ${rcp}_${wb}_endfiles.txt
       		gdalbuildvrt -input_file_list ${rcp}_${wb}_endfiles.txt -b 1 -separate ${rcp}_${wb}_sd_endcent_CONUS.vrt
       		gdal_translate ${rcp}_${wb}_sd_endcent_CONUS.vrt ${rcp}_${wb}_EVI_endcent_sd_${xoff}_${yoff}.tif -co TILED=YES -co COMPRESS=LZW -co COPY_SRC_OVERVIEWS=YES -co BIGTIFF=YES
		rm ${rcp}_${wb}_sd_endcent_CONUS.vrt
		rm ${rcp}_${wb}_endfiles.txt

	done
done
