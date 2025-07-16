#!/bin/bash
#script to combine raw MODIS data into  multibanded tif by year, then tile into 15 separate files
#call via: source /home/ubuntu/data/src/bash_scripts/modis_yr_stack_tile.sh


VIS="SAVI EVI  Rely"
xtilesize=10110
ytilesize=2611
height=13054
width=30329

for vi in $VIS
do
	cd /home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/{$vi}
	for yrs in {2000..2023}
	do
		find "$PWD" |grep "${yrs}" |sort  > ${vi}_${yrs}_files.txt
		gdalbuildvrt -separate ${vi}${yrs}.vrt -input_file_list /home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/${vi}/${vi}_${yrs}_files.txt
		#stack MODIS files per year and save as multibanded tif
		gdal_translate /home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/$vi/$vi$yrs.vrt /home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/$vi/$vi$yrs.tif -co TILED=YES -co COMPRESS=LZW -co COPY_SRC_OVERVIEWS=YES -co BIGTIFF=YES
		#separate into 15 separate tiles for smaller file sizes
		for yoff in $(seq 0 $ytilesize $height)
		do
			for xoff in $(seq 0 $xtilesize $width)
			do
				gdal_translate -srcwin $xoff $yoff $xtilesize $ytilesize /home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/$vi/$vi$yrs.tif /home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/$vi/tiled/${vi}_${yrs}_${xoff}_$yoff.tif
			done
		done
	done
done