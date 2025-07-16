'''
Takes an nc4 file and reprojects, resamples to match MODIS VI data, calling gdal directly
(espg 4326, pixels size ~ 250m)

using python 3.5.5

script originally by Mike Tercek, last edited by Carolyn Lober 3/12/24
'''

import sys,os,datetime, glob

####
#find "$PWD" |sort  > SAVI_files.txt
#gdalbuildvrt -separate SAVI.vrt -input_file_list /home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/SAVI/SAVI_files.txt

def crop_BLCA(fn, out_fn):
    crop_command  = 'gdal_translate -projwin -107.82987636004 38.6408349291459 -107.618160237389 38.5035082343929 ' + fn +' ' +out_fn
    os.system(crop_command)
    return out_fn

VI_files = sorted(glob.glob('/home/ubuntu/data/results/CONUS/wb_proj/wateryear/*.tif'))
out_files = [x.replace("CONUS","BLCA") for x in VI_files]
run = [crop_BLCA(VI_files[i],out_files[i]) for i in range(0,len(VI_files))]



VI = "EVI"
VI_files = sorted(glob.glob('/home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/'+VI+'/MOD13Q1_'+VI+'_*.tif'))
VI_files = sorted(glob.glob('/home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/Rely/MOD13Q1_Rely_*.tif'))
out_files = [x.replace("CONUS","BLCA") for x in VI_files]

VI_files = sorted(glob.glob('/home/ubuntu/data/results/CONUS/gs_dates/projected/gs_*.tif'))
out_files = [x.replace("CONUS","BLCA") for x in VI_files]
run = [crop_BLCA(VI_files[i],out_files[i]) for i in range(0,len(VI_files))]
