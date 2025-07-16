import rioxarray
import xarray as xr
import  glob
import rasterio
import os

vi_files = sorted(glob.glob('/home/ubuntu/data/raw_data/MODIS/CONUS/VI_16Days_250m_v61/EVI/tiled/*.tif'))

#tiling each years first gs doy and  last gs doy layers
for vi_file in vi_files:
    vi_data = rioxarray.open_rasterio(vi_file)
    with rasterio.open(vi_file) as src:
            meta = src.meta
    meta.update({"dtype": 'int16', "count":2})
    yr = vi_file.split('_')[-3]
    gs_data = rioxarray.open_rasterio('/home/ubuntu/data/results/CONUS/gs_dates/window_5days/gs_'+yr+'.tif')
    gs_data = gs_data.rio.reproject_match(vi_data)
    output_file = '/home/ubuntu/data/results/CONUS/gs_dates/window_5days/tiled/gs_'+yr+'_'+vi_file.split('_')[-2]+'_'+vi_file.split('_')[-1]
    with rasterio.open(output_file, 'w', **meta) as dst:
        dst.write(gs_data)
    


# just tiling the median rasters
park = "CONUS"
vi_files = sorted(glob.glob('/home/ubuntu/data/raw_data/MODIS/'+park+'/VI_16Days_250m_v61/EVI/tiled/EVI_2000_*.tif'))

for vi_file in vi_files:
    vi_data = rioxarray.open_rasterio(vi_file)
    with rasterio.open(vi_file) as src:
            meta = src.meta
    meta.update({"dtype": 'int16', "count":2})
    gs_data = rioxarray.open_rasterio('/home/ubuntu/data/results/CONUS/gs_dates/window_14days/cutoff_2.5/gs_all_doys.tif')
    gs_data = gs_data[-2:,:,:]
    gs_data = gs_data.rio.reproject_match(vi_data, resampling = rasterio.enums.Resampling.cubic)
    output_file = '/home/ubuntu/data/results/'+park+'/gs_dates/window_14days/cutoff_2.5/tiled/cubic/gs_medians_'+vi_file.split('_')[-2]+'_'+vi_file.split('_')[-1]
    with rasterio.open(output_file, 'w', **meta) as dst:
        dst.write(gs_data)
    


# projecting without tiles the median rasters
park = "BLCA"
vi_file = sorted(glob.glob('/home/ubuntu/data/raw_data/MODIS/'+park+'/VI_16Days_250m_v61/EVI/*.tif'))[0]
vi_data = rioxarray.open_rasterio(vi_file)
with rasterio.open(vi_file) as src:
        meta = src.meta
meta.update({"dtype": 'int16', "count":2})
gs_data = rioxarray.open_rasterio('/home/ubuntu/data/results/CONUS/gs_dates/window_14days/cutoff_2.5/gs_all_doys.tif')
gs_data = gs_data[-2:,:,:]
gs_data = gs_data.rio.reproject_match(vi_data)
output_file = '/home/ubuntu/data/results/'+park+'/gs_dates/window_14days/cutoff_2.5/gs_medians.tif'
with rasterio.open(output_file, 'w', **meta) as dst:
        dst.write(gs_data)
        


import rioxarray
import xarray as xr
import  glob
import rasterio
import os

#tiling the (reprojected) wateryear waterbalance variables
park = "CONUS"
base_path = '/home/1024/ma/sherman/NPS SIP/pivot_pts/'
#base_path = '/home/ubuntu/data/'
vi_files = sorted(glob.glob(base_path+'raw_data/MODIS/'+park+'/VI_16Days_250m_v61/SAVI/tiled/SAVI_2000_*.tif'))
wb_files = glob.glob(base_path+'results/'+park+'/wb_proj/wateryear/bilinear/V_1_5_*.tif')
wb_vars = [x.split('_')[8].split('.')[0] for x in wb_files]
for wb_var in wb_vars:
        print(wb_var)
        for vi_file in vi_files:
                print(vi_file)
                vi_data = rioxarray.open_rasterio(vi_file)
                with rasterio.open(vi_file) as src:
                        meta = src.meta
                meta.update({'dtype':'int16', 'count':24})
                yr = vi_file.split('_')[-3]
                output_file = base_path+'results/CONUS/wb_proj/wateryear/bilinear/tiled/V_1_5_wateryears_gridmet_historical_'+wb_var+'_resampled'+'_'+vi_file.split('_')[-2]+'_'+vi_file.split('_')[-1]
                if os.path.isfile(output_file):
                        print("Already exists")
                        continue
                wb_data = rioxarray.open_rasterio(base_path+'results/CONUS/wb_proj/wateryear/bilinear/V_1_5_wateryears_gridmet_historical_'+wb_var+'_resampled.tif')
                # just tile the years corresponding to MODIS data
                wb_data = wb_data[-24:,:,:]
                wb_data = wb_data.rio.reproject_match(vi_data, resampling = rasterio.enums.Resampling.bilinear)
                with rasterio.open(output_file, 'w', **meta) as dst:
                        dst.write(wb_data)
        



#tiling the (future) wateryear waterbalance variables
park = "CONUS"
base_path = '/home/1024/ma/sherman/NPS SIP/pivot_pts/'
#base_path = '/home/ubuntu/data/'
vi_files = sorted(glob.glob(base_path+'raw_data/MODIS/'+park+'/VI_16Days_250m_v61/SAVI/tiled/SAVI_2000_*.tif'))
wb_files = glob.glob(base_path+'raw_data/WB/futures/bilinear/V_1_5_*.tif')
wb_vars = np.unique([x.split('_')[-1].split('.')[0] for x in wb_files])
rcps = np.unique([x.split('_')[-2] for x in wb_files])
mods = np.unique([x.split('_')[-3] for x in wb_files])
for wb_var in wb_vars:
        print(wb_var)
        for vi_file in vi_files:
                print(vi_file)
                vi_data = rioxarray.open_rasterio(vi_file)
                with rasterio.open(vi_file) as src:
                        meta = src.meta
                meta.update({'dtype':'int16', 'count':24})
                yr = vi_file.split('_')[-3]
                output_file = base_path+'results/CONUS/wb_proj/futures/bilinear/tiled/V_1_5_'+wb_var+'_resampled'+'_'+vi_file.split('_')[-2]+'_'+vi_file.split('_')[-1]
                if os.path.isfile(output_file):
                        print("Already exists")
                        continue
                wb_data = rioxarray.open_rasterio(base_path+'raw_data/WB/futures/annualnetcdfs/V_1_5_wateryears_gridmet_historical_'+wb_var+'_resampled.tif')
                # just tile the years corresponding to MODIS data
                wb_data = wb_data.rio.reproject_match(vi_data, resampling = rasterio.enums.Resampling.bilinear)
                with rasterio.open(output_file, 'w', **meta) as dst:
                        dst.write(wb_data)

