'''
This script is used for un-tiled file structures, such as for specific parks
Otherwise, it is the same as vi_smooth_tiles
From raw MODIS VI data for each year, processes and computes the integrated VI value over the growing season
 - masks with MODIS reliability data
 - substitutes given value (500) over snow
 - smooths with savgol filter
 - integrates over the growing season
   - median of calculated growing seasons
   - gs defined by aet>cutoff for "win" consecutive days
 - saves result
'''
import numpy as np
import xarray as xr
import glob
import os
import rioxarray
import rasterio
import time
import gc
from scipy.signal import savgol_filter

def integrate_gs(vi_vals, dates):
    gs_bounds = vi_vals[0:2]
    vi_vals = vi_vals[2:]
    dates = np.array(dates)[~np.isnan(vi_vals)]
    if len(dates)==0:
        return np.nan
    gs_dates = np.concatenate(([gs_bounds[0]],dates,[gs_bounds[1]]))
    # get VI vals for gs_dates with linear interpolation for the start and end dates
    gs_vals = np.interp(gs_dates, dates, vi_vals[~np.isnan(vi_vals)])
    del vi_vals, dates, gs_bounds
    # calculate and return the integral using trapezoidsl
    return np.trapz(gs_vals, gs_dates)

VI = "EVI"
park = "BLCA"
win  = "window_14days/cutoff_2.5"
base_path = '/home/ubuntu/data/'
#base_path = '/home/1024/ma/sherman/NPS_SIP/pivot_pts/'


yday_path = base_path+'raw_data/MODIS/'+park+'/VI_16Days_250m_v61/'+VI+'/'
vi_tile_path = base_path+'raw_data/MODIS/'+park+'/VI_16Days_250m_v61/'+VI+'/'
rely_tile_path = base_path+'raw_data/MODIS/'+park+'/VI_16Days_250m_v61/Rely/'
gs_dates_tile_path = base_path+'results/'+park+'/gs_dates/'+win+'/'
gs_dates_median_tile_path = base_path+'results/'+park+'/gs_dates/'+win+'/'
tile_path_out = base_path+'results/'+park+'/integrated_vi_gs/'+win+'/'+VI+'/'

for yr in range(2000,2024):
    print(yr)
    yday_files = glob.glob(os.path.join(yday_path, '*_{yr}_*.tif'.format(yr = yr)))
    #Find the dates of the bands
    vi_yday = sorted([int(yday_file.split('_')[7][:-4]) for yday_file in yday_files])
    #FOR MATH COMPUTERS - SAVE AS TEXT AND READ IN
    #np.savetxt(yday_path+"ydays"+str(yr)+".csv",vi_yday, delimiter =",")    
    #vi_yday = np.loadtxt(yday_path+"ydays"+str(yr)+".csv", delimiter=",")
    #vi_yday = vi_yday.tolist()
    #List all the tiled files for relevant year
    #gs_date_tiles = sorted(glob.glob(os.path.join(gs_dates_tile_path, '*_{yr}_*.tif'.format(yr = yr))))
    vi_tiles = sorted(glob.glob(os.path.join(vi_tile_path, '*_{yr}_*.tif'.format(yr = yr))))
    rely_tiles = sorted(glob.glob(os.path.join(rely_tile_path, '*_{yr}_*.tif'.format(yr = yr))))
    #Use the first file as a template to read metadata
    with rasterio.open(vi_tiles[1]) as src:
        meta = src.meta
    meta.update({"dtype": 'float64', "count":1})
    #Open VI data
    vi_rast = xr.concat([rioxarray.open_rasterio(i) for i in vi_tiles], dim = 'band')
    rely_rast = xr.concat([rioxarray.open_rasterio(i) for i in rely_tiles], dim = 'band')
    if len(vi_yday) != vi_rast.shape[0]:
        raise ValueError("ydays must have the same length as the number off bands in the raster")
    print("Masking")
    #Mask based on Reliability [clouds -> NA]
    vi_rast = xr.where(rely_rast==3,  np.nan,vi_rast)
    #Winter Minimum [snow -> 500]
    vi_rast = xr.where(rely_rast==2,  500, vi_rast)
    del rely_rast
    #No value values
    vi_rast = xr.where(vi_rast<=-3000,  np.nan, vi_rast)
    #Interpolate nan values to avoid errors in Savgol smoother
    vi_rast['band'] = vi_yday
    n = vi_rast.shape[0]
    print("Interpolating over nans")
    vi_rast_nonenan = vi_rast.interpolate_na(dim = 'band', method = "linear",fill_value = "extrapolate")
    del vi_rast
    #Check that any remaining nans are entirely nan cells (or only one non-nan) 
    if len(np.unique(np.sum(np.isnan(vi_rast_nonenan),axis = 0)))>3:
        raise ValueError("Nans exist apart from fully nan cells, this will be an issue in future steps")
    #Save which cells are nan for re-masking out later
    nan_cells = xr.where(np.sum(np.isnan(vi_rast_nonenan),axis = 0)>(n-2), True, False)
    #Fill with -9999, so that savgol doesn't throw an error
    vi_rast_nonenan = vi_rast_nonenan.fillna(-9999) 
    #Smooth using savgolay_filter
    print("Savgol smoothing")
    vi_rast_sm = savgol_filter(vi_rast_nonenan, window_length = 5, polyorder = 1, deriv=0, axis = 0)
    del vi_rast_nonenan
    #Read in and quality check the growing season start and end dates
    #gs_date_rast = rioxarray.open_rasterio(gs_date_tiles[t])
    gs_date_rast = rioxarray.open_rasterio(os.path.join(gs_dates_median_tile_path+ 'gs_medians.tif'))
    if np.array((gs_date_rast[0,]<0).any()):
        raise ValueError("some cell has invalid growing season start date")
    if np.array((gs_date_rast[1,]>366).any()):
        raise ValueError("some cell has invalid growing season end date")
    #Mask the dates outside of the growing season
    vi_yday_arr = np.broadcast_to(np.array(vi_yday)[:,  np.newaxis, np.newaxis],vi_rast_sm.shape)    
    vi_rast_sm = xr.where(vi_yday_arr<np.broadcast_to(gs_date_rast[0],vi_yday_arr.shape), np.nan, vi_rast_sm)
    vi_rast_sm = xr.where(vi_yday_arr>np.broadcast_to(gs_date_rast[1],vi_yday_arr.shape), np.nan, vi_rast_sm)
    stacked = np.concatenate((gs_date_rast,vi_rast_sm), axis = 0)
    del gs_date_rast, vi_rast_sm
    #Integrate over growing season
    print("Integrating")
    result = np.apply_along_axis(integrate_gs, axis = 0, arr = stacked, dates = vi_yday)
    #Mask back out the fully-nan cells
    result = np.where(nan_cells, np.nan, result)
    #Save file
    print("Saving")
    output_file = os.path.join(tile_path_out, 'gs_summary_{yr}'.format(yr = yr))
    with rasterio.open(output_file, 'w', **meta) as dst:
        dst.write(result,1)
    del stacked, result
    en = time.time()
    print("finished a year!")
    gc.collect()
