"""
Takes tiled waterbalance or gridmet data and regresses onto integrated, centered vegetation indicies for every pixel
Returns tiled pivot points and other values from the simple linear regression

Computes quickly using vectorized functions - based on https://stackoverflow.com/questions/52108417/how-to-apply-linear-regression-to-every-pixel-in-a-large-multi-dimensional-array
Utilizes multiple cores via parallelization
!! WARNING !! Parallelization doesn't work on AWS (yet) - instead use manual loops at the bottom of the script

Created on Mon Mar 18 2024
Modified April 3, 2024
author: js
"""

import numpy as np
import xarray as xr
import glob
import os
import rioxarray
import rasterio
import time
import multiprocessing
from scipy.stats import t
from itertools import repeat

def fast_pivs(vi_gs, wb_wy):
    #pass in file names of growing season integrated VI and wateryear summed waterbalance variable
    if "wb_proj"  in wb_wy:
        wb_var = wb_wy.split("_")[-4] #for waterbalance variables
    elif  "gm_proj" in wb_wy:
        wb_var = wb_wy.split("/")[-1].split("_")[0] #for gridmet variables
    print(wb_var)
    st = time.time()
    #Use the first file as a template to read metadata
    with rasterio.open(wb_wy) as src:
        meta = src.meta
    meta.update({"dtype": 'float64', "count":7})
    #Open VI data
    vi_rast = rioxarray.open_rasterio(vi_gs)
    wb_rast = rioxarray.open_rasterio(wb_wy)
    if wb_rast.shape != vi_rast.shape:
        raise ValueError("vi and wb rasters must have the same dimensions")
    print("Imputing nans")
    vi_rast.values = np.where(wb_rast.values==-32767, np.nan, vi_rast.values)
    wb_rast = xr.where(wb_rast==-32767, np.nan, wb_rast)
    num_nans = np.broadcast_to(np.isnan(vi_rast).sum(axis = 0), vi_rast.shape)
    vi_rast = xr.where(num_nans>20, np.nan,vi_rast)
    #make sure vi and wb have the same nans
    wb_rast.values = np.where(np.isnan(vi_rast.values), np.nan, wb_rast.values)
    del num_nans
    print("Centering")
    cell_means = np.nanmean(vi_rast, axis = 0)
    cell_means = np.broadcast_to(cell_means, vi_rast.shape)
    vi_rast = vi_rast - cell_means
    del cell_means
    print("Scaling")
    vi_rast = vi_rast/10000
    wb_rast = wb_rast/10
    print("Calculating pivot points")
    #fast coefficient estimates
    n = np.count_nonzero(np.logical_not(np.isnan(wb_rast.values)) & np.logical_not(np.isnan(vi_rast.values)), axis = 0)
    sdx = np.nanstd(wb_rast, axis = 0)
    sdy = np.nanstd(vi_rast, axis = 0)
    meanx = np.nanmean(wb_rast, axis = 0)
    meany = np.nanmean(vi_rast, axis = 0) #should be 0
    cor = (np.nansum(vi_rast.values*wb_rast.values, axis = 0) - meany*meanx*n)/(sdx*sdy*n)
    rsq = np.power(cor,2)
    slope = cor*sdy/sdx
    intercept = meany - slope*meanx
    xintercept = -intercept/slope #should be meanx
    tstats = cor*np.sqrt(n-2)/np.sqrt(1-np.power(cor,2))
    stderr = slope/tstats
    p_val = (t.sf(np.abs(tstats),n-2))*2
    result = np.stack((intercept, slope, xintercept, rsq, p_val, stderr, n), axis = 0)
    #Save file
    print("Saving")
    output_file = os.path.join(pivs_tile_path, 'pivs_{vi}_{wb}_{xoff}_{yoff}.tif'.format(vi = VI,wb = wb_var, xoff = xoff, yoff = yoff))
    with rasterio.open(output_file, 'w', **meta) as dst:
        dst.descriptions = tuple(['y-intercept', 'slope', 'x-intercept','rsqaured','pvalue','stderr','n'])
        dst.write(result)
    en = time.time()
    print(en-st)
    return 0


VI = "SAVI"
park = "CONUS"
win  = "window_14days/cutoff_2.5/bilinear"
base_path = '/home/ubuntu/data/'
base_path = '/home/1024/ma/sherman/NPS_SIP/pivot_pts/'

vi_tile_path = base_path+'results/'+park+'/integrated_vi_gs/'+win+'/'+VI+'/tiled/'
wb_tile_path = base_path+'results/'+park+'/wb_proj/wateryear/tiled/'
wb_tile_path = base_path+'results/'+park+'/gm_proj/bilinear/tiled/'
pivs_tile_path = base_path+'results/'+park+'/pivot_pts/'+win+'/tiled/'

#hard code the tile coordinates
xoffs = range(0,30330,10110)
yoffs = range(0,13055,2611)

for xoff in xoffs:
    for yoff in yoffs:
            print(xoff)
            print(yoff)
            vi_gs = vi_tile_path+VI+"_int_"+str(xoff)+"_"+str(yoff)+".tif"
            wb_wy = sorted(glob.glob(os.path.join(wb_tile_path, '*_{xoff}_{yoff}.tif'.format(xoff = xoff, yoff = yoff))))
            if "1024/ma/sherman" in base_path: #can use parallelization on the math servers
                pool = multiprocessing.Pool()
                results = pool.starmap(fast_pivs, zip(repeat(vi_gs),wb_wy))
                pool.close()
            elif "ubuntu/data" in base_path: #For using AWS, use simple loop instead of parallelization
                for w in range(0,len(wb_wy)):
                    #pass in file names of growing season integrated VI and wateryear summed waterbalance variable
                    if "wb_proj"  in wb_wy:
                        wb_var = wb_wy.split("_")[-4] #for waterbalance variables
                    elif  "gm_proj" in wb_wy:
                        wb_var = wb_wy.split("/")[-1].split("_")[0] #for gridmet variables
                    print(wb_var)
                    st = time.time()
                    #Use the first file as a template to read metadata
                    with rasterio.open(wb_wy) as src:
                        meta = src.meta
                    meta.update({"dtype": 'float64', "count":7})
                    #Open VI data
                    vi_rast = rioxarray.open_rasterio(vi_gs)
                    wb_rast = rioxarray.open_rasterio(wb_wy)
                    if wb_rast.shape != vi_rast.shape:
                        raise ValueError("vi and wb rasters must have the same dimensions")
                    print("Imputing nans")
                    vi_rast.values = np.where(wb_rast.values==-32767, np.nan, vi_rast.values)
                    wb_rast = xr.where(wb_rast==-32767, np.nan, wb_rast)
                    num_nans = np.broadcast_to(np.isnan(vi_rast).sum(axis = 0), vi_rast.shape)
                    vi_rast = xr.where(num_nans>20, np.nan,vi_rast)
                    #make sure vi and wb have the same nans
                    wb_rast.values = np.where(np.isnan(vi_rast.values), np.nan, wb_rast.values)
                    del num_nans
                    print("Centering")
                    cell_means = np.nanmean(vi_rast, axis = 0)
                    cell_means = np.broadcast_to(cell_means, vi_rast.shape)
                    vi_rast = vi_rast - cell_means
                    del cell_means
                    print("Scaling")
                    vi_rast = vi_rast/10000
                    wb_rast = wb_rast/10
                    print("Calculating pivot points")
                    #fast coefficient estimates
                    n = np.count_nonzero(np.logical_not(np.isnan(wb_rast.values)) & np.logical_not(np.isnan(vi_rast.values)), axis = 0)
                    sdx = np.nanstd(wb_rast, axis = 0)
                    sdy = np.nanstd(vi_rast, axis = 0)
                    meanx = np.nanmean(wb_rast, axis = 0)
                    meany = np.nanmean(vi_rast, axis = 0) #should be 0
                    cor = (np.nansum(vi_rast.values*wb_rast.values, axis = 0) - meany*meanx*n)/(sdx*sdy*n)
                    rsq = np.power(cor,2)
                    slope = cor*sdy/sdx
                    intercept = meany - slope*meanx
                    xintercept = -intercept/slope #should be meanx
                    tstats = cor*np.sqrt(n-2)/np.sqrt(1-np.power(cor,2))
                    stderr = slope/tstats
                    p_val = (t.sf(np.abs(tstats),n-2))*2
                    result = np.stack((intercept, slope, xintercept, rsq, p_val, stderr, n), axis = 0)
                    #Save file
                    print("Saving")
                    output_file = os.path.join(pivs_tile_path, 'pivs_{vi}_{wb}_{xoff}_{yoff}.tif'.format(vi = VI,wb = wb_var, xoff = xoff, yoff = yoff))
                    with rasterio.open(output_file, 'w', **meta) as dst:
                        dst.descriptions = tuple(['y-intercept', 'slope', 'x-intercept','rsqaured','pvalue','stderr','n'])
                        dst.write(result)
                    en = time.time()
                    print(en-st)