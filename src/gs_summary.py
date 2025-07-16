#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 08:52 2024

@author: cl
"""
import numpy as np
import os
import glob
import rasterio
import xarray as xr
import rioxarray
import xrscipy.signal as dsp
import dask.array as da
import datetime
from scipy.signal import savgol_filter
from pyprojroot.here import here
from dask.diagnostics import ProgressBar

def integrate_smooth_gs(vals, dates):
    
    if vals.size == 0:
        raise ValueError("vals is length 0")
    else:
        print(vals.size)
        print(type(vals))
    
    # get gs dates from vals
    gs_start = int(vals[0]); vals = vals[1:]
    gs_end = int(vals[0]); vals = vals[1:]
    
    vals = np.array(vals)
    dates = np.array(dates)
    
    # smooth with savgol
    vals_smth = savgol_filter(vals, window_length = 5, polyorder = 1, deriv = 0)
    
    # check that x and y have the same length
    if len(vals) != len(dates):
        raise ValueError("x and y must have the same length.")
    
    # get dates of MODIS pixels in growing season and add start and end dates
    gs_dates = np.concatenate([[gs_start], dates[(dates > gs_start) * (dates < gs_end)], [gs_end]])
    
    # get VI vals for gs_dates with linear interpolation for the start and end dates
    gs_vals = np.interp(gs_dates, dates, vals)
    
    # calculate the integral using trapezoids
    gs_integral = np.trapz(gs_vals, gs_dates)
    
    return gs_integral

var = "EVI"

path_in = here('raw_data/MODIS/CONUS/VI_16Days_250m_v61/{var}/'.format(var = var))
path_rely = here('raw_data/MODIS/CONUS/VI_16Days_250m_v61/Rely/')
path_gs_dates = here('results/CONUS/gs_dates/projected/')
path_out = here('results/CONUS/gs_summary/{var}/'.format(var = var))

i = 0

for yr in range(2001,2022):
    print(yr)
    print(datetime.datetime.now())
    
    # get MODIS files for the relevant year
    file_pattern = os.path.join(path_in, '*_{yr}_*.tif'.format(yr = yr))
    yr_files = glob.glob(file_pattern)
    vi_yday = [int(yr_file.split('_')[7][:-4]) for yr_file in yr_files]
    rely_pattern = os.path.join(path_rely, '*_{yr}_*.tif'.format(yr = yr))
    rely_files = glob.glob(rely_pattern)
    rely_yday = [int(yr_file.split('_')[7][:-4]) for yr_file in yr_files]
    
    if vi_yday != rely_yday:
        print('Warning: VI and Rely yday do not match')
    
    # get gs_dates file for the relevant year
    gs_pattern = os.path.join(path_gs_dates, '*_{yr}*.tif'.format(yr = yr))
    gs_dates = rioxarray.open_rasterio(glob.glob(gs_pattern)[0])
    
    # load in MODIS tiffs
    vi_da = xr.concat([rioxarray.open_rasterio(i, chunks = "auto") for i in yr_files],
                      dim = xr.Variable('time',vi_yday)).sortby('time')
    rely_da = xr.concat([rioxarray.open_rasterio(i, chunks = "auto") for i in rely_files],
                      dim = xr.Variable('time',rely_yday)).sortby('time')
    
    # get metadata
    with rasterio.open(yr_files[0]) as src:
        meta = src.meta
    
    # to dataset
    vi_ds = vi_da.to_dataset('band').rename({1: 'vi'})
    rely_ds = rely_da.to_dataset('band').rename({1: 'rely'})
    
    # mask reliability (clouds -> NA, winter min -> 0.05)
    vi_masked = da.where(rely_da == 3, np.nan, vi_da)
    vi_masked = da.where(rely_da == 2, 500, vi_masked)
    vi_masked = vi_masked[:,0,:,:] # remove extra dimension
    
    # Create a list of delayed objects
    stacked = da.vstack((gs_dates, vi_masked))
    stacked = stacked.rechunk({25, 5, 30})
    
    result = da.apply_along_axis(integrate_smooth_gs, axis=0, arr=stacked, dates=vi_yday, dtype = np.float64, shape = stacked.shape)
    output_file = os.path.join(path_out, 'gs_summary_{yr}.tif'.format(yr = yr))
    with ProgressBar():
        with rasterio.open(output_file, 'w', **meta) as dst:
            dst.write(result.compute(), 1)
    
    i += 1
    if i > 0:
        break

    
    
    
    
    
    
    
