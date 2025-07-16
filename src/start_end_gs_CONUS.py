#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 18:57:25 2024
Modified Mar 14

@author: mt, edits by js
"""
# Demo code to show a vectorized method for finding start and end of growing season, marked as condition AET > 0 in a 3D array with time axis as first dimension and 2D spatial array as next two axes.
# This is meant to be run on daily AET netCDF files with from the NPS water balance dataset.
# Uses numpy vector operations that go straight to C instead of python loops. The entire annual data file would be read in as one step before this processing.
# IMPORTANT. Python indexing begins at zero, unlike R, which begins at 1. R is weird.

#js: Edits include using dask arrays to apply to all of CONUS
#js: output is PYTHON index of growing season gs_dates
#js: the first day is not checked for growing season completeness. Could do this by adding zero at the beginning of the cumulative_sums array and changing some of the indices, but its overlooked for now.
#js: saved as csv files for convenience, combined into geotif using supplementary R file
#js: for loop certainly not optimized. Each year takes approx 7 mins... could definitely improve

import xarray as xr
import numpy as np
import dask.array as da

for x in range(2022,2024):
    print(x)
    fn = '/home/ubuntu/data/raw_data/WB/daily/aet/v_1_5_'+str(x)+'_gridmet_historical_aet.nc4'
    a1 = xr.open_dataset(fn, drop_variables = "tbnds", chunks = "auto")
    a1 = a1.drop_dims("tbnds")
    a = a1.variables["aet"]
    a = da.where(a > 0.25, 1,0) # change array values to just one or zero, all values greater than zero are now 1
    a = a.astype("int16")
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # MODIFICATION TO SCREEN EARLY SEASON ANOMALIES
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------
    # The code above works fine if you don't care about catching early season spurious values.
    # For example in [0,2] of y0 above there is the value 0.5 that does not show up in later days
    # y4 has another 0.5 at [0,2] but then as days progress [0,2] goes back to zero.
    # We really want to ignore these early 0.5 values and find and only find the first index of non-zero valeus that are recurring at [0,2], which happens starting in y13
    # So... starting after the line defining the array a above, i.e. starting atfer a = np.where(a > 0, 1,0) add this:
    time_window = 14 # Number of days in a row you want nonzero values to be consistent
    cumulative_sums = da.cumsum(a, axis = 0) # take running sums on time axis, makes a 3D array, with 2D slices of 2 x 20, as time increases cell values show running total of nonzero values to date
    differences = cumulative_sums[time_window:,:,:] - cumulative_sums[:-time_window,:,:] # Time shift the 2D slices by t = time_window and find rate of increase from one time period to next for each cell.
    differences = differences / time_window # # Think of differences now as a rolling mean of the rate of increase in from one T-day time window to next. Since all values are either 0 or 1, anything less than 1 in this variable indicates at least one zero in time_window of length t.
    differences = da.where(differences < 1, 0, differences) # Replace all the difference values < 1 with 0 and keep only the one values. Now the first 1 values that we encounter in time will be from time windows that contained all 1 with no zeros.
    # Now calculate first and last indices using original methods from above but applied to the new difference array.
    new_first_inds = da.argmax(differences, axis = 0) + 1 # See not below on why need to add 1.
    new_reversed_array = differences[::-1,:,:]
    new_last_inds = da.argmax(new_reversed_array, axis = 0)
    new_last_inds = differences.shape[0] - new_last_inds + time_window - 1 # See note below on why add this offset.
    # The differencing step makes the difference array shorter than the original array a on the time axis by t = timewindow amount.
    # So the indices need to be adjusted here to reflect the shorter length of the array that is applied to argmax.
    # However since we are also taking the average differences and screening out values < 1, the actual needed offset is + 1 for first index and time_window - 2 for last index
    # This is because we lose length time window, but then screen out time_window -1 when we take values below 1 in the average.
    # The final products would be two 2D geotiffs. One has new_first_inds, one has new_last_inds
    fnf = "/home/ubuntu/data/results/CONUS/gs_dates/window_14days/cutoff_2.5/first_inds_"+str(x)+".csv"
    fnl = "/home/ubuntu/data/results/CONUS/gs_dates/window_14days/cutoff_2.5/last_inds_"+str(x)+".csv"
    #a2.to_netcdf(fnf)
    np.savetxt(fnf,np.array(new_first_inds.compute()), delimiter = ",")
    np.savetxt(fnl,np.array(new_last_inds.compute()), delimiter = ",")
