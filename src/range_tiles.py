"""
Takes tiled water balance or gridmet data and gets x min and x max, returning the range
(domain) of water balance for each pixel

Utilizes multiple cores via parallelization
!! WARNING !! Parallelization doesn't work on AWS (yet) - instead use manual loops at the bottom of the script

Created on Mon Mar 18 2024
Modified June 12, 2024
author: js, cl
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

def fast_pivs(wb_wy):
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
    
    wb_rast = rioxarray.open_rasterio(wb_wy)
    
    print("Imputing nans")
    wb_rast = xr.where(wb_rast==-32767, np.nan, wb_rast)
    
    print("Scaling")
    wb_rast = wb_rast/10
    
    print("Calculating min and max")
    minx = np.nanmin(wb_rast, axis = 0)
    maxx = np.nanmax(vi_rast, axis = 0) #should be 0
    
    result = np.stack((maxx - minx), axis = 0)
    
    #Save file
    print("Saving")
    output_file = os.path.join(pivs_tile_path, 'range_{wb}_{xoff}_{yoff}.tif'.format(wb = wb_var, xoff = xoff, yoff = yoff))
    with rasterio.open(output_file, 'w', **meta) as dst:
        dst.descriptions = tuple(['range'])
        dst.write(result)
    en = time.time()
    print(en-st)
    return 0


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
            wb_wy = sorted(glob.glob(os.path.join(wb_tile_path, '*_{xoff}_{yoff}.tif'.format(xoff = xoff, yoff = yoff))))
            if "1024/ma/sherman" in base_path: #can use parallelization on the math servers
                pool = multiprocessing.Pool()
                results = pool.map(fast_pivs, wb_wy)
                pool.close()
                
