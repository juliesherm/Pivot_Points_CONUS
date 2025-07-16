"""
Computes mean historical iEVI and saves using tiled data by averaging over bands (years) 
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
from itertools import repeat

def vi_mean(vi_gs):
    #pass in file names of growing season integrated VI
    st = time.time()
    #Use the first file as a template to read metadata
    with rasterio.open(vi_gs) as src:
        meta = src.meta
    meta.update({"dtype": 'int64', "count":1, "compress":"LZW"})
    #Open VI data
    vi_rast = rioxarray.open_rasterio(vi_gs)
    print("Averaging")
    cell_means = np.nanmean(vi_rast, axis = 0)
    #print("Scaling")
    #cell_means = cell_means/10000
    #Save file
    print("Saving")
    output_file = os.path.join(vi_tile_path, 'mean_{vi}_{xoff}_{yoff}.tif'.format(vi = VI, xoff = xoff, yoff = yoff))
    with rasterio.open(output_file, 'w', **meta) as dst:
        dst.descriptions = tuple(['average_int_gs_vi'])
        dst.write(cell_means,1)
    en = time.time()
    print(en-st)
    return 0


VI = "EVI"
park = "CONUS"
win  = "window_14days/cutoff_2.5/bilinear"
base_path = '/home/ubuntu/data/'
base_path = '/home/1024/ma/sherman/NPS_SIP/pivot_pts/'
base_path = '/home/scratch/sherman/'
backup_path = '/home/16t-usb/sherman2/NPS_SIP/pivot_pts/'


vi_tile_path = base_path+'results/'+park+'/integrated_vi_gs/'+win+'/'+VI+'/tiled/'
vi_tile_path = base_path+'results/'+park+'/integrated_vi_gs/'
#hard code the tile coordinates
xoffs = range(0,30330,10110)
yoffs = range(0,13055,2611)

for xoff in xoffs:
    for yoff in yoffs:
            print(xoff)
            print(yoff)
            vi_gs = vi_tile_path+VI+"_int_"+str(xoff)+"_"+str(yoff)+".tif"
            if "1024/ma/sherman" in base_path: #can use parallelization on the math servers
                #pool = multiprocessing.Pool()
                #results = pool.starmap(vi_mean, zip(repeat(vi_gs)))
                #pool.close()
                print("Don't work here")
            elif "ubuntu/data" in base_path or "scratch/sherman" in base_path: #For using AWS, use simple loop instead of parallelization
                vi_mean(vi_gs)
