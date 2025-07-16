"""
Created on Fri May 10 2024
Modified May 10, 2024

@author: js
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

VIS = ["EVI", "SAVI"]
WBS = ["accumswe", "aet", "deficit", "pet", "rain", "runoff", "soilwater"]
GMS = ["tmin", "tmax", "tavg", "precip"]

park = "CONUS"
win  = "window_14days/cutoff_2.5/bilinear"
base_path = '/home/ubuntu/data/'
#base_path = '/home/1024/ma/sherman/NPS SIP/pivot_pts/'

pivs_path = base_path+'results/'+park+'/pivot_pts/'+win+'/'

pval_file = pivs_path + VIS[0]+"_"+WBS[0]+"_pval.tif"
sig_rast = rioxarray.open_rasterio(pval_file)
sig_rast = xr.where(sig_rast==-32767, np.nan, sig_rast)
sig_rast = xr.where(sig_rast==32767, np.nan, sig_rast)
sig_rast = xr.where(np.isnan(sig_rast), np.nan, 0)

for vi in VIS:
  for wb in WBS:
    #Open VI data
    pval_file = pivs_path + vi+"_"+wb+"_pval.tif"
    pivs_rast = rioxarray.open_rasterio(pval_file)    
    pivs_rast = xr.where(pivs_rast==-32767, np.nan, pivs_rast)
    pivs_rast = xr.where(pivs_rast==32767, np.nan, pivs_rast)
    pivs_rast = xr.where(pivs_rast<0.1, 1, 0)
    sig_rast = sig_rast+pivs_rast
    #what percent of pixels have a significant relationship for given variable
    print(vi)
    print(wb)
    print(np.nanmean(pivs_rast))

sig_rast = xr.where(sig_rast>0, 1, 0)
#what percent of pixels have a significant relationship for ANY variable
print(np.nanmean(sig_rast))
#0.4363654845846714

'''
### plot??
#pyplot.imshow(sig_rast, cmap='pink')
from rasterio.plot import show
import matplotlib.pyplot as plt
rasterio.plot(sig_rast)
plt.savefig(base_path+'test.png')
'''


'''
EVI:
WB:
0.0952915731119564 - accumswe
0.30103088419824164 - aet
0.3125913862361464 - deficit
0.21086436533759945 - pet
0.29288041425590056 - rain
0.2284634777690521 - runoff
0.28943566384477754 - soilwater
GM:
0.10591111964736677 - tmin
0.22549179502934844 - tmax
0.18433815543385076 - tavg
0.3035498802082998 - precip

SAVI:
WB:
0.09746315109311426
0.30911179628453156
0.31876504625551966
0.21870810336900504
0.2993811315665557
0.2338031557122893
0.2966256080345865
GM:
0.110541755486028
0.23216578570235263
0.19129510782894607
0.30991999785835955

Total:
WB:
0.4363654845846714
GM:
0.5217418297129253
'''