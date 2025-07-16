'''
Designed for tiled CONUS future projection data (use future_tile.sh to create from resampled WB data)
From yearly future data for model, compare to the calculated pivot point (x-intercept)
 - calculate total years below (above - deficit) pivot 
 - calculate intensity (above - deficit) pivot
 - calculate longest stretch of years below (above - deficit) pivot
 - calculate the last year below pivot
 - saves result
'''

import numpy as np
import xarray as xr
import glob
import os
import rioxarray
import rasterio
import time
import  gc
import multiprocessing
from itertools import repeat

def comp_mod_piv(mod_file, metric):
    params = mod_file.split("_")
    mod = params[-5]
    rcp = params[-4]
    wb = params[-3]
    yoff = params[-1].split(".")[0]
    xoff = params[-2]
    piv_file = os.path.join(pivs_path, 'EVI_{wb}_{xoff}_{yoff}_xint.tif'.format(wb = wb,xoff=xoff,yoff=yoff))
    piv_file = os.path.join(pivs_path, 'pivs_EVI_{wb}_{xoff}_{yoff}.tif'.format(wb = wb,xoff=xoff,yoff=yoff))
    out_path = base_path+'results/'+park+'/future_comparison/tiled/'+metric+'/'
    mid_output_file = os.path.join(out_path, 'comp_EVI_midcent_{metric}_{mod}_{rcp}_{wb}_{xoff}_{yoff}.tif'.format(mod=mod,metric = metric, rcp=rcp, wb=wb,xoff=xoff,yoff=yoff))
    end_output_file = os.path.join(out_path, 'comp_EVI_endcent_{metric}_{mod}_{rcp}_{wb}_{xoff}_{yoff}.tif'.format(mod=mod,metric = metric, rcp=rcp, wb=wb,xoff=xoff,yoff=yoff))
    if os.path.exists(end_output_file):
        print("Already completed this one!")
        #print(output_file)
        return 0
    st = time.time()
    #Use the first file as a template to read metadata
    with rasterio.open(mod_file) as src:
        meta = src.meta
    #Open future projected data
    mod_rast = rioxarray.open_rasterio(mod_file)
    if mod_rast.shape[0]>70:
        mod_rast = mod_rast[-70:,:,:]
    #Open pivot points
    piv_rast = rioxarray.open_rasterio(piv_file)
    slp = piv_rast[1,:,:]
    xint = piv_rast[2,:,:]
    yint = piv_rast[0,:,:]
    #pval = piv_rast[4,:,:]
    #Save nans
    print("Finding nans")
    nan_cells = xr.where(np.isnan(xint),np.nan,0)
    #TO DO
    #Open p-values
    #Mask by significant p-values
    #print("Comparing to x-intercepts")
    piv_rast = np.broadcast_to(xint*10,mod_rast.shape)
    #Mask based on above/below pivot point
    #True is "less water",  False is "more water"
    if wb=="deficit":
        mod_rast_less = piv_rast<mod_rast
    else:
        mod_rast_less = piv_rast>mod_rast
    if  metric=='last':
        last_output_file = os.path.join(out_path, 'last_EVI_{mod}_{rcp}_{wb}_{xoff}_{yoff}.tif'.format(mod=mod, rcp=rcp, wb=wb,xoff=xoff,yoff=yoff))
        ######################## LAST GOOD WATER YEAR ##########################
        #year after which all following years are "less water" years
        #find  the last False
        mod_rast_less = ~mod_rast_less[::-1] #convert problem to finding  first "True"
        #last_result = 2030 + mod_rast_less.shape[0]- mod_rast_less.argmax('time') - 1
        last_result = 2030 + mod_rast_less.shape[0]- mod_rast_less.argmax(mod_rast_less.dims[0]) - 1
        #Save file
        meta.update({"dtype": "int16","count": 1,"compress":"lzw"})
        print("Saving")
        with rasterio.open(last_output_file, 'w', **meta) as dst:
            dst.write(last_result,1)
        del last_result
    if metric=='mag' or metric=='sd':
        ######################## MAGNITUDE OF WATER STRESS ##########################
        print("Summed magnitude")
        mid_cent_mod = mod_rast[0:35,:,:]/10
        end_cent_mod = mod_rast[35:,:,:]/10
        #Integrated amount below pivot (intensity)
        slp_rast = np.broadcast_to(slp,mid_cent_mod.shape)
        yint_rast = np.broadcast_to(yint,mid_cent_mod.shape)
        mid_mag = np.multiply(mid_cent_mod,slp_rast)+yint_rast
        end_mag = np.multiply(end_cent_mod,slp_rast)+yint_rast
        del slp_rast,  yint_rast
        if metric=='sd':
            mid_result = np.std(mid_mag, axis = 0)
            end_result = np.std(end_mag, axis = 0)
        if metric=='mag':
            mid_result = mid_mag.sum(mid_mag.dims[0])
            end_result = end_mag.sum(end_mag.dims[0])
        print("Saving")
        meta.update({"dtype": "float64","count": 1,"compress":"lzw"})
        with rasterio.open(mid_output_file, 'w', **meta) as dst:
            dst.descriptions = tuple([metric])
            dst.write(mid_result,1)
        del mid_mag
        with rasterio.open(end_output_file, 'w', **meta) as dst:
            dst.descriptions = tuple([metric])
            dst.write(end_result,1)
        del end_mag
    if metric=='total':
        ######################## FREQ AND DURATION OF WATER STRESS ##########################
        mid_cent_mod = mod_rast_less[0:35,:,:]
        end_cent_mod = mod_rast_less[35:,:,:]
        #Total below pivot point (frequency)
        mid_tot = mid_cent_mod.sum(mid_cent_mod.dims[0])
        #
        end_tot = end_cent_mod.sum(end_cent_mod.dims[0])
        #Longest consecutive stretch below pivot point (duration)
        mid_cmod = mid_cent_mod.cumsum(mid_cent_mod.dims[0])
        mid_cumbelow = mid_cmod-mid_cmod.where(np.logical_not(mid_cent_mod)).ffill(mid_cent_mod.dims[0]).fillna(0)
        mid_cumbelowmax = mid_cumbelow.max(mid_cumbelow.dims[0])
        #
        end_cmod = end_cent_mod.cumsum(end_cent_mod.dims[0])
        end_cumbelow = end_cmod-end_cmod.where(np.logical_not(end_cent_mod)).ffill(end_cent_mod.dims[0]).fillna(0)
        end_cumbelowmax = end_cumbelow.max(end_cumbelow.dims[0])
        #Stack  and  save
        mid_result = np.stack((mid_tot, mid_cumbelowmax))
        del mid_tot, mid_cumbelow, mid_cmod
        #Mask back out the fully-nan cells
        nan_cells = np.broadcast_to(nan_cells,mid_result.shape)
        mid_result = xr.where(np.isnan(nan_cells), np.nan, mid_result)
        #Save file
        print("Saving")
        meta.update({"dtype": "int16","count": 2,"compress":"lzw"})
        with rasterio.open(mid_output_file, 'w', **meta) as dst:
            dst.descriptions = tuple(['total', 'longest'])
            dst.write(mid_result)
        del mid_result
        end_result = np.stack((end_tot, end_cumbelowmax))
        del end_tot, end_cumbelow, end_cmod 
        #Mask back out the fully-nan cells
        #nan_cells = np.broadcast_to(nan_cells,end_result.shape)
        end_result = xr.where(np.isnan(nan_cells), np.nan, end_result)
        #Save file
        print("Saving")
        meta.update({"dtype": "int16","count": 2,"compress":"lzw"})
        with rasterio.open(end_output_file, 'w', **meta) as dst:
            dst.descriptions = tuple(['total', 'longest'])
            dst.write(end_result)
        del end_result
    en = time.time()
    print(en-st)
    return 0


park = "CONUS"
win  = "window_14days/cutoff_2.5/bilinear"
#base_path = '/home/ubuntu/data/'
base_path = '/home/1024/ma/sherman/NPS_SIP/pivot_pts/'
base_path = '/home/scratch/sherman/'

pivs_path = base_path+'results/'+park+'/pivot_pts/'+win+'/tiled/EVI/'
mod_path = base_path+'results/'+park+'/wb_proj/futures/tiled/'
out_path = base_path+'results/'+park+'/future_comparison/tiled/'

wbs = ["accumswe", "aet", "deficit","pet", "runoff","soilwater"]
#all_mod_files = []
#for wb in wbs:
#    print(wb)
#    all_mod_files = all_mod_files+glob.glob(os.path.join(mod_path, '*_{wb}_*.tif'.format(wb = wb)))

#import random
#all_mod_files = random.shuffle(all_mod_files)
#pool = multiprocessing.Pool()
#results = pool.starmap(comp_mod_piv, zip(all_mod_files))
#pool.close()   
#for fils in all_mod_files:
#    comp_mod_piv(fils)

#####################################################################33
## Need  to move data from USB to scratch by tile
park = "CONUS"
win  = "window_14days/cutoff_2.5/bilinear"
#base_path = '/home/ubuntu/data/'
base_path = '/home/1024/ma/sherman/NPS_SIP/pivot_pts/'
base_path = '/home/scratch/sherman/'

pivs_path = base_path+'results/'+park+'/pivot_pts/'+win+'/tiled/EVI/'
mod_path = base_path+'results/'+park+'/wb_proj/futures/tiled/'

#wbs = ["accumswe", "aet", "deficit","pet", "runoff","soilwater"]
xoffs = ['0', '10110','20220']
xoffs = ['0']
yoffs = ['0','2611','5222','7833','10444']
yoffs = ['7833','10444']
#wbs = ["deficit","accumswe"]
wbs = ["deficit","accumswe","aet","soilwater","runoff", "pet"]
#wbs = ["aet","soilwater","runoff", "pet"]
for wb in wbs:
    print(wb)
    for xoff in xoffs:
        print(xoff)
        for yoff in yoffs:
            print(yoff)
            #copy over model projections
            cp_from_usb = 'cp -n /home/16t-usb/sherman2/NPS_SIP/pivot_pts/results/CONUS/wb_proj/futures/tiled/*'+wb+'_'+xoff+'_'+yoff+'.tif /home/scratch/sherman/results/CONUS/wb_proj/futures/tiled'
            print(cp_from_usb)
            os.system(cp_from_usb)
            #copy over pivot points
            #piv_from_usb = 'cp -n /home/16t-usb/sherman2/NPS_SIP/pivot_pts/results/CONUS/pivot_pts/window_14days/cutoff_2.5/bilinear/tiled/EVI/*'+wb+'_'+xoff+'_'+yoff+'.tif /home/scratch/sherman/results/CONUS/pivot_pts/window_14days/cutoff_2.5/bilinear/tiled/EVI'
            #print(piv_from_usb)
            #os.system(piv_from_usb)
            all_mod_files = []
            all_mod_files = all_mod_files+glob.glob(os.path.join(mod_path, '*_{wb}_*.tif'.format(wb = wb)))
            for fil in all_mod_files:
                comp_mod_piv(fil, metric = "sd")
                os.remove(fil)


