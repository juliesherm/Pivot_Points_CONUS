# Pivot_Points_CONUS
Scripts used to create and analyze pivot points for the continental United States. 

Data and methods for this work is heavily inspired by  
 -  Landscape pivot points and responses to water balance in national parks of the southwest US by Thoma, Munson, and Witwicki [https://doi.org/10.1111/1365-2664.13250](https://doi.org/10.1111/1365-2664.13250)
 -  Robust projections and consequences of an expanding bimodal growing season in the western United States by Tercek, Gross, and Thoma [https://doi.org/10.1002/ecs2.4530](https://doi.org/10.1002/ecs2.4530)


The workflow for this project is approximately as follows:

## Download and process MODIS data

- Download MODIS VI and Rely data files
- [modis_yr_stack_tile.sh](src/modis_yr_stack_tile.sh) converts raw data to multibanded tiff 'tiles'. Each tiff file covers 1 year and 1/15 spatial tile over CONUS. Each band is one of the 16-day composite values.


## Download and process water balance data

- gridmet_processor.py resamples waterbalance data to match MODIS resolution
- 

## Download and process temperature and precipitation data

- gridmet_processor.py resamples gridMET data to match to MODIS resolution/CRS
- [gm_wateryear_average.py](src/gm_wateryear_average.py) calculates the wateryear (October-September) average for gridmet data
- [gm_yr_stack_tile.sh](src/gm_yr_stack_tile.sh) converts projected wateryear average to multibanded tiff 'tiles'. Each tiff file covers all years and 1/15 spatial tile over CONUS. Each band is one year.



## Define growing season

- [start_end_gs_CONUS.py](src/start_end_gs_CONUS.py) calculates the start and end date for each year based on AET values.
- [combine_gs_mats.R](src/combine_gs_mats.R) computes the median over the years
- [gs_tile_tiffs.sh](src/gs_tile_tiffs.sh) cuts up the growing season tiffs into 15 tiles 




## Calculate cummulative annual vegetation production

- [vi_smooth_tiles.py](src/vi_smooth_tiles.py) smooths, gap-fills and integrates over the growing season for each year/tile to obtain annual vegetation summary.
- [int_vi_yr_stack.sh](src/int_vi_yr_stack.sh) combines yearly integrated VI into single file. Each file is one tile with bands corresponding to years (2001--2023)





## Calculate pivot points

- [pivot_tiles.py](src/pivot_tiles.py) calculates pivot points from tiled yearly waterbalance and iEVI data (which is also centered during this script)


## Predict future VI using pivot points and water balance projections

Note that future, yearly VI projections are never explictly calculated/saved for reasons of storage space. Instead, future values of water balance variables are compared to the pivot points directly to produced summary values. Using the `magnitude' option in the compare_future file allows direct computation of VI projections, but it is subsequently summed over the mid- or end-century periods

- [compare_future.py](src/compare_future.py) creates metrics to summarize future VI projections in relation to the pivot points   


## Analysis of results

### Assess results by land classification
- [project_nlcd.R](src/project_nlcd.R) projects NLCD data to match MODIS CRS/resolution. Also contains plotting functionality for the NLCD data 
- [sig_table.R](src/sig_table.R) creates a table of the proportion of significant pixels by independent WB predictor and land class (from NLCD) 



### Plotting

- [plot_gs.R](src/plot_gs.R) plots the growing season for CONUS or for BLCA

## Misc Functions
- [avg_tmin_tmax.sh](src/avg_tmin_tmax.sh) calculates the average between tmin and tmax daily temps (although I believe the T_avg variable was ultimately dropped from analysis)
- [range_tiles.py](src/range_tiles.py) calculates the range of the given waterbalance data
- There are several shell functions that were used to tile/untile various variables. For instance:
  - [pivs_CONUS.sh](src/pivs_CONUS.sh) combines pivot point output tiles for a given VI/WB pair and output variable of interest (intercept, slope, rsquared, etc)

