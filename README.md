# Pivot_Points_CONUS
Scripts used to create and analyze pivot points for the continental United States. 

The workflow for this project is approximately as follows:

## Download and process MODIS data

1. Download MODIS VI and Rely data files

2. [modis_yr_stack_tile.sh](src/modis_yr_stack_tile.sh) converts raw data to multibanded tiff 'tiles'. Each tiff file covers 1 year and 1/15 spatial tile over CONUS. Each band is one of the 16-day composite values.


## Download and process water balance data

1. gridmet_processor.py resamples waterbalance data to match MODIS resolution

## Download and process temperature and precipitation data

1. gridmet_processor.py resamples waterbalance data to match MODIS resolution



## Define growing season

1. [start_end_gs_CONUS.py](src/start_end_gs_CONUS.py) calculates the start and end date for each year based on AET values.

2. 



## Calculate cummulative annual vegetation production

1. [vi_smooth_tiles.py](src/vi_smooth_tiles.py) smooths, gap-fills and integrates over the growing season for each year/tile to obtain annual vegetation summary.




## Calculate pivot points

1. [pivot_tiles.py](src/pivot_tiles.py) calculates pivot points from tiled yearly waterbalance and iEVI data (which is also centered during this script)


## Predict future VI using pivot points and water balance projections

Note that future, yearly VI projections are never explictly calculated/saved for reasons of storage space. Instead, future values of water balance variables are compared to the pivot points directly to produced summary values as described in the next section. Using the `magnitude' option in the compare_future file allows direct computation of VI projections, but it is subsequently summed over the mid- or end-century periods


## Create metrics to summarize future VI projections

1. [compare_future.py](src/compare_future.py)



## Plotting

1. [plot_gs.R](src/plot_gs.R) plots the growing season for CONUS or for BLCA

## Misc Functions
