# Pivot_Points_CONUS
Scripts used to create and analyze pivot points for the continental United States. 

The workflow for this project is approximately as follows:

## 1. Download and process MODIS data

A. Download MODIS VI and Rely data files

B. [modis_yr_stack_tile.sh](src/modis_yr_stack_tile.sh) converts raw data to multibanded tiff 'tiles'. Each tiff file covers 1 year and 1/15 spatial tile over CONUS. Each band is one of the 16-day composite values.


## 2. Download and process water balance data

A. gridmet_processor.py resamples waterbalance data to match MODIS resolution

## 3. Define growing season



## 4. Calculate cummulative annual vegetation production

A. [vi_smooth_tiles.py](src/vi_smooth_tiles.py) smooths, gap-fills and integrates over the growing season for each year/tile to obtain annual vegetation summary.


## 3. Calculate pivot points

A. [pivot_tiles.py](src/pivot_tiles.py) calculates pivot points from tiled yearly waterbalance and iEVI data (which is also centered during this script)


## 4. Predict future VI using pivot points and water balance projections



## 5. Create metrics to summarize future VI projections

A. [compare_future.py](src/compare_future.py)

## 6. Plotting
