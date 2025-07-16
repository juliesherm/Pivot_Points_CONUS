## ---------------------------
##
## Script name: calculate_pivotpts.R
##
## Purpose of script: Calculates pivot points optionally for specified park, default is CONUS
##
## Author: Julie Sherman, Carolyn Lober
##
## Date Created: 2/8/2024
## Last edited: 3/15/2024
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

# load packages
library(terra)
library(tidyverse)
library(here)
library(Kendall)
library(xts)
library(signal)
library(stringr)
library(pracma) # for trapz integration

# load functions
source(here("src", "getlm.R")) 
source(here("src", "wateryear.R"))
source(here("src", "read_in_data.R"))
source(here("src", "stk_trend.R")) 
source(here("src", "get_lm_trend.R")) # used in stk_trend()
source(here("src", "get_mk_trend.R")) # used in stk_trend()
source(here("src", "smooth_ts.R")) 
source(here("src", "get_dates.R")) 
#source(here("src", "read_in_data_local.R")) 
source(here("src", "vi_summary_gs.R"))
source(here("src", "integrate_gs.R")) # used in vi_summary_gs() 

park <- "GYE"
# should work with any park or I&M network standard 4 letter code; for I&M networks,
# the script will use a bounding box around all parks in the network 
# for entire CONUS, use "CONUS" or set to NULL

#########################################
# Read in files, assign variable names
#########################################

# read in MODIS variables
vi_stk <- read_in_data(park, "MODIS")
cloud_stk <- read_in_data(park, "Rely")

vi_gs <- read_in_data(park, "MODIS_gs")
names(vi_gs) <- str_split_i(names(vi_gs), "/", 11)
vi_vars <- names(vi_gs)
vi_years <- seq(2000, len = nlyr(vi_gs[[1]]))
vi_dates <- get_dates("MODIS", vi_stk)

# read in GRIDmet variables
gm_stk <- read_in_data(park, source = "GRIDMET")
gm_vars <- names(gm_stk)
gm_dates <- get_dates("GRIDMET", gm_stk)
gm_years <- seq()

# read in water balance variables
wb_wy <- read_in_data(park, source = "WB")
names(wb_wy) <- str_split_i(names(wb_wy), "_", 9)
wb_vars <- names(wb_wy)
wb_years <- seq(1980, len = nlyr(wb_wy[[1]]))
#wb_dates <- get_dates("WB", wb_stk)

# check if WB projection matches MODIS; if not, print warning message 
# (recommended to project with python/gdal script as it is much faster than terra/R)
if (crs(wb_stk[[1]]) != crs(vi_stk[[1]])) {
  print("Water balance data is not projected, must be projected to match crs of MODIS")
}

#########################################
# Reliability procedure (mask cloudy and burned VI pixels)
#########################################
# TO-DO
# - remove burned areas

# create mask and mask all vi stacks
cloud_mask <- cloud_stk == 3 # create logical matrix with T = cloudy pixels
vi_stk <- lapply(vi_stk, terra::mask, cloud_mask, maskvalues = TRUE)

# mask snowy pixels with base value of 500 (0.05)
snow_mask <- cloud_stk == 2 
vi_stk <- lapply(vi_stk, terra::mask, snow_mask, maskvalues = TRUE, updatevalue = 500)
rm(cloud_mask, snow_mask, cloud_stk)

#########################################
# Fix winter VI values? # lines 805-824
#########################################
# get maximum winter (Jan) value for each pixel
# set winter minimum to max * 0.6 so any values below the min become the min
# fill gaps with linear interpolation 

# get random pixel ts to test with
set.seed(25809)
rand_ndx <- sample(dim(vi_stk[[1]])[1] * dim(vi_stk[[1]])[2], 1)
vec <- as.data.frame(c(vi_stk[[1]]),cells=TRUE)
vec <- vec[which(vec$cell == rand_ndx),] %>% select(!cell) %>%
  t() %>% as.numeric 

# note that this is only applied to EVI for now
#vi_smth <- terra::app(vi_stk[[1]],\(vec) smooth_ts(vec,vi_dates))

overwrite <- TRUE
vi_smth <- sapply(vi_vars, function(var) {
  fname <- here("results", park, "vi_smoothed", paste0(var, "_smth.tif"))
  if (!overwrite & file.exists(fname)) {
    print(paste0("reading in smoothed data for ", var))
    sm <- terra::rast(fname)
  } else { 
    print("smoothed vegetation index data not found, smoothing now")
    print(var)
    sm <- terra::app(vi_stk[[var]], \(vec) smooth_ts(vec,vi_dates))
    writeRaster(sm, fname, overwrite = overwrite)
  }
  return(sm)
})

# remove snowy pixels

#########################################
# Summarize over wateryear and growing season
#########################################
# growing season NDVI vs water year avg

### new: get gs_summary
vi_gs_cs <- sapply(vi_vars[1], function(var) {
  vi_summary_gs(vi_smth[[var]], vi_dates, gs_dates_cs, var)
})

### old:
vi_smth <- lapply(vi_smth, \(stk) stk[[which(year(vi_dates) == 2023)]])
vi_dates <- vi_dates[-550]
gs_dates <- terra::rast(here::here("results","CONUS","gs_dates","gs_all_doys_proj.tif"))

# only use years with wb data
vi_gs <- lapply(vi_gs, function(stk) {
  return(stk[[vi_years %in% wb_years]])
})
vi_gs_anom <- lapply(vi_gs, terra::app, scale, center = TRUE, scale = FALSE)
names(vi_gs_anom) <- names(vi_gs)

# only use wb years with vi data
wb_wy_small <- lapply(wb_wy, function(stk) {
  return(stk[[wb_years %in% vi_years]])
})


# get water year mean or sum
wb_sum_vars <- c("aet", "pet", "deficit", "rain")
wb_mean_vars <- c("soilwater")
wb_max_vars <- c("accumswe")
wb_wy_sum <- lapply(wb_stk[wb_sum_vars], terra::tapp, index = wateryear(wb_dates), 
                  fun = sum, na.rm = TRUE)
wb_wy_mean <- lapply(wb_stk[wb_mean_vars], terra::tapp, index = wateryear(wb_dates), 
                   fun = mean, na.rm = TRUE)
wb_wy_max <- lapply(wb_stk[wb_max_vars], terra::tapp, index = wateryear(wb_dates), 
                  fun = max, na.rm = TRUE)
wb_wy <- c(wb_wy_sum, wb_wy_mean, wb_wy_max)
rm(wb_wy_sum, wb_wy_mean, wb_wy_max)

# data starts in January, so the first water year is incomplete and removed
# data ends in December, so the last water year is incomplete and removed 
wb_wy <- lapply(wb_wy, function(x) subset(x, -c(1,nlyr(x))))
wb_years <- unique(wateryear(wb_dates))[-c(1,length(unique(wateryear(wb_dates))))]

return_na <- function(vec) {
  return(rep(NA, times = length(vec)))
}
# get data going back 2 and 3 years
wb_wy2 <- lapply(wb_wy, function(stk) {
  my_rast <- rast(sapply(wb_years, function(yr) {
    if (yr < 2002) {
      return(terra::app(stk[[which(wb_years == yr)]], return_na))
    } else {
      return(terra::app(stk[[which(wb_years %in% c(yr, yr - 1))]], mean))
    }
  }))
  names(my_rast) <- wb_years
  return(my_rast)
})
names(wb_wy2) <- paste0(names(wb_wy), "2")
wb_wy3 <- lapply(wb_wy, function(stk) {
  my_rast <- rast(sapply(wb_years, function(yr) {
    if (yr < 2003) {
      return(terra::app(stk[[which(wb_years == yr)]], return_na))
    } else{
      return(terra::app(stk[[which(wb_years %in% c(yr, yr - 1, yr - 2))]], mean))
    }
  }))
  names(my_rast) <- wb_years
  return(my_rast)
})
names(wb_wy3) <- paste0(names(wb_wy),"3")

wb_wy <- c(wb_wy,wb_wy2,wb_wy3)
rm(wb_wy2, wb_wy3)

#########################################
# Scale data
#########################################

# scale vi data by 10000 
vi_gs_anom <- lapply(vi_gs_anom, function(x) {
  x / 10000
})

# scale wb data by 10 to get mm
wb_wy <- lapply(wb_wy, function(x) {
  x / 10
}) 


#########################################
# Calculate temporal trend in VI 
#########################################

if (!overwrite & file.exists(here("results", park, "vi_trends"))) {
  print("reading existing files")
  vi_mean_annual_trend <- sapply(
    list.files(here("results", park, "vi_trends"), pattern = "mean", full.names = TRUE),
    terra::rast
  )
  vi_max_annual_trend <- sapply(
    list.files(here("results", park, "vi_trends"), pattern = "max", full.names = TRUE),
    terra::rast
  )
  vi_gs_annual_trend <- sapply(
    list.files(here("results", park, "vi_trends"), pattern = "gs", full.names = TRUE),
    terra::rast
  )
} else {
  print("calculating trends and writing to files")
  # get lm trend and MannKendall trend for vegetation indices
  vi_mean_annual_trend <- lapply(vi_smth, \(stk) stk_trend(stk, vi_dates, "mean"))
  vi_max_annual_trend <- lapply(vi_smth, \(stk) stk_trend(stk, vi_dates, "max"))
  vi_gs_annual_trend <- lapply(vi_gs, \(stk) stk_trend(stk, NULL, "gs"))
  
  # write rasters out to folder
  for (var in names(vi_smth)) {
    writeRaster(vi_mean_annual_trend[[var]], here("results", park, "vi_trends", paste0(var, "_mean_trend.tif")),
                overwrite = overwrite)
    writeRaster(vi_max_annual_trend[[var]], here("results", park, "vi_trends", paste0(var, "_max_trend.tif")), 
                overwrite = overwrite)
    writeRaster(vi_gs_annual_trend[[var]], here("results", park, "vi_trends", paste0(var, "_gs_trend.tif")),
                overwrite = overwrite)
  }
}

# set names of list
names(vi_mean_annual_trend) <- names(vi_smth)
names(vi_max_annual_trend) <- names(vi_smth)
names(vi_gs_annual_trend) <- names(vi_smth)

#########################################
# Calculate pivot points, responses, correlations with climate vars (lines 1467+ in ndvi_3_1)
#########################################

# lm method below inspired by https://gist.github.com/bennyistanto/1a691954c68d1bc11d2fb1c4f2df2b6d

# get all combos of VIs and WB vars
vi_wb_prs <- expand.grid(names(vi_gs), names(wb_wy))

# for each combo of VI and WB var, apply the lm to every pixel
overwrite <- FALSE
start <- Sys.time()
pivotpts <- apply(vi_wb_prs, 1, function(x) {
  # check if pivot points have already been calculated either for the park or CONUS
  if (!overwrite & file.exists(here("results", park, "pivotpts", paste0(x[1], "_", x[2], ".tif")))) {
    print("Pivot points have already been calculated, and overwrite is not TRUE. Reading in existing files")
    my_pivotpts <- terra::rast(here("results", park, "pivotpts", paste0(x[1], "_", x[2], ".tif")))
    
  } else if (!overwrite & file.exists(here("results", "CONUS", "pivotpts", paste0(x[1], "_", x[2], ".tif")))) {
    print("Pivot points have already been calculated for CONUS and overwrite is not TRUE. Reading in and cropping existing files.")
    
    # load CONUS pivot points
    my_pivotpts <- terra::rast(here("results", "CONUS", "pivotpts", paste0(x[1], "_", x[2], ".tif")))
    
    # load park boundary and subset
    park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                 layer = "CONUS_Park_AOAs"
    )
    park_boundary <- subset(park_boundary, park_boundary$UNITCODE == park) %>%
      terra::project(., crs(dat[[1]]))
    my_pivotpts <- lapply(my_pivotpts, function(lyr) {
      terra::crop(lyr, park_boundary)
    })
  } else {
    print(paste0("Calculating pivot points for ", x[1], " and ", x[2]))
    # calculate pivot points from scratch
    my_vi <- vi_gs_anom[[x[1]]]
    my_wb <- wb_wy_small[[x[2]]]
    
    # check that dimensions match
    if (dim(my_vi)[3] != dim(my_wb)[3]) {
      stop("Number of layers must match")
    }
    
    # stack rasters
    my_stk <- c(my_vi, my_wb)
    
    # calculate pivot points for each pixel - add std error for x intercept
    my_pivotpts <- terra::app(my_stk, getlm)
    names(my_pivotpts) <- c("y.intercept", "slope", "r.squared", "x.intercept",
                            "p-value.int", "p-value.slp", "x.min", "x.max")
    
    # when calculating new pivot points, write rasters
    terra::writeRaster(my_pivotpts, here("results", park, "pivotpts", paste0(x[1], "_", x[2], ".tif")), overwrite = overwrite)
  }
  
  # return calculated or read files
  return(my_pivotpts)
})
end <- Sys.time()
print(paste("Calculating pivot points took", round(difftime(end, start, units = "secs")/60, 3), "minutes to finish"))

# name pivot pts with descriptive name
names(pivotpts) <- paste(vi_wb_prs$Var1, vi_wb_prs$Var2, sep = "_")

