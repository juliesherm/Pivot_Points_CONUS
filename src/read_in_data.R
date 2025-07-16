## ---------------------------
##
## Script name: read_in_data.R
##
## Purpose of script: reads in data from either MODIS, GRIDMET, WB (water balance model), DEM,
## or vegmap
##
## Author: Julie Sherman, Carolyn Lober
##
## Date Created:
## Last edited: 4/3/2024
##
## ---------------------------
##
## Notes:
##  - park should be supplied as official four letter abbreviation for the park or 
##    I&M network (e.g. "NCPN") or "CONUS"
##  - data must be stored in the file structure described in the project's README
##  - currently hard-coded that WB always gets monthly
##  - valid sources are "MODIS", "MODIS_gs", "Rely", "GRIDMET", "WB", "vegmap", "pivotpts", 
##    "DEM", and "vi_trends"
##
## ---------------------------

read_in_data <- function(park = "CONUS", source) {
  # print current task
  print(paste0("Getting ", source, " variables for ", ifelse(is.null(park), "CONUS", park)))
  
  # load modis data
  if (source == "MODIS") {

    # identify the VI variables available
    vi_vars <- list.files(path = here("raw_data", "MODIS", "CONUS", "VI_16Days_250m_v61"))[str_detect(list.files(path = here("raw_data", "MODIS", "CONUS", "VI_16Days_250m_v61")), "VI")]
    if (length(vi_vars) == 0) {
      stop("Failed to find MODIS vegetation index data")
    }
    
    # load CONUS raster stack for each VI variable
    dat <- sapply(vi_vars, function(var) {
      terra::rast(list.files(
        path = here("raw_data", "MODIS", "CONUS", "VI_16Days_250m_v61", var, "tiled"),
        pattern = "*_0_0.tif$",
        full.names = TRUE
      ))
    })
  } else if (source == "MODIS_gs") {
    
    if (dir.exists(here("results", "CONUS", "integrated_vi_gs"))) {
      vi_vars <- list.dirs(path = here("results", "CONUS", "integrated_vi_gs", "window_14days", "cutoff_2.5", "bilinear"), 
                           full.names = FALSE, 
                           recursive = FALSE)
      
      dat <- sapply(vi_vars[1], function(var) {
        fnames <- list.files(
          path = here("results", "CONUS", "integrated_vi_gs", "window_14days", "cutoff_2.5", "bilinear", var, "tiled"), 
          pattern = "*_0_0.tif$",
          full.names = TRUE
        )
        sapply(fnames, terra::rast)
      })
    }
    
    
    # load cloud mask layer (also from MODIS)
  } else if (source == "Rely") {
    dat <- terra::rast(list.files(
      path = here("raw_data", "MODIS", "CONUS", "VI_16Days_250m_v61", "Rely", "tiled"),
      pattern = "*_0_0.tif$",
      full.names = TRUE
    ))
    
    # load GRIDmet data
  } else if (source == "GRIDMET") {

    # if projected data (matching MODIS, EPSG 4326) exists, use it instead
    if (file.exists(here("results", "CONUS", "gm_proj"))) {
      print("Reading projected wb data")
      
      # get gm vars
      gm_vars <- list.dirs(here("results", "CONUS", "gm_proj"),
                           full.names = FALSE, recursive = FALSE
      )
      gm_vars <- gm_vars[-which(gm_vars == "tiled")]
      
      # load raster stack for each water balance variable
      dat <- sapply(gm_vars, function(var) {
        terra::rast(list.files(
          path = here("results", "CONUS", "gm_proj", var),
          pattern = "tif$",
          full.names = TRUE
        ))
      })
    } else {
      
      # identify wb vars available
      gm_vars <- list.dirs(here("raw_data", "GRIDmet", "wateryear"),
                           full.names = FALSE, recursive = FALSE
      )
      
      # load raster stack for each water balance variable
      dat <- sapply(gm_vars, function(var) {
        terra::rast(list.files(
          path = here("raw_data", "GRIDmet", "wateryear", var),
          pattern = "nc4$",
          full.names = TRUE
        ))
      })
    }
    
    # read water balance data
  } else if (source == "WB") {
    
    # if projected data (matching MODIS, EPSG 4326) exists, use it instead
    if (file.exists(here("results", "CONUS", "wb_proj", "wateryear"))) {
      print("Reading projected wb data (nn)")
      
      # get wb vars
      files <- list.files(here("results", "CONUS", "wb_proj", "wateryear", "nearest_neighbor"),
                           pattern = "tif$",
                           full.names = TRUE, recursive = FALSE
      )
      
      # load raster stack for each water balance variable
      dat <- sapply(files, terra::rast)
    } else {
      
      # identify wb vars available
      wb_vars <- list.dirs(here("raw_data", "WB", "wateryear"),
                           full.names = FALSE, recursive = FALSE
      )
      
      # load raster stack for each water balance variable
      dat <- sapply(wb_vars, function(var) {
        terra::rast(list.files(
          path = here("raw_data", "WB", "wateryear", var),
          pattern = "nc4$",
          full.names = TRUE
        ))
      })
    }
    
    # load vegmap data
  } else if (source == "vegmap") {
    # looks for a shp file in the vegmaps/park folder
    fname <- list.files(
      path = here("raw_data", "vegmaps", park),
      pattern = "shp$",
      full.names = TRUE
    )
    
    if (file.exists(fname)) {
      dat <- terra::vect(fname)
    } else {
      stop("Can't find vegmap for specified park")
    }
  } else if (source == "pivotpts") {
    # check if pivotpts are saved for specified park (or CONUS)
    if (dir.exists(here("results", park, "pivotpts"))) {
      # read in data
      fname <- list.files(
        path = here("results", park, "pivotpts"),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if pivot points exist for the park, return early without cropping step below
      dat <- sapply(fname, terra::rast)
      
      # get descriptive name from file name without extension
      names(dat) <- str_split_i(list.files(
        path = here("results", park, "pivotpts"),
        pattern = "tif$",
        full.names = FALSE
      ), "\\.", 1)
      
      return(dat)
      
    } else if (park != "CONUS" & dir.exists(here("results", "CONUS", "pivotpts"))) {
      # if pivot points don't exist for the park but do for CONUS, read from CONUS
      fname <- list.files(
        path = here("results", "CONUS", "pivotpts"),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if pivot points are being read from CONUS for a smaller specified park, don't return early
      # so that they are cropped in the next step
      dat <- sapply(fname, terra::rast)
    } else {
      stop("Can't find pivot points.")
    }
    
  } else if (source == "DEM") {
    
    # check if pivotpts are saved for specified park (or CONUS)
    if (dir.exists(here("raw_data", "DEM", park))) {
      # read in data
      fname <- list.files(
        path = here("raw_data", "DEM", park),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if pivot points exist for the park, return early without cropping step below
      dat <- sapply(fname, terra::rast)
      
      return(dat)
      
    } else if (park != "CONUS" & dir.exists(here("raw_data", "DEM", "CONUS"))) {
      # if pivot points don't exist for the park but do for CONUS, read from CONUS
      fname <- list.files(
        path = here("raw_data", "DEM", "CONUS"),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if pivot points are being read from CONUS for a smaller specified park, don't return early
      # so that they are cropped in the next step
      dat <- sapply(fname, terra::rast)
    } else {
      stop("Can't find dem.")
    }
    
  } else if (source == "vi_trends") {
    
    # check if vi_trends are saved for specified park (or CONUS)
    if (dir.exists(here("results", park, "vi_trends"))) {
      # read in data
      fname <- list.files(
        path = here("results", park, "vi_trends"),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if vi_trends exist for the park, return early without cropping step below
      dat <- sapply(fname, terra::rast)
      
      return(dat)
      
    } else if (park != "CONUS" & dir.exists(here("results", "CONUS", "vi_trends"))) {
      # if vi_trends don't exist for the park but do for CONUS, read from CONUS
      fname <- list.files(
        path = here("results", "CONUS", "vi_trends"),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if vi_trends are being read from CONUS for a smaller specified park, don't return early
      # so that they are cropped in the next step
      dat <- sapply(fname, terra::rast)
    } else {
      stop("Can't find vi_trends")
    }
    
  } else if (source == "vi_gs") {
    
    # check if vi_gs are saved for specified park (or CONUS)
    if (dir.exists(here("results", park, "integrated_vi_gs"))) {
      # read in data
      fname <- list.files(
        path = here("results", park, "integrated_vi_gs"),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if vi_gs exist for the park, return early without cropping step below
      dat <- sapply(fname, terra::rast)
      
      return(dat)
      
    } else if (park != "CONUS" & dir.exists(here("results", "CONUS", "integrated_vi_gs"))) {
      # if vi_gs don't exist for the park but do for CONUS, read from CONUS
      fname <- list.files(
        path = here("results", "CONUS", "integrated_vi_gs"),
        pattern = "tif$",
        full.names = TRUE
      )
      
      # if vi_gs are being read from CONUS for a smaller specified park, don't return early
      # so that they are cropped in the next step
      dat <- sapply(fname, terra::rast)
    } else {
      stop("Can't find existing folder integrated_vi_gs") 
    }
    
  } else {
    stop("Source not defined, can't read data")
  }
  
  # if park or I&M network is supplied, subset CONUS data to only the park
  if (park == "GYE" & source != "vegmap") {
    gye_boundary <- terra::vect(here("raw_data", "GYE_boundary", "GYE_boundary_dd", "GYE_boundary_dd.shp"))


  } else if (park != "CONUS" & source != "vegmap") {
    network_codes <- unique(terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                        layer = "NPS_IMD_EnvironmentalSettingsParks")$NetworkCode)
    
    if (park == "GYE") {
      park_boundary <- terra::vect(here("raw_data", "GYE_boundary", "GYE_boundary_dd", "GYE_boundary_dd.shp"))
    } else if (park %in% network_codes) {
      park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                   layer = "CONUS_Park_AOAs"
      )
      park_boundary <- subset(park_boundary,park_boundary$NetworkCode == park) %>%
        terra::project(., crs(dat[[1]]))
    } else {
      park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                   layer = "CONUS_Park_AOAs"
      )
      park_boundary <- subset(park_boundary, park_boundary$UNITCODE == park) %>%
        terra::project(., ifelse(is.na(crs(dat[[1]])), crs(dat[[1]][[1]]), crs(dat[[1]])))
    }
    
    if(is.list(dat)){
      print("Cropping to park extent")
      dat <- lapply(dat, function(lyr) {
        terra::crop(lyr, terra::ext(park_boundary))
      })
      
    }else{
      # for single variable, avoid lapply for desired output 
      print("Cropping to park extent")
      dat <- terra::crop(dat, terra::ext(park_boundary))
    }
  }
  
  # return specified data
  return(dat)
}
