## ---------------------------
##
## Script name: read_in_data.R
##
## Purpose of script: reads in data from either MODIS, GRIDMET, WB (water balance model),
## or vegmap
##
## Author: 
##
## Date Created: 
##
## Copyright (c) 
## Email: 
##
## ---------------------------
##
## Notes: 
##  - park should be supplied as official four letter abbreviation for the park
##  - data must be in a uniform file set up as described below:
##  - currently hard-coded that WB always gets monthly
##  - there must be only one .shp file in the relevant vegmap folder
## 
## ---------------------------

read_in_data <- function(park = NULL, source = "MODIS"){
  
  print(paste("Getting ", source, " variables for ",park))
  
  if(source =="MODIS"){
    
    dirname <- list.files(path=here("raw_data","MODIS",park),pattern="VI") 
    
    vi_vars <- list.files(path=here("raw_data","MODIS",park,dirname))[str_detect(list.files(path=here("raw_data","MODIS",park,dirname)),"VI")]
    dat <- sapply(vi_vars,function(var){1
      terra::rast(list.files(path=here("raw_data","MODIS",park,dirname,var), 
                             pattern = "tif$",
                             full.names = TRUE))
    })
    
  } else if(source=="Rely"){
    
    dirname <- list.files(path=here("raw_data","MODIS",park),pattern="VI") 
    
    dat<-terra::rast(list.files(path=here("raw_data","MODIS",park,dirname,"Rely"), 
                                pattern = "tif$",
                                full.names = TRUE))
    
  } else if(source =="GRIDMET"){
    
    GM_vars<-list.dirs(here("raw_data","GRIDmet"), 
                       full.names = FALSE, recursive = FALSE)
    
    dat <- sapply(GM_vars,function(var){
      terra::rast(list.files(path = here("raw_data","GRIDmet",var),
                             pattern = "nc$",
                             full.names = TRUE))
    })
    
  } else if(source =="WB"){
    
    WB_vars<-list.dirs(here("raw_data","WB","monthly"), 
                       full.names = FALSE, recursive = FALSE)
    # read original data
    dat<-sapply(WB_vars,function(var){
      terra::rast(list.files(path = here("raw_data","WB","monthly",var),
                             pattern = "nc4$",
                             full.names = TRUE))
    })
    
  } else if(source =="vegmap"){
    
    fname <- list.files(path = here("raw_data","vegmaps",park),
                        pattern="shp$|gdb$",
                        full.names=TRUE)
    
    if (file.exists(fname)) {
      
      if (str_detect(fname, "gdb$")) {
        dat <- terra::vect(fname,layer = "ROMO_VegPolys")
      } else {
        dat <- terra::vect(fname)
      }
      
      
    } else {
      
      stop("Can't find vegmap for specified park")
      
    }
    
  }
  return(dat)
}
