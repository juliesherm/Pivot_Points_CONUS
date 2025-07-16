## ---------------------------
##
## Script name: piv_future.R
##
## Purpose of script: Project pivot points and other analysis for specified park
##
## Author: Julie Sherman
##
## Date Created: 2/13/2024
## Last edited: 2/13/2024
##
## ---------------------------
##
## Notes:
##  - 
##
## ---------------------------

#################################
## load packages and functions
#################################
library(terra)
library(here)
library(stringr)
source(here("src","best_lm.R"))

park = "GLAC"

#################################
## Read in WB future projections
#################################
full_names<-list.files(path = here("raw_data","WB","futures","annualnetcdfs"),
                        pattern = "nc4$",
                        full.names = FALSE, recursive = FALSE)
models<-stringr::str_split_i(full_names, "_",5)
rcps<-stringr::str_split_i(full_names, "_",6)
future_vars<-stringr::str_sub(stringr::str_split_i(full_names, "_",7),end = -5)

wb_future<- sapply(list.files(path = here("raw_data","WB","futures","annualnetcdfs"),
                            pattern = "nc4$",
                            full.names = TRUE), terra::rast)
names(wb_future)<-paste(models,rcps,future_vars, sep = "_")

#################################
## Read in pivot point relationships
#################################
vi_vars <- list.files(path = here("raw_data", "MODIS", "CONUS", "VI_16Days_250m_v61"))[str_detect(list.files(path = here("raw_data", "MODIS", "CONUS", "VI_16Days_250m_v61")), "VI")]
wb_vars <- list.dirs(here("raw_data", "WB", "monthly"),
                     full.names = FALSE, recursive = FALSE)
vi_wb_prs <- expand.grid(vi_vars, wb_vars)
pivotpts <- apply(vi_wb_prs, 1, function(x) {
    my_pivotpts <- terra::rast(here("results", park, "pivotpts", paste0(x[1], "_", x[2], "_wbsum.tif")))
  return(my_pivotpts)
})
names(pivotpts)<-apply(vw_wb_prs,1,function(x) paste0(x[1], "_", x[2]))
best_mods<-best_lm(pivotpts)  

#################################
## Project and scale water balance futures
#################################
# scale data by 10 to get mm
sapply(names(wb_future), function(x) {
    writeRaster(terra::project(wb_future[[x]], best_mods)/10,
                here("results",park,"wb_proj","future_proj",paste0(x,".tif")))
})

#################################
## Predict VI at every pixel according to the best relationship
#################################
m<-models[1]
rcp<-rcps[1]

for(m  in unique(models)[5:6]){
  print(m)
  for(rcp in unique(rcps)){
    print(rcp)
    
    if(file.exists(here("results",park,"vi_futures",paste0(m, rcp,".tif")))){ next }
    
    wb_future<-lapply(unique(future_vars), function(x) {
      rast(here("results",park,"wb_proj","future_proj",
                paste0(m,"_", rcp,"_", x,".tif")))})
    names(wb_future)<-paste0(m,"_", rcp,"_", tolower(unique(future_vars)))
    #create blank raster
    vi_future<-wb_future[[1]]  
    values(vi_future)<-NA
    
    for(i in unique(best_mods$which.max)$which.max){
      print(i)
      #name of the water balance future to use
      nm<-paste(m, rcp,stringr::str_split_i(names(pivotpts)[i], "_",2),sep = "_")
      if(is.null(wb_future[[nm]])){ next }
      for(j in 1:nlyr(wb_future[[nm]])){
        #    wb_interest[[j]][best_mods$which.max == i]<-wb_future[[nm]][[j]][best_mods$which.max == i]
        vi_future[[j]][best_mods$which.max == i]<- best_mods$slope[best_mods$which.max == i]*wb_future[[nm]][[j]][best_mods$which.max == i]+
          best_mods$y.intercept[best_mods$which.max == i]                  
        
      }
    }
    writeRaster(vi_future, 
                here("results",park,"vi_futures",paste0(m, rcp,".tif")))
    
  }
}


#####################
## Parallel Version
#####################
# currently doesn't work because variables are not stored in memory
clust<-parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl = clust)

foreach(m =unique(models), .packages = c("here","terra"), 
        .export = c("park","rcps","future_vars","best_mods","pivotpts")) %dopar% {
  
  for(rcp in unique(rcps)){
    print(rcp)
    if(file.exists(here("results",park,"vi_futures",paste0(m, rcp,".tif")))){ next }
      
      wb_future<-lapply(unique(future_vars), function(x) {
        rast(here("results",park,"wb_proj","future_proj",
                  paste0(m,"_", rcp,"_", x,".tif")))})
      names(wb_future)<-paste0(m,"_", rcp,"_", tolower(unique(future_vars)))
      #create blank raster
      vi_future<-wb_future[[1]]  
      values(vi_future)<-NA
      
      for(i in unique(best_mods$which.max)$which.max){
        print(i)
        #name of the water balance future to use
        nm<-paste(m, rcp,stringr::str_split_i(names(pivotpts)[i], "_",2),sep = "_")
        if(is.null(wb_future[[nm]])){ next }
        for(j in 1:nlyr(wb_future[[nm]])){
          #    wb_interest[[j]][best_mods$which.max == i]<-wb_future[[nm]][[j]][best_mods$which.max == i]
          vi_future[[j]][best_mods$which.max == i]<- best_mods$slope[best_mods$which.max == i]*wb_future[[nm]][[j]][best_mods$which.max == i]+
            best_mods$y.intercept[best_mods$which.max == i]                  
          
        }
      }
      writeRaster(vi_future, 
                  here("results",park,"vi_futures",paste0(m, rcp,".tif")))
      
    }
  }
parallel::stopCluster(cl= clust)