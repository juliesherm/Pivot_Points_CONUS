library(terra)
library(here)
library(stringr)

## this script combines the csvs's output from start_end_CONUS.py into .tif raster files
## probably should just fix how they're saved in that script...
## also calculates the median per-pixel growing season across the years
## saves as yearly files and all together
#wb_proj_template<-terra::rast(here("raw_data","WB","daily","aet","v_1_5_2021_gridmet_historical_aet.nc4"))
wb_proj_template<-terra::rast(here("raw_data","WB","monthly","aet","v_1_5_2000_gridmet_historical_aet_monthly.nc4"))
#loop over the years first and last and save as tif with same projection info as wb_proj_template
#(slow, but should only do this once)
st<-Sys.time()
win = "window_14days"
wb_ext<-terra::ext(wb_proj_template)
wb_crs<-terra::crs(wb_proj_template)
for(yr in 2022:2023){
  print(yr)
  fname<-here("results","CONUS","gs_dates",win,"cutoff_2.5",paste0(c("first","last"),"_inds_",yr,".csv"))
  gs<-c(rast(as.matrix(read.table(fname[1],sep = ",", header = FALSE)),
                 crs = wb_crs,  extent = wb_ext),
        rast(as.matrix(read.table(fname[2], sep = ",", header = FALSE)),
                crs = wb_crs,  extent = wb_ext))
  names(gs)<-paste0(c("first","last"),yr)
  terra::writeRaster(gs, here("results","CONUS","gs_dates",win,"cutoff_2.5",paste0("gs_",yr,".tif")), overwrite = TRUE)
  
}
Sys.time()-st #12  minutes, not bad


#read in all years and combine as one tif file with layer medians
gs_tifs<-list.files(path = here("results","CONUS","gs_dates",win,"cutoff_2.5"),
                    pattern = ".tif$",
                    full.names = TRUE)
gs<-terra::rast(gs_tifs)
#calculate median first and last growing season date
st<-Sys.time()
gs<-c(gs,tapp(gs, index = str_sub(names(gs), end = -5), median))
names(gs)[nlyr(gs)]<-"medianlast"
names(gs)[nlyr(gs)-1]<-"medianfirst"
terra::writeRaster(gs, here("results","CONUS","gs_dates",win,"cutoff_2.5","gs_all_doys.tif"), overwrite = TRUE)
Sys.time()-st # 1 minute!

#gs<-terra::rast(here::here("results","CONUS","gs_dates",win,"gs_all_doys.tif"))

#project to MODIS resolution 
#(or use gdal - so much faster)
#vi_proj_res<-terra::rast(here("raw_data","MODIS","CONUS","VI_16Days_250m_v61","EVI","MOD13Q1_EVI_2022_305.tif"))

#st<-Sys.time()
#terra::project(gs, vi_proj_res,method = "near",  
#               filename = here("results","CONUS","gs_dates",win,"gs_all_doys_proj.tif"),
#               overwrite = TRUE)
#Sys.time()-st


