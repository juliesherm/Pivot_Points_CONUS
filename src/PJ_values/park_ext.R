park_ext <- function(park, out_format = "extent") {
  #get the park extent in the crs of MODIS
  raw_vi<-lapply(list.files(here("raw_data","MODIS","CONUS","VI_16Days_250m_v61","EVI","tiled"),
                            pattern = paste0("EVI","_2000_*"), full.names = T),
                 terra::rast)
  
  network_codes <- unique(terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                      layer = "NPS_IMD_EnvironmentalSettingsParks")$NetworkCode)
  
  if (park %in% network_codes) {
    park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                 layer = "CONUS_Park_AOAs"
    )
    park_boundary <- subset(park_boundary,park_boundary$NetworkCode == park) %>%
      terra::project(., crs(raw_vi[[1]]))
  } else{
    park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                 layer = "CONUS_Park_AOAs"
    )
    park_boundary <- subset(park_boundary, park_boundary$UNITCODE == park) %>%
      terra::project(., crs(raw_vi[[1]]))
  }
  if(out_format =="extent"){
    park_extent<-terra::ext(park_boundary)
    return(park_extent)  
  }else if(out_format=="boundary"){
    return(park_boundary)
  }
}
  