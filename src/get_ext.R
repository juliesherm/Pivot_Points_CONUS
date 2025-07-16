get_ext<-function(interest_ext, boundary = FALSE){
  if(is.character(interest_ext)){ 
    #if you pass a park name, return extent of that park
    park = interest_ext
    if(interest_ext=="CONUS"){
      return(NULL)
    }
    #get the park boundary and extent in the crs of MODIS
    ref_crs<-lapply(list.files(here("raw_data","MODIS","CONUS","VI_16Days_250m_v61","EVI","tiled"),
                              pattern = paste0("EVI","_2000_*"), full.names = T),
                    terra::rast)
    if(length(ref_crs)==0){
      ref_crs <- terra::rast(list.files(here("results","CONUS","wb_proj","wateryear","bilinear","tiled"),
                             pattern = paste0("*_resampled_*"), full.names = T)[1])
    }else{
      ref_crs = ref_crs[[1]]
    }
    network_codes <- unique(terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                        layer = "NPS_IMD_EnvironmentalSettingsParks")$NetworkCode)
    
    if (park %in% network_codes) {
      park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                   layer = "CONUS_Park_AOAs")
      park_boundary <- subset(park_boundary,park_boundary$NetworkCode == park) %>%
        terra::project(., crs(ref_crs))
    } else{
      park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                                   layer = "CONUS_Park_AOAs")
      park_boundary <- subset(park_boundary, park_boundary$UNITCODE == park) %>%
        terra::project(., crs(ref_crs))
    }
    if(boundary){
      return(park_boundary)
    }else{
      park_extent<-terra::ext(buffer(park_boundary,500))
      # the buffer was 500 (m),  for the previous results, can make larger
      # to show more PJ (~2000)
      return(park_extent)
    }
  }
  else{ 
    #if you want to return the tile of the given park extent
    xoffs = seq(0,30329,by = 10110)
    yoffs = seq(0,13054,by = 2611)
    int_list<-list()
    for(xoff in xoffs){
      for(yoff in yoffs){
        tile_ext<-terra::ext(terra::rast(list.files(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear","tiled"),
                                                    pattern = paste0(xoff,"_",yoff,".tif"),full.names = TRUE)[1]))
        #print(tile_ext)
        intersect<-tryCatch(!is.null(terra::intersect(interest_ext,tile_ext)), error=function(e) return(FALSE)) 
        if(intersect){
          int_list[[length(int_list)+1]]<-c(xoff, yoff)
          #print(c(xoff, yoff))
        }
      }
    }
    return(int_list)
  }
}