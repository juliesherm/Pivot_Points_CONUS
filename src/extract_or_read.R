## ---------------------------
##
## Script name: extract_or_read.R
##
## Purpose of script: Looks for existing extraction and reads in data, or extracts from dat_raster provided
##
## Author: Carolyn Lober
##
##
## ---------------------------
##
## Notes:
##  - vegmap must have column "Hectares" (could change to calculate the area if this is a problem)
##
## ---------------------------

extract_or_read <- function(park, vegmap, dat_type, dat_raster, overwrite) {
  
  # set filename that extracted data is expected to be saved with
  if (dat_type == "dem") {
    fname = "poly_dem.csv"
  } else if (dat_type == "pivotpts") {
    fname = "poly_pivotpts.csv"
    
    # put into one spatraster stack (faster)
    dat_raster <- rast(dat_raster)
    names(dat_raster) <- basename(sources(dat_raster,bands=TRUE)$source) %>%
      sub(pattern = "(.*)\\..*$", replacement = "\\1", .) %>% 
      paste(.,sources(dat_raster,bands=TRUE)$bands,sep="_")
    
  } else if (dat_type == "vi_trends") {
    fname = "poly_vi_trends.csv"
    
    # put into one spatraster stack
    dat_raster <- rast(dat_raster)
    my_names <- basename(sources(dat_raster, bands = TRUE)$source)
    my_names <- paste(str_split_i(my_names,"_",1),str_split_i(my_names,"_",2),sep="_")
    my_names <- paste(my_names, sources(dat_raster,bands=TRUE)$bands,sep="_")
    names(dat_raster) <- my_names
  } else if (dat_type == "burned_areas") {
    fname = "poly_burned.csv"
  } else {
    stop("dat_type not recognized.")
  }
  
  # either read or extract data
  if (!overwrite & file.exists(here("results", park, fname))) {
    print(paste0("Data already extracted for ", park, ", reading from file results/", park, "/", fname))
    dat <- read.csv(here("results", park, fname))
  } else { 
    
    # extract polygons 
    start <- Sys.time()
    if (dat_type == "pivotpts") {
      dat_shp <- terra::extract(dat_raster,
                                vegmap_8.3,
                                fun = table, 
                                weights = TRUE,
                                bind = TRUE 
      ) #function(x, w) {
        #frac_nona <- sum(!is.na(x) * w, na.rm = TRUE) / sum(w)
        #if (frac_nona >= 0.5) {
        #  mean(x, na.rm = TRUE)
        #}
        #else {
        #  NA
        #}
      #}
    } else {
      dat_shp <- terra::extract(dat_raster,
                                vegmap_8.3,
                                fun = mean, 
                                weights = TRUE,
                                bind = TRUE 
      )
    }
     
    end <- Sys.time()
    print(paste("Extract took", round(difftime(end, start, units = "secs"), 3), "seconds to finish"))
    
    # turn shapefile into dataframe for writing to csv depending on data type
    if (dat_type == "dem") {
      dat <- dat_shp %>% as.data.frame %>%
        select(OBJECTID, elevation, slope, aspect)
    } else if (dat_type == "pivotpts") {
      dat <- dat_shp %>% as.data.frame %>%
        pivot_longer(cols = (ncol(vegmap)+1):ncol(dat_shp), names_to = "wb_var") %>%
        separate_wider_delim(cols = wb_var, delim = "_", names = c("vi_var", "wb_var", "type")) %>%
        mutate(type = case_match(type, 
                                 "1" ~ "y.intercept", 
                                 "2" ~ "slope", 
                                 "3" ~ "r.squared", 
                                 "4" ~ "x.intercept",
                                 "5" ~ "p-value.int",
                                 "6" ~ "p-value.slp",
                                 "7" ~ "x.min",
                                 "8" ~ "x.max")) 
    } else if (dat_type == "vi_trends") {
      dat <- dat_shp %>% as.data.frame %>%
        pivot_longer(cols = (ncol(vegmap)+1):ncol(dat_shp), names_to = "vi_var") %>%
        separate_wider_delim(cols = vi_var, delim = "_", names = c("vi_var", "max_or_mean", "type")) %>%
        mutate(type = case_match(type, 
                                 "1" ~ "y.intercept", 
                                 "2" ~ "slope", 
                                 "3" ~ "r.squared", 
                                 "4" ~ "tau",
                                 "5" ~ "p-value")) 
    } else if (dat_type == "burned_areas") {
      dat <- dat_shp %>% as.data.frame %>%
        select(OBJECTID, USGS_Wildfires_Most_Recent_Year_Burned_Raster)
    }
    
    # write to csv if data was just extracted
    write.csv(dat, here("results", park, fname), row.names = FALSE)
    
  }
  
  return(dat)
}
