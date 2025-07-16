##plot GLAC pivot points
park = "GLAC" # glac is in tile 0_0 only - may need to change below for other parks

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

park_extent<-terra::ext(park_boundary)

vi = "EVI"
wb = "aet"

pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear","tiled",
                       paste0("pivs_",vi,"_",wb,"_0_0.tif")),
                  win = glac_extent)


pdf(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
         paste0(vi,wb,"_pivot_plots.pdf")))
par(mfrow = c(3,2))
plot(pivs,3, main = paste("x-intercept",vi,wb))
plot(pivs,2, main = paste("slope",vi,wb))
plot(pivs,4, main = paste("R2",vi,wb))
plot(pivs,5, main = paste("p-value",vi,wb))
plot(pivs,6, main = paste("standard error",vi,wb))
plot(pivs,7, main = paste("n",vi,wb))
dev.off()
