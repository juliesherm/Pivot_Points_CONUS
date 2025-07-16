library(here)
library(terra)
library(dplyr)

#dt_pivs <- sf::st_read(dsn = here("results", "NCPN", "NCPN_Alliance_PivotPts_GDB", "NCPN_Alliance_PivotPts.gdb"), layer = "tbl_ClimateVariableData")

##plot GLAC steps
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
} else {
  park_boundary <- terra::vect(here("raw_data", "NPS_IMD_Units", "IMD_BND_NPSUnitAOAs.gdb"),
                               layer = "CONUS_Park_AOAs"
  )
  park_boundary <- subset(park_boundary, park_boundary$UNITCODE == park) %>%
    terra::project(., crs(raw_vi[[1]]))
}

glac_extent<-terra::ext(park_boundary)

## read in data for several key "steps" in the pivot-point calculation pipeline
vi = "EVI"
wb = "aet"


raw_vi<-terra::rast(here("raw_data","MODIS","CONUS","VI_16Days_250m_v61",vi,"tiled",
                         paste0(vi,"_2000_0_0.tif")),
                    win = glac_extent)
gs_dates<-terra::rast(here("results","CONUS","gs_dates","window_14days","cutoff_2.5","tiled","bilinear",
                           "gs_medians_0_0.tif"),
                      win = glac_extent)
int_vi_gs<-terra::rast(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","bilinear",vi,"tiled",
                            paste0(vi,"_int_0_0.tif")),
                       win = glac_extent)
raw_wb<-terra::rast(here("results","CONUS","wb_proj","wateryear","bilinear","tiled",
                         paste0("V_1_5_wateryears_gridmet_historical_",wb,"_resampled_0_0.tif")),
                    win = glac_extent)
pivs<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear","tiled",
                       paste0("pivs_",vi,"_",wb,"_0_0.tif")),
                  win = glac_extent)



pdf(here("results","GLAC",paste0(vi,wb,"_bilinear_plots.pdf")))
par(mfrow = c(5,2))
plot(raw_vi,1,  main = paste0(vi, " Feb 18, 2000"))
plot(raw_vi,11,  main = paste0(vi, " July 27, 2000"))
plot(gs_dates,1 , main = "growing season start")
plot(gs_dates,2 , main = "growing season end")
plot(int_vi_gs,1, main = paste("Integrated growing season", vi, "2000"))
plot(int_vi_gs,11, main = paste("Integrated growing season",vi,"2010"))
plot(raw_wb,1, main = paste(wb,"2000"))
plot(raw_wb,11, main = paste(wb,"2010"))
plot(pivs,3, main = paste("pivot",vi,wb))
plot(pivs,4, main = paste("R2",vi,wb))
plot(pivs,2, main = paste("slope",vi,wb))
plot(pivs,5, main = paste("p-value",vi,wb))
plot(pivs,7, main = paste("n",vi,wb))
dev.off()
