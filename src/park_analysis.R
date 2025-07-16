## ---------------------------
##
## Script name: park_analysis.R
##
## Purpose of script: Calculates pivot points and other analysis for specified park
##
## Author: Julie Sherman, Carolyn Lober
##
## Date Created: 2/8/2024
## Last edited: 2/9/2024
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

#### General TO DO -
## how much/what type of temporal smoothing or averaging is required
## what is the temporal extent we care about
## adjust file structure/data retrieval for AWS
## compare "pixel-scale" (this script) with "polygon-scale" (already published) for Arches, Capitol Reef

# load packages
library(terra)
library(tidyverse)
library(here)
library(stringr)
library(tidyterra) # for plotting

# load functions
source(here("src", "read_in_data.R"))
source(here("src", "extract_or_read.R"))

park <- "GLAC"
# should work with any park or I&M network standard 4 letter code; for I&M networks,
# the script will use a bounding box around all parks in the network 
# for entire CONUS, use "CONUS" or set to NULL

#########################################
# Read in pivot points
#########################################

# read in park boundary, polygons, and veg types
vegmap <- read_in_data(park, "vegmap")
vegmap_8.3 <- subset(vegmap, vegmap[["Hectares"]] > 8.3)

# read in pivot points
pivotpts <- read_in_data(park, "pivotpts")

# read in vi trends
vi_trends <- read_in_data(park, "vi_trends")

# read in burned areas
burned_areas <- terra::rast(here("raw_data", "burned_areas", "Fire_Summary_Rasters_GeoTiffs",
                                 "USGS_Wildfires_Most_Recent_Year_Burned_Raster.tif")) 
burned_areas <- terra::project(burned_areas, pivotpts[[1]]) 

# read in DEM
if (!file.exists(here("results", park, "poly_dem.csv"))) {
  dem <- read_in_data(park, "DEM")
  dem <- terra::mosaic(terra::sprc(dem))
  
  # calculate and add slope and aspect
  add(dem) <- terra::terrain(dem, "dem.slope") 
  add(dem) <- terra::terrain(dem[[1]], "dem.aspect")
  names(dem)[1] <- "dem.elevation"
} else {
  dem <- NA
}

# get vi and wb vars from pivotpts names
vi_vars <- names(pivotpts) %>%
  str_split_i(., "_", 1) %>%
  unique
wb_vars <- names(pivotpts) %>%
  str_split_i(., "_", 2) %>%
  unique

#########################################
# Calculate proportion of pivot points with significant (p < 0.1) relationships between WB and VI
#########################################

# count number of cells with p-value.slp < 0.1 (about 95%)
pivotpts_pval <- rast(pivotpts) 
pivotpts_pval <- subset(pivotpts_pval,str_detect(names(pivotpts_pval), "6")) # layer 6 is p-value.slp
global(app(pivotpts_pval, function(vec) any(vec < 0.1)), sum, na.rm = TRUE) / ncell(pivotpts[[1]])

# which variables have the highest number of significant relationships?
pivotpts_sig <- lapply(pivotpts, \(stk) { stk[["p-value.slp"]] < 0.1 })
pivotpts_nsig <- lapply(pivotpts_sig, \(lyr) { global(lyr, sum, na.rm = TRUE) }) %>% unlist
names(pivotpts)[order(pivotpts_nsig, decreasing = TRUE)] # which var pairs have the most significant pixels? 
pivotpts_psig <- (pivotpts_nsig / ncell(pivotpts[[1]]))[order(pivotpts_nsig, decreasing = TRUE)]

#########################################
# Extract dem information for vegetation polygons 
#########################################

poly_dem <- extract_or_read(park, vegmap_8.3, "dem", dem, overwrite = FALSE)

#########################################
# Extract burned areas information for vegetation polygons 
#########################################

poly_burned <- extract_or_read(park, vegmap_8.3, "burned_areas", burned_areas, overwrite = FALSE)

#########################################
# Extract pivot points for vegetation polygons 
#########################################

poly_pivotpts <- extract_or_read(park, vegmap_8.3, "pivotpts", pivotpts, overwrite = FALSE)

# add dem, burned areas info
poly_pivotpts <- merge(poly_pivotpts, poly_dem, by = "OBJECTID", sort = FALSE)
poly_pivotpts <- merge(poly_pivotpts, poly_burned, by = "OBJECTID", sort = FALSE)

#########################################
# Simplify vegmap names
#########################################

#unique(poly_pivotpts$NVC2_L6X) 
poly_pivotpts$NVC2_L6_simple <- str_replace(poly_pivotpts$NVC2_L6X, "Rocky Mountain", "RM")

#########################################
# Plot pivot points
#########################################

# pivot points vs vegetation type
vi_plot <- "NDVI"
wb_plot <- "accumswe" 

# faceted lm plots
poly_pivotpts %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::filter(vi_var == vi_plot, wb_var == wb_plot, `p-value.slp` < 0.1) %>% 
  ggplot() + 
  geom_segment(aes(x = x.min, 
                   y = x.min*slope + y.intercept, 
                   xend = x.max, 
                   yend = x.max*slope + y.intercept,
                   col = dem.elevation)) + 
  facet_wrap(~ NVC2_L6X) + 
  xlab(paste0("Wateryear ", wb_plot, " (mm)")) + 
  ylab(paste0("Growing season ", vi_plot, " anomaly")) +
  geom_hline(yintercept=0)

# rsq bar plots
poly_pivotpts %>%
  pivot_wider(names_from = "type", values_from = "value") %>% 
  dplyr::filter(vi_var == vi_plot, wb_var %in% c("accumswe", "aet", "deficit", "pet", "rain", "soilwater")) %>% 
  group_by(NVC2_L6X, wb_var) %>%
  summarize(rsq_mean = mean(r.squared), rsq_se = sd(r.squared)/sqrt(length(r.squared)), 
            corr_sign = mode(slope > 0)) %>%
  ggplot(aes(x = wb_var, y = rsq_mean)) + 
  geom_bar(aes(fill = corr_sign), stat = "identity") + 
  scale_fill_manual(values = c("blue","red"), labels=c("Positive","Negative"), name = "Slope") +
  geom_errorbar(aes(y = rsq_mean, ymin = rsq_mean - rsq_se, ymax = rsq_mean + rsq_se), width = 0.25) +
  facet_wrap(~ NVC2_L6X) +
  xlab(paste0("Water balance variable")) +
  ylab(expression("Mean "~r^2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# pivot points spatial
poly_pivotpts %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::filter(vi_var == vi_plot, wb_var == wb_plot, `p-value.slp` < 0.1) %>%
  terra::merge(vegmap_8.3, ., by = "OBJECTID") %>% 
  ggplot() + 
  geom_spatvector(aes(fill = x.intercept), color = NA) +
  scale_fill_continuous(name = "Pivot point (mm)")

poly_pivotpts %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::filter(vi_var == vi_plot) %>%
  ggplot(aes(y = r.squared, x = wb_var, col = USGS_Wildfires_Most_Recent_Year_Burned_Raster)) +
  geom_point(position = "jitter") +
  scale_color_continuous(na.value = "#aaaaaa55")

poly_pivotpts %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::filter(vi_var == vi_plot) %>%
  ggplot(aes(y = r.squared, x = USGS_Wildfires_Most_Recent_Year_Burned_Raster)) +
  geom_point() +
  facet_wrap(~wb_var) + 
  geom_smooth() 
 
#########################################
# Extract vi trend for vegetation polygons 
#########################################
# get mean pivot point within veg polygons (and mean)

poly_vi_trends <- extract_or_read(park, vegmap_8.3, "vi_trends", vi_trends, overwrite = FALSE)

# add dem info
poly_vi_trends <- merge(poly_vi_trends, poly_dem, by = "OBJECTID", sort = FALSE)

#########################################
# Plot vi trends
#########################################

vi_plot <- "EVI"
type <- "max"

poly_vi_trends %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::filter(vi_var == vi_plot, max_or_mean == type, `p-value` < 0.1) %>% 
  ggplot(aes(x = reorder(NVC2_L6X, slope, na.rm=TRUE), y = slope)) +
  geom_boxplot() + 
  xlab("Vegetation type") + 
  ylab(paste0("Trend in ", type, " annual ", vi_plot)) +
  geom_hline(yintercept = 0) + 
  coord_flip() + 
  theme_minimal() +
  theme(legend.position = "right")

# trend spatial
poly_vi_trends %>%
  pivot_wider(names_from = "type", values_from = "value") %>%
  dplyr::filter(vi_var == vi_plot, max_or_mean == type, `p-value` < 0.1) %>% 
  terra::merge(vegmap_8.3, ., by = "OBJECTID") %>%
  ggplot() + 
  geom_spatvector(aes(fill = slope), color = NA) + 
  scale_fill_gradient(low = "red", high = "green")

#########################################
# Extract vegtypes by rasterizing
#########################################

# rasterize vegmap
vrastname <- here("results", park, "vegrast.tif")

if (file.exists(vrastname)) {
  vegrast <- terra::rast(vrastname)
} else {
  print("No vegetation raster yet. Creating from veg polygons")
  vegrast <- terra::project(vegmap, pivotpts[[1]]) %>%
    rasterize(., pivotpts[[1]], field = "NCV2_L6X", fun = max)
  terra::writeRaster(vegrast, vrastname)
}

# find the average intercept, slope, rsq by vegetation type
lm_by_veg <- lapply(pivotpts, terra::zonal, z = vegrast, na.rm = TRUE)
all_by_veg <- bind_rows(lm_by_veg)
all_by_veg$wb_var <- rep(names(lm_by_veg), each = length(unique(all_by_veg[NVC2_L6X])))
# subset to only the highest rsquared
best_by_veg <- all_by_veg %>%
  group_by(NVC2_L6X) %>%
  slice_max(r.squared)
best_by_veg$wb_var

# intraclass correlation
# checks how similar the veg-class values are compared to overall similarity
# all seem to be less than 0.2 - indicating bad clustering
# veg types not much more similar to self at different location than another type
slps <- pivotpts$EVI_deficit[, , 2]
clusts <- vegrast[, , 1]
summary_aov <- summary(aov(slps ~ as.factor(clusts)))
summary_aov[[1]][1, 2] / sum(summary_aov[[1]][, 2])

#########################################
# Completeness assessment (lines 437+ in ndvi_3_1)
#########################################
#### TO DO - not sure how necessary this is rn
# Extract date info
# determine file names and convert to date that can be used in xts later
nms <- list.files(paste0(getwd(), "/VI_16Days_250m_v61", "/SAVI"),
                  pattern = "tif$",
                  full.names = FALSE
)
length(nms) # number of dates analyzed
yr_doy <- substr(nms, start = 14, stop = 21)
yr <- as.numeric(substr(yr_doy, start = 1, stop = 4))
doy <- as.numeric(substr(yr_doy, start = 6, stop = 9))
dates_df <- data.frame(dayyear = doy, year = yr)
dates_df
dates_df$origin <- as.Date(paste0(dates_df$year, "-01-01"), tz = "MST") - 1
head(df)
ndvi_dates <- as.Date(dates_df$dayyear, origin = dates_df$origin)
# make wb and ndvi time stamps both posix so they can be merged
mydate <- as.POSIXlt(ndvi_dates, format = "%d/%b/%Y:%H:%M:%S")
head(mydate)
ndvi_dates <- as.POSIXct(mydate, format = "%m/%d/%Y")
head(ndvi_dates) # intermediate step to drop time of day
mnths <- lubridate::month(ndvi_dates)

# annually
pct_complete_ann <- rowMeans(!is.na(SAVI_df[, -c(1, 2)]))

# growing season
gs_mos <- c(3:10)
non_gs_cols <- 2 + which(!(mnths %in% gs_mos)) # plus 2 from the x and y cols
pct_complete_gs <- rowMeans(!is.na(SAVI_df[, -c(1, 2, non_gs_cols)]))