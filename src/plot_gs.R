# Plot the growing season for CONUS or a specific park

setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")
library(here)
library(terra)
library(tidyterra)

#################################
# for CONUS
gs_first<-terra::rast(here("results","CONUS","gs_dates","window_14days","cutoff_2.5","bilinear","gs_first.tif"))
gs_last<-terra::rast(here("results","CONUS","gs_dates","window_14days","cutoff_2.5","bilinear","gs_last.tif"))

states<-terra::vect(here("raw_data","states","s_18mr25.shp"))
states<-crop(states,  gs_last)

#mask out non-conus
EVI_CONUS<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                            "EVI_deficit_slope.tif"))

gs_first<-mask(gs_first, EVI_CONUS,  maskvalues = NA)
gs_last<-mask(gs_last, EVI_CONUS,  maskvalues = NA)
gs_len<-gs_last-gs_first

#### Plot
#plot(gs_first, range = c(0,365), decreasing = TRUE,main = "Median start of growing season", 
#     axes = FALSE)
#plot(gs_last, range = c(0,365),main = "Median end of growing season",
#     axes = FALSE)
##plot(gs_first, range = c(365,0),main = "Median end of growing season", axis = FALSE)
#plot(gs_len,  main = "Growing Season Length",  range =  c(0,365), 
#     axes = FALSE)
ggplot()+
  geom_spatraster(data = gs_len)+
  scale_fill_distiller(palette = "YlGn",direction = 1,
                      limit = c(0,365), name = "Days",
                      na.value = "white")+
  geom_spatvector(data = states, col = "black", 
                  show.legend = FALSE,fill = NA, lwd = 0.5)+
  theme_void()+labs(subtitle = "Growing Season Length")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

###################################################
# For BLCA

source(here("src","get_ext.R"))

park <- "BLCA"

#### Read in park boundary
park_extent<-get_ext(park)

#get growing season
gs_first<-terra::rast(here("results","CONUS","gs_dates","window_14days","cutoff_2.5","bilinear","gs_first.tif"),
                      win = park_extent)
gs_last<-terra::rast(here("results","CONUS","gs_dates","window_14days","cutoff_2.5","bilinear","gs_last.tif"),
                     win = park_extent)
gs_len<-gs_last-gs_first

#plot for BLCA
plist<-list()
#first
plist[[1]]<-ggplot()+ geom_spatraster(data = gs_first, aes(fill = gs_first))+
  geom_spatvector(data = park_boundary, fill = NA, col = "black",  lwd = 1)+
  scale_fill_gradient(low="forestgreen", high = "white",
                       limit = c(77,116), name = "day of year") +
  theme_void()+labs(subtitle = "Median Start of Growing Season")
#last
plist[[2]]<-ggplot()+ geom_spatraster(data = gs_last, aes(fill = gs_last))+
  geom_spatvector(data = park_boundary, fill = NA, col = "black",  lwd = 1)+
  scale_fill_gradient(high="forestgreen",low = "white",
                       limit = c(277,296), name = "day of year") +
  theme_void()+labs(subtitle = "Median End of Growing Season")
#length
plist[[3]]<-ggplot()+ geom_spatraster(data = gs_len, aes(fill = gs_last))+
  geom_spatvector(data = park_boundary, fill = NA, col = "darkgrey",  lwd = 1)+
  scale_fill_gradient(low= "white", high="forestgreen",
                       limit = c(169,217), name = "# days") +
  theme_void()+labs(subtitle = "Growing Season Length")+
  geom_spatvector(data = value_polys, aes(col = ValueRate),show.legend = FALSE, fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))+
  geom_spatvector_text(aes(label = poly_num),data = value_polys,
                       fontface = "bold",
                       color = "black")
  

print(plist[[3]])

ggarrange(plist[[1]],plist[[2]],plist[[3]], nrow = 1, ncol = 3)

#################################################
#for COLM
plot(gs_first, range = c(64,85),main = "Median start of growing season",
     axes = FALSE)
lines(park_boundary)
plot(gs_last, range = c(250,310),main = "Median end of growing season",
     axes = FALSE)
lines(park_boundary)
plot(gs_len,  main = "Growing Season Length",  range =  c(180,235), 
     axes = FALSE)
lines(park_boundary)

