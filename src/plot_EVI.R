# Plot the growing season for CONUS or a specific park

setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")
library(here)
library(terra)
library(tidyterra)
library(scales)
source(here("src","get_ext.R"))


################################## 
park = "CONUS" # park can be a four-letter abreviation or "CONUS"

#### Read in park boundary
if(park=="CONUS"){
  park_extent = NULL
}else{
  park_extent<-get_ext(park)
  park_boundary<-get_ext(park, boundary = TRUE)
  #determine  which tile  the  park is in
  tile<-get_ext(park_extent)[[1]]
  xoff = tile[1]
  yoff = tile[2]
}

#### Read in mean iEVI
mean_EVI<-terra::rast(here("results","CONUS","integrated_vi_gs","mean_vi","meanEVI.tif"),
                      win = park_extent)
mean_EVI <-mean_EVI/10000
#### Read in states
states<-terra::vect(here("raw_data","states","s_18mr25.shp"))
states<-crop(states,  mean_EVI)

#mask out non-conus
EVI_CONUS<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                            "EVI_deficit_slope.tif"))

mean_EVI<-mask(mean_EVI, EVI_CONUS,  maskvalues = NA)

#ggplot
#set up plotting things
lims<-switch(park,
             BLCA = c(10,116),
             CONUS = c(0,200))
brks<-switch(park,
             BLCA = c(1, 2,3),
             CONUS = seq(lims[1], lims[2], by = 50))
num_brks<-length(brks)
brks_labs<-c(paste0("<", brks[1]), brks[2:(num_brks-1)], paste0(">", brks[num_brks]))
plist<-list()

#plot
plist[[1]]<-ggplot()+ geom_spatraster(data = mean_EVI, aes(fill = meanEVI))+
  scale_fill_distiller(palette="Greens",direction = "1",
                      limit = lims, name = NULL,
                        na.value = "white",oob = squish,
                        breaks =brks, labels = brks_labs)+
  theme_void()+labs(subtitle = "Historical Mean iEVI")+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

if(park =="BLCA"){
  plist[[1]]<-plit[[1]]+
      scale_fill_distiller(high="forestgreen", low = "white",
                          limit = lims, name = NULL,
                          na.value = "white") +
    geom_spatvector(data = value_polys, aes(col = ValueRate),show.legend = FALSE, 
                    fill = NA,lwd = 1) + 
    scale_color_manual(values=c("orange", "red", "#56B4E9"))+
    geom_spatvector(data = park_boundary, fill = NA, col = "darkgray",  lwd = 1)
  print(plist[[1]])
}else if(park=="CONUS"){
  plist[[1]]<-plist[[1]]+geom_spatvector(data = states, col = "black", 
                                         show.legend = FALSE,
                                         fill = NA, lwd = 0.5)
  print(plist[[1]])
}else{
  print(plist[[1]])
}
#


###########################################################
## plot the variance

all_EVI<-terra::rast(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","bilinear","EVI","tiled",
                          paste0("EVI_int_", xoff, "_",yoff,".tif")),
                      win = park_extent)
all_EVI <-all_EVI/10000
sd_EVI<-stdev(all_EVI)

#ggplot
plist<-list()
plist[[1]]<-ggplot()+ geom_spatraster(data = sd_EVI, aes(fill = std))+
  geom_spatvector(data = park_boundary, fill = NA, col = "black",  lwd = 1)+
  scale_fill_gradient(high="purple4", low = "white",
                      limit = c(0.5,22.6), name = NULL) +
  theme_void()+labs(subtitle = "Standard Deviation of Growing Season EVI")+
  geom_spatvector(data = value_polys, aes(col = ValueRate),show.legend = FALSE, fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))

print(plist[[1]])
