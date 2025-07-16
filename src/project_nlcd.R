setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")
library(here)
library(terra)
library(tidyterra)
library(ggplot2)
source(here("src","get_ext.R"))



#############################################3
# For all of CONUS
land_class<-terra::rast(here("raw_data","nlcd_2021","nlcd_2021_land_cover_l48_20230630.img"))
plot(land_class)

#read in or reproject projected land cover
project = FALSE
if(project){
  #project to MODIS  pixels
  conus_pix<-terra::rast(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","bilinear",
                              "EVI_int_CONUS.tif"))
  land_class<-terra::project(land_class,conus_pix,  method = "mode",
                             filename = here("raw_data","nlcd_2021","nlcd_2021_proj","land_cover_2021_proj.tif"),
                             overwrite =  TRUE)
}else{
  land_class_proj<-as.factor(terra::rast(here("raw_data","nlcd_2021","nlcd_2021_proj",
                                              "land_cover_2021_proj.tif")))
  coltab(land_class_proj)<-coltab(land_class)
#  land_class_proj<-mask(land_class_proj,  land_class_proj==0,
#                        maskvalues=TRUE, updatevalue= NA)
}


states<-terra::vect(here("raw_data","states","s_18mr25.shp"))
states<-crop(states,  land_class)

plot(land_class_proj,  axes = FALSE)
plot(states, add = TRUE)
#plot(land_class,  axes = FALSE,  mar = c(3.1,1.1,2.1,12.1))
par(mfrow = c(1,1))

#filter out developed land, open water, ...
veg_mask<-land_class_proj$`NLCD Land Cover Class`>40 & 
  (land_class_proj$`NLCD Land Cover Class`<80 | land_class_proj$`NLCD Land Cover Class`>89 ) 

####################################
#  Create and save a tif that has the  "cultivated land" etc that gets  masked
#cult_mask<-land_class_proj$`NLCD Land Cover Class`>11 & 
#  (land_class_proj$`NLCD Land Cover Class`<40 | land_class_proj$`NLCD Land Cover Class`>80 ) &
#  land_class_proj$`NLCD Land Cover Class`<89

#plot(cult_mask)
#writeRaster(cult_mask, filename=here("results","CONUS","nlcd_masked.tif"))

#######################################
nlcd_labels = c("Water", "Barren","Deciduous","Evergreen",
           "Shrub/Scrub","Grassland","Hay/Pasture",
           "Cultivated")
land_class_veg<-mask(land_class_proj, veg_mask, 
                     maskvalues=FALSE, 
                     updatevalue=NA,overwrite = TRUE,
                     filename=here("raw_data","nlcd_2021","nlcd_2021_proj","veg_cover_more_2021_proj.tif"))
  
  tryCatch(mask(land_class_proj, veg_mask, 
                     maskvalues=FALSE, 
     updatevalue=NA,overwrite = TRUE,
     filename=here("raw_data","nlcd_2021","nlcd_2021_proj","veg_cover_more_2021_proj.tif")),
     error = function(e) rast(here("raw_data","nlcd_2021","nlcd_2021_proj","veg_cover_more_2021_proj.tif")))
plot(land_class_veg)

# shading looks bad - the brown like fades out and looks pink...
#ggplot()+geom_spatraster(data = land_class_proj, alpha = 0.5)+
#  geom_spatraster(data = land_class_veg, alpha = 1)+
#  #scale_fill_coltab(data = land_class_proj,name = NULL)+ 
#  guides(alpha = "none", fill = "none") +
#  theme_void()+labs(subtitle = "NLCD Land Classification")

#############################################3
#just BLCA

park = "BLCA" # park can be a four-letter abreviation or "CONUS"

#### Read in park boundary
if(park=="CONUS"){
  park_extent = NULL
}else{
  park_extent<-get_ext(park)
  park_boundary<-get_ext(park, boundary = TRUE)
}

#### Read in value polygons (from park managers)
value_polys <- terra::vect(here("raw_data", "NCPN_PJ_Values_2024.gdb"))
value_polys<- subset(value_polys, value_polys$Park==park)
value_polys$poly_num<-1:9

#project to MODIS crs
raw_vi<-lapply(list.files(here("raw_data","MODIS","CONUS","VI_16Days_250m_v61","EVI","tiled"),
                          pattern = paste0("EVI","_2000_*"), full.names = T),
               terra::rast)
value_polys <- terra::project(value_polys,crs(raw_vi[[1]]))

#### Read in NPS vegetation map
nps_class<-terra::vect(here("raw_data","vegmaps","blcageodata",
                            "blcageodata.gdb"),
                       layer = "fcl_Veg_Polys")

nps_class <- terra::project(nps_class,crs(raw_vi[[1]]))
rm(raw_vi)
nps_class<-terra::intersect(nps_class, park_extent)



#get the top 10-ish (in terms of area) veg classes
ar_classes<-nps_class %>% group_by(Map_Unit_ID)%>% 
  summarize(Ar = sum(Hectares))%>%arrange(Ar)
top_classes<-subset(nps_class,nps_class$Map_Unit_ID %in%ar_classes[-c(1:20)]$Map_Unit_ID)

top_classes$Map_Unit_Common_Name<-substr(top_classes$Map_Unit_Common_Name,
                                         1,25)


# plotting!!

# plot all ~major~ types
ggplot()+
  geom_spatvector(aes(fill = Map_Unit_Common_Name), data = top_classes)+
  theme_void()+
  theme(
    legend.text = element_text(size = 10,vjust = 0.1), 
    legend.position = "right",
    legend.key.width = unit(0.2, "cm"),
    legend.key.size = unit(0.2, "cm"),
    legend.title = element_blank())+
  geom_spatvector(data = park_boundary, fill = NA, col = "darkgray",  lwd = 1)+
  geom_spatvector(data = value_polys, aes(col = ValueRate),show.legend = FALSE, fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))


# Julie-defined PJ types from NPS VEGMAP
#pj_classes<-subset(nps_class,  nps_class$Map_Unit_ID %in% c(69,71,72,74))
pj_classes<-subset(nps_class,  nps_class$Map_Unit_ID %in% c(69,71,72,74,78,87))

# plot just the pj types
ggplot()+
  geom_spatvector(aes(fill = Map_Unit_ID), data = pj_classes)+
  theme_bw()+  geom_spatvector(data = park_boundary, fill = NA, col = "darkgray",  lwd = 1)+
  geom_spatvector(data = value_polys, aes(col = ValueRate),show.legend = FALSE, fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))

##########
#NLCD
#### Read in NLCD land classification
land_class_proj<-terra::rast(here("raw_data","nlcd_2021","nlcd_2021_proj","land_cover_2021_proj.tif"),
                               win = park_extent)
  
land_class_proj<-as.factor(land_class_proj)
coltab(land_class_proj)<-coltab(land_class)

pj_masked<-mask(land_class_proj, pj_classes)
names(pj_masked)<-"NLCD"

ggplot()+geom_spatraster(data = land_class_proj, alpha = 0.4)+
  geom_spatraster(data = pj_masked,alpha = 1)+
  scale_fill_coltab(data = land_class_proj,
                    labels = c("Water", "Barren","Deciduous","Evergreen",
                               "Shrub/Scrub","Grassland","Hay/Pasture",
                               "Cultivated"), name = NULL)+ 
  guides(alpha = "none") +
  geom_spatvector(data = park_boundary, fill = NA, col = "darkgray",  lwd = 1)+
  theme_void()+labs(subtitle = "NLCD Land Classification")+
  geom_spatvector(data = value_polys, aes(col = ValueRate),show.legend = FALSE, fill = NA,lwd = 1) + 
  scale_color_manual(values=c("orange", "red", "#56B4E9"))

