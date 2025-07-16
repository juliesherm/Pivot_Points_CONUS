
library(here)
library(terra)
library(stringr)
library(ggplot2)
source(here("src","read_in_data.R"))
source(here("src","best_lm.R"))
source(here("src","get_ext.R"))

park <- "BLCA"

#### Read in park extent
park_extent<-get_ext(park)
#### Read in park boundary
park_boundary<-get_ext(park, boundary = TRUE)
#### Read in park veg types
#veg_types<-read_in_data(park, source = "vegmap")



#### Read in rsquares and p-values, and find the best fits for each pixel
#piv_files<-list.files(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear"), 
#                      pattern = "_rsq",full.names =  TRUE)
#pivotpts<-lapply(piv_files,terra::rast, win = park_extent)
#best<- best_lm(pivotpts)
best<-terra::rast(here("results","CONUS","pivot_pts","window_14days","cutoff_2.5","bilinear",
                       "best_r2.tif"),win = park_extent)
plot(best, col = rep(RColorBrewer::brewer.pal(n = 11, name = "Set3"),2),
     legend = FALSE, main = "Predictor with largest R2")
legend("bottomright",
       cex=0.7,inset = c(0,0.25),xpd = TRUE,
       legend = c("accumswe", "aet", "deficit", "precip","pet",
                  "rain","runoff","soilwater","tavg","tmax","tmin"), 
       fill = RColorBrewer::brewer.pal(n = 11, name = "Set3"))#### Choose pivot point to look at and read in
hist(best)

#### Assess significance of the regressions
print(paste(global(is.na(best["p-value.slp"]), "mean")$mean*100,"% of pixels have no significant relationships"))
plot(best["p-value.slp"] < 0.1, main = "Significant relationships")

#### Which pairs gave the best regressions
hist(best["which.max"])
plot(best["which.max"])

#### Determine sensitivity as the steepness of slope
plot(abs(best["slope"]), main = "Sensitivity")


#### Read in value polygons (from park managers)
value_polys <- terra::vect(here("raw_data", "NCPN_PJ_Values_2024.gdb"))
plot(value_polys, "ValueRate", col = c("green","yellow","red"))
value_polys <- terra::project(value_polys,crs(pivotpts[[1]]))
hist(pivotpts, 2)
plot(value_polys, "ValueRate",add = TRUE, col = c("green","yellow","red"))

#### Create combined value-sensitivity score
#mean sensitivity per polygon
value_polys$best_model_poly<-terra::extract(best['which.max'], value_polys, "mode")$which.max
value_polys$avsensitivity<-terra::extract(best["slope"], value_polys, fun = "mean", na.rm = TRUE)$slope
value_polys$numcells<-terra::extract(best["slope"], value_polys, fun = "mean", na.rm = TRUE)$slope
value_polys$avrsq<-terra::extract(best["r.squared"], value_polys, fun = "mean", na.rm = TRUE)$slope
value_polys$ValueRate<-as.numeric(stringr::str_split_i(value_polys$ValueRate," ",1))

value_polys$score<-(value_polys$ValueRate)^1*(1-value_polys$sensitivity/(value_polys$ValueRate/2.5)^(1/10))*3

#### Plot in value-sensitivity plane
vulnerability <-rep(seq(0,1, length.out = 1000),3)
value<-rep(1:3, each = 1000)
ggplot(mapping = aes(vulnerability, value, fill = (value)^1*(1-vulnerability/(value/2.5)^(1/10))*3))+
  geom_tile()+theme_bw()+ylab("value score")+
  scale_x_continuous("vulnerability", breaks = c(0,0.5,1),labels =  c("least","","most"))+
  scale_fill_gradient(low = "white", high = "red",breaks = c(1,4.5,8),
                      labels = c("accept", "direct", "resist"))+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10,
                                title = ""))
