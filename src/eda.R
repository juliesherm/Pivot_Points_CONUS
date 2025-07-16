# Exploratory Data Analysis
setwd("C:/Users/Julie/Documents/MATH/NPS SIP/pivot_pts")
#setwd("/home/1024/ma/sherman/NPS_SIP/pivot_pts")

library(here)
library(terra)
library(dplyr)
library(ggpubr)
library(tidyterra)
library(scales)


#view histogram of integrated VI
vi_gs<-terra::rast(here("results","CONUS","integrated_vi_gs","window_14days","cutoff_2.5","bilinear","EVI","tiled","EVI_int_0_5222.tif"))
vi_gs


wb_wy<-terra::rast(here("results","CONUS","wb_proj","wateryear","bilinear","tiled","V_1_5_wateryears_gridmet_historical_accumswe_resampled_0_5222.tif"))
wb_wy<-terra::rast(here("results","CONUS","wb_proj","wateryear","bilinear","tiled","V_1_5_wateryears_gridmet_historical_soilwater_resampled_0_5222.tif"))
#choose a cell to look at that is in one of the polygons
set.seed(32233)
######fix best!!
cell_choices<-cells(vi_gs)
rand_cell<-sample(cell_choices,1)
rw<-rowFromCell(vi_gs, rand_cell)
cl<-colFromCell(vi_gs, rand_cell)
cell_vi<-unlist(vi_gs[rw,cl,]/10000)
cell_wb<-unlist(wb_wy[rw,cl,]/10)
hist(cell_vi)
hist(log(cell_vi))
plot(cell_vi~cell_wb)
plot(log(cell_vi)~log(cell_wb))
summary(lm(log(cell_vi)~log(cell_wb)))
summary(lm(cell_vi~cell_wb))

wb_wy
#plot(vi_gs)